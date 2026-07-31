"""Reusable differential core — domain-agnostic.

Oracle strategies + per-op dispatch + the diff loop + build/run helpers. No
OCCT, no floating-point assumption, no domain paths baked in. A new domain
supplies only: a C/C++ driver, a Rust runner, an ops declaration, and a corpus.
Both the geometry and number-theory harnesses import THIS unchanged —
the proof the core is not geometry-specific.
"""

import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


# --- oracle strategies (pluggable; op picks one via its ops declaration) ------

def approx_eq(x, y, p):
    if x == y:
        return True
    d = abs(x - y)
    return d <= p["abs"] or d <= p["rel"] * max(abs(x), abs(y))


def struct_eq(g, c, p):
    if isinstance(g, dict) and isinstance(c, dict):
        return g.keys() == c.keys() and all(struct_eq(g[k], c[k], p) for k in g)
    if isinstance(g, list) and isinstance(c, list):
        return len(g) == len(c) and all(struct_eq(a, b, p) for a, b in zip(g, c))
    if isinstance(g, (int, float)) and isinstance(c, (int, float)):
        return approx_eq(float(g), float(c), p)
    return g == c


def s_bit_equal(g, c, _p):
    """Byte-for-byte identical (deterministic: hashing, parsers, integer math)."""
    return g.strip() == c.strip()


def s_numeric_tol(g, c, p):
    """Flat number list within absolute OR relative tolerance."""
    gv, cv = [float(x) for x in g.split()], [float(x) for x in c.split()]
    return len(gv) == len(cv) and all(approx_eq(a, b, p) for a, b in zip(gv, cv))


def s_struct_tol(g, c, p):
    """JSON object, same shape, numeric leaves within tolerance."""
    return struct_eq(json.loads(g), json.loads(c), p)


STRATEGIES = {
    "bit_equal": s_bit_equal,
    "numeric_tol": s_numeric_tol,
    "struct_tol": s_struct_tol,
}

_DEFAULT_ABS_REL = {"abs": 1e-9, "rel": 1e-12}  # hoisted constant (#10)


def op_policy(ops, op):
    """(strategy_fn, params) declared for op in the ops dict, or None.

    Strategy name is validated against STRATEGIES — a typo is caught at
    declaration time, not silently during diff.  NameError contains the
    offending op name for fast diagnosis (#6)."""
    d = ops.get(op)
    if not isinstance(d, dict) or "oracle" not in d:
        return None
    oracle = d["oracle"]
    if oracle not in STRATEGIES:                                          # (#6)
        raise NameError(f"op {op!r}: unknown oracle strategy {oracle!r}") # (#6)
    dflt = ops.get("_defaults", _DEFAULT_ABS_REL)
    if not isinstance(dflt, dict):
        raise TypeError(f"op {op!r}: _defaults must be a dict, got {type(dflt).__name__}")
    params = {"abs": d.get("abs", dflt["abs"]), "rel": d.get("rel", dflt["rel"])}
    return STRATEGIES[oracle], params


def diff(inputs, golden, cand, ops):
    """inputs: corpus lines; golden/cand: two runners' output lines. Returns
    (mismatches, ok_matched, err_matched, per_op). Classification (ok/err/noop)
    is universal; payload comparison uses the op's declared strategy. per_op maps
    op-id -> {'ran','mismatch'} so verify results flow back to the index.

    Undeclared ops are logged as mismatches (not a fatal exit) so write_results()
    still runs and reconcile gets a complete picture (#5)."""
    if not (len(inputs) == len(golden) == len(cand)):
        sys.exit(f"line count mismatch: in={len(inputs)} a={len(golden)} b={len(cand)}")
    mismatches, oks, errs = [], 0, 0
    per_op = {}

    # Precompute strategy for every declared op — one lookup per op, not per
    # input line (#10)
    policy_cache = {}
    for opin in ops:
        if not opin.startswith("_"):
            try:
                policy_cache[opin] = op_policy(ops, opin)
            except NameError:
                policy_cache[opin] = None   # will be reported as undeclared later

    for i, (line, g, c) in enumerate(zip(inputs, golden, cand)):
        parts = line.split(maxsplit=1)                          # (#4, #14 — single split)
        op = parts[0].strip() if parts else ""                  # (#4 — strip leading whitespace)
        if op not in per_op:                                    # (#11 — no setdefault alloc)
            per_op[op] = {"ran": 0, "mismatch": 0}
        per_op[op]["ran"] += 1

        gh, _, gp = g.partition(" ")
        ch, _, cp = c.partition(" ")
        if gh != ch:
            mismatches.append((i, line, g, c))
            per_op[op]["mismatch"] += 1
            continue
        if gh != "ok":
            errs += 1
            continue

        pol = policy_cache.get(op)  # None for undeclared ops (not pre-cached at L94)
        if pol is None:
            # Undeclared op — log as mismatch, don't kill the harness (#5)
            per_op[op]["mismatch"] += 1
            mismatches.append((i, line, g, c))
            continue
        strat, params = pol
        if strat(gp, cp, params):
            oks += 1
        else:
            mismatches.append((i, line, g, c))
            per_op[op]["mismatch"] += 1
    return mismatches, oks, errs, per_op


def write_results(path, per_op):
    """Persist per-op pass/fail for `migindex reconcile`.
    Appends a run record to the `runs` array so verify history accumulates.
    Backward-compatible: migrates old flat {op: {ran, pass}} to {runs: [...]} on first write."""
    entry = {
        "at": datetime.now(timezone.utc).isoformat(),
        "ops": {op: {"ran": r["ran"], "pass": r["mismatch"] == 0} for op, r in per_op.items()},
    }
    p = Path(path)
    if p.exists():
        try:
            data = json.loads(p.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            data = {"runs": []}  # treat corrupted file like missing — start fresh
        if "runs" in data:
            data["runs"].append(entry)
        else:
            # Migrate old flat format: wrap existing data as run 0
            data = {"runs": [{"at": None, "ops": data}, entry]}
    else:
        data = {"runs": [entry]}
    p.write_text(json.dumps(data, indent=2), encoding="utf-8")


# --- build/run helpers (generic: cargo + subprocess over the line protocol) ---

_EXE_SUFFIX = os.environ.get("EXE_SUFFIX", ".exe" if sys.platform == "win32" else "")


def build(manifest_dir, timeout=300):
    manifest_dir = Path(manifest_dir)
    subprocess.run(
        ["cargo", "build", "--quiet", "--manifest-path",
         str(manifest_dir / "Cargo.toml")],
        check=True,
        timeout=timeout,                             # (#7)
    )
    return manifest_dir / "target/debug" / (manifest_dir.name + _EXE_SUFFIX)  # (#16)


def run(exe, corpus, extra_path_dirs=(), timeout=120):
    env = dict(os.environ)
    if extra_path_dirs:
        env["PATH"] = os.pathsep.join(str(d) for d in extra_path_dirs) + os.pathsep + env["PATH"]
    p = subprocess.run(
        [str(exe)], input=corpus, capture_output=True, text=True,
        env=env, timeout=timeout,                    # (#7)
    )
    if p.returncode != 0:
        sys.exit(f"{Path(exe).name} failed (rc={p.returncode}):\n{p.stderr}")
    return p.stdout.splitlines()


# --- selftest ----------------------------------------------------------------

def selftest():
    import tempfile

    # -- approx_eq
    assert approx_eq(1.0, 1.0, {"abs": 1e-9, "rel": 1e-12})
    assert approx_eq(1.0, 1.0 + 1e-10, {"abs": 1e-9, "rel": 1e-12})
    assert approx_eq(1000.0, 1000.0 + 1e-10, {"abs": 1e-9, "rel": 1e-12})  # rel: 1e-10 <= 1e-12*1e3=1e-9
    assert not approx_eq(1.0, 1.0 + 2e-9, {"abs": 1e-9, "rel": 1e-15})    # outside both

    # -- bit_equal (3rd arg ignored, unified strategy interface)
    assert s_bit_equal("42  \n", "42\n", {})
    assert not s_bit_equal("42", "43", {})

    # -- numeric_tol
    p = {"abs": 1e-9, "rel": 1e-12}
    assert s_numeric_tol("1 2 3", "1 2 3", p)
    assert not s_numeric_tol("1 2 3", "1 2 4", p)
    assert not s_numeric_tol("1 2", "1 2 3", p)  # length mismatch

    # -- struct_tol
    assert s_struct_tol('{"a": 1.0, "b": [2.0]}', '{"a": 1.0, "b": [2.0]}', p)
    assert not s_struct_tol('{"a": 1.0}', '{"a": 2.0}', p)
    assert not s_struct_tol('{"a": 1.0}', '{"b": 1.0}', p)           # key mismatch
    assert not s_struct_tol('{"a": [1.0]}', '{"a": [1.0, 2.0]}', p) # list length

    # -- op_policy
    ops = {"echo": {"oracle": "bit_equal"}, "calc": {"oracle": "numeric_tol"}}
    strat, params = op_policy(ops, "echo")
    assert strat is STRATEGIES["bit_equal"]
    assert strat is not None
    assert op_policy(ops, "nonexistent") is None
    try:
        op_policy({"bad": {"oracle": "typoed_strategy"}}, "bad")
        assert False, "should have raised NameError"
    except NameError:
        pass

    # -- diff: all ok
    inputs = ["echo 42", "calc 1 2"]
    golden = ["ok 42", "ok 1 2"]
    cand = ["ok 42", "ok 1 2"]
    mismatches, oks, errs, per_op = diff(inputs, golden, cand, ops)
    assert mismatches == []
    assert oks == 2
    assert errs == 0
    assert per_op["echo"]["ran"] == 1 and per_op["echo"]["mismatch"] == 0
    assert per_op["calc"]["ran"] == 1 and per_op["calc"]["mismatch"] == 0

    # -- diff: error classification match
    inputs = ["echo 42", "calc 1 2"]
    golden = ["err", "err"]
    cand = ["err", "err"]
    mismatches, oks, errs, _ = diff(inputs, golden, cand, ops)
    assert oks == 0
    assert errs == 2  # both classified as err on both sides
    assert mismatches == []

    # -- diff: undeclared op → mismatch
    inputs = ["unknown_op 1"]
    golden = ["ok 1"]
    cand = ["ok 1"]
    mismatches, oks, errs, _ = diff(inputs, golden, cand, ops)
    assert len(mismatches) == 1

    # -- diff: mismatched payload
    inputs = ["calc 1 2"]
    golden = ["ok 3.0"]
    cand = ["ok 4.0"]
    mismatches, _, _, _ = diff(inputs, golden, cand, ops)
    assert len(mismatches) == 1

    # -- diff: ok vs err classification mismatch
    inputs = ["calc 1 2"]
    golden = ["ok 3.0"]
    cand = ["err"]
    mismatches, _, _, _ = diff(inputs, golden, cand, ops)
    assert len(mismatches) == 1

    # -- write_results
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "results.json"

        # first write
        write_results(path, {"echo": {"ran": 3, "mismatch": 0}})
        data = json.loads(path.read_text(encoding="utf-8"))
        assert "runs" in data and len(data["runs"]) == 1
        assert data["runs"][0]["ops"]["echo"]["pass"] is True

        # append second run
        write_results(path, {"echo": {"ran": 3, "mismatch": 1}})
        data = json.loads(path.read_text(encoding="utf-8"))
        assert len(data["runs"]) == 2
        assert data["runs"][1]["ops"]["echo"]["pass"] is False

        # backwards compat: old flat format gets migrated
        path2 = Path(d) / "old.json"
        path2.write_text('{"echo": {"ran": 1, "pass": true}}', encoding="utf-8")
        write_results(path2, {"echo": {"ran": 2, "mismatch": 0}})
        data2 = json.loads(path2.read_text(encoding="utf-8"))
        assert len(data2["runs"]) == 2 and data2["runs"][0]["at"] is None

    print("selftest ok")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "selftest":
        selftest()
