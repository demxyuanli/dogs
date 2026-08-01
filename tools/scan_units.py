#!/usr/bin/env python3
"""scan_units — rebuild migration-index.json from actual occt-core code.

Scans crates/occt-core/src for all .rs files (excluding lib.rs + target/),
treats each file as one migration unit:
  - src/gp/xyz.rs  -> unit "gp_xyz"  (struct GpXyz)
  - src/bvh/mod.rs -> unit "bvh_mod" (mod-level types like BvhNode)
  - src/numeric.rs -> unit "numeric" (top-level file)
status = done (code exists), deps extracted from `use crate::` refs,
target/target_path from the file's first pub struct/enum/trait.

Usage: python3 tools/scan_units.py [--write]
"""
import os
import re
import json
import sys
from datetime import datetime, timezone

CRATE_SRC = os.path.join(os.path.dirname(__file__), "..", "crates", "occt-core", "src")
OUT = os.path.join(os.path.dirname(__file__), "..", "migration-index.json")

# region by top-level module (assembly-line-ability class)
REGION = {
    "gp": "wide", "kernel": "repl", "precision": "repl", "bspl": "deep",
    "bvh": "deep", "bnd": "deep", "geom": "deep", "elib": "deep",
    "poly": "deep", "toploc": "deep", "int": "deep", "io": "deep",
    "convert": "deep", "cslib": "deep", "gprop": "deep", "quantity": "deep",
    "message": "periph", "numeric": "deep", "hull": "deep", "gcpnts": "deep",
}


def now():
    return datetime.now(timezone.utc).isoformat()


def scan_files():
    """Return list of (abs_path, unit_id, modname, basename, is_mod)."""
    files = []
    for dp, dn, fn in os.walk(CRATE_SRC):
        if "target" in dp:
            continue
        for f in sorted(fn):
            if not f.endswith(".rs") or f == "lib.rs":
                continue
            rel = os.path.relpath(dp, CRATE_SRC)
            mod = rel.replace(os.sep, "/") if rel != "." else "."
            base = f[:-3]
            if base == "mod":
                uid = f"{mod}_mod"
            elif mod == ".":
                uid = base
            else:
                uid = f"{mod}_{base}"
            files.append((os.path.join(dp, f), uid, mod, base, base == "mod"))
    return files


def extract_file(path):
    """Return (pub_symbols, use_crate_tokens) for a file."""
    src = open(path, encoding="utf-8", errors="replace").read()
    syms = set(re.findall(r"pub\s+(?:struct|enum|trait|type)\s+([A-Za-z_]\w*)", src))
    deps = set(re.findall(r"use\s+crate::([A-Za-z_]\w*(?:::[A-Za-z_]\w*)*)", src))
    return syms, deps


def first_symbol(path):
    """Best name for a file: struct/enum/trait > type alias > fn > const."""
    src = open(path, encoding="utf-8", errors="replace").read()
    for pat in (
        r"pub\s+(?:struct|enum|trait)\s+([A-Za-z_]\w*)",
        r"pub\s+type\s+([A-Za-z_]\w*)",
        r"pub\s+fn\s+([A-Za-z_]\w*)",
        r"pub\s+const\s+([A-Za-z_]\w*)",
    ):
        m = re.search(pat, src)
        if m:
            return m.group(1)
    return None


def main():
    write = "--write" in sys.argv
    files = scan_files()

    # per-file target (first pub symbol, prefer struct/enum)
    target_of = {}
    syms_of = {}      # unit -> {pub symbols}
    deps_of = {}      # unit -> set of crate:: tokens
    for path, uid, mod, base, is_mod in files:
        syms, deps = extract_file(path)
        syms_of[uid] = syms
        deps_of[uid] = deps
        target_of[uid] = first_symbol(path)

    # module symbol index: {mod}::{Symbol} -> unit
    # (mod.rs re-exports and inline types)
    mod_symbol = {}
    for path, uid, mod, base, is_mod in files:
        if not is_mod:
            continue
        for sym in syms_of[uid]:
            mod_symbol[f"{mod}::{sym}"] = uid
        # mod.rs may also re-export child symbols: read its pub use lines
        src = open(path, encoding="utf-8", errors="replace").read()
        for pm in re.finditer(r"pub\s+use\s+(?:self::)?(\w+)::([A-Za-z_]\w*);", src):
            fname = pm.group(1)
            sym = pm.group(2)
            child = f"{mod}_{fname}"
            if child in syms_of:
                mod_symbol[f"{mod}::{sym}"] = child

    # filename index for resolution: {mod}_{name} or {name}
    unit_ids = set(uid for _, uid, _, _, _ in files)

    def resolve(dep_token):
        """crate::gp::xyz::GpXyz -> unit id. Returns None if unresolvable."""
        parts = dep_token.split("::")
        mod = parts[0]
        rest = parts[1:]
        # exact module::Symbol (mod.rs inline or re-export)
        if rest:
            key = f"{mod}::{rest[-1]}"
            if key in mod_symbol:
                return mod_symbol[key]
        # walk path: crate::a::b::c  -> try a_b, a_b_c
        if mod != ".":
            acc = mod
            for p in rest:
                cand = f"{acc}_{p}"
                if cand in unit_ids:
                    return cand
                acc = cand
        else:
            for p in rest:
                if p in unit_ids:
                    return p
        # module-level fallback: dep on the module's own mod.rs
        if f"{mod}_mod" in unit_ids:
            return f"{mod}_mod"
        return None

    # build units
    now_ts = now()
    units = {}
    for path, uid, mod, base, is_mod in files:
        deps = sorted({resolve(d) for d in deps_of[uid] if resolve(d)})
        deps = [d for d in deps if d != uid]  # drop self-refs
        u = {
            "module": mod,
            "layer": 0,
            "status": "done",
            "region": REGION.get(mod.split("/")[0], "unit"),
            "deps": deps,
            "target": target_of[uid],
            "target_path": os.path.relpath(path, os.path.dirname(CRATE_SRC)).replace(os.sep, "/"),
            "updated_at": now_ts,
        }
        units[uid] = u

    idx = {"version": 1, "units": units}

    if write:
        os.makedirs(os.path.dirname(OUT), exist_ok=True)
        with open(OUT, "w", encoding="utf-8") as f:
            json.dump(idx, f, indent=2, ensure_ascii=False)
        print(f"Wrote {len(units)} units to {OUT}")
    else:
        print(f"Dry-run: {len(units)} units")
        n_tgt = sum(1 for u in units.values() if u["target"])
        n_dep = sum(1 for u in units.values() if u["deps"])
        print(f"  with target struct: {n_tgt}/{len(units)}")
        print(f"  with deps: {n_dep}/{len(units)}")
        # per-module counts
        from collections import Counter
        c = Counter(u["module"].split("/")[0] for u in units.values())
        for m, n in sorted(c.items()):
            print(f"  {m:12s} {n} units")
        # sample
        for uid in ["gp_xyz", "bvh_mod", "bspl_bezier", "bnd_b2b3", "geom_delaunay", "kernel_handle"]:
            if uid in units:
                print(f"\n  {uid}: target={units[uid]['target']} deps={units[uid]['deps']}")


if __name__ == "__main__":
    main()
