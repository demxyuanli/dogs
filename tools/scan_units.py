#!/usr/bin/env python3
"""scan_units — rebuild migration-index.json from actual workspace code.

Scans crates/*/src for all .rs files (excluding lib.rs + target/), one file
per migration unit. Unit id = crate-prefixed snake name:
  - crates/occt-core/src/gp/xyz.rs  -> "gp_xyz"      (occt-core: no prefix, keeps legacy ids)
  - crates/occt-geom/src/curve.rs   -> "geom_curve"
  - crates/occt-geom2d/src/curve.rs -> "geom2d_curve"
  - crates/occt-math/src/bfgs.rs    -> "math_bfgs"
  - crates/occt-topo/src/abs.rs     -> "topo_abs"
Deps extracted from `use crate::` (intra-crate) and `use occt_*::` (cross-crate).
status = done (code exists).

Usage: python3 tools/scan_units.py [--write]
"""
import os
import re
import json
import sys
from collections import Counter
from datetime import datetime, timezone

WORKSPACE = os.path.join(os.path.dirname(__file__), "..", "crates")
OUT = os.path.join(os.path.dirname(__file__), "..", "migration-index.json")

# crate dir -> unit id prefix (empty = no prefix, legacy occt-core ids)
CRATE_PREFIX = {
    "occt-core": "",
    "occt-geom": "geom",
    "occt-geom2d": "geom2d",
    "occt-math": "math",
    "occt-topo": "topo",
}
# crate dir -> crate name in use statements
CRATE_NAME = {
    "occt-core": "occt_core",
    "occt-geom": "occt_geom",
    "occt-geom2d": "occt_geom2d",
    "occt-math": "occt_math",
    "occt-topo": "occt_topo",
}
# crate-name -> prefix, for cross-crate dep resolution
NAME_TO_PREFIX = {v: k for k, v in CRATE_PREFIX.items() if k}
NAME_TO_PREFIX.update({"occt_core": ""})
# crate-name -> crate dir
NAME_TO_CRATE = {v: k for k, v in CRATE_NAME.items()}

REGION = {
    "gp": "wide", "kernel": "repl", "precision": "repl", "bspl": "deep",
    "bvh": "deep", "bnd": "deep", "geom": "deep", "elib": "deep",
    "poly": "deep", "toploc": "deep", "int": "deep", "io": "deep",
    "convert": "deep", "cslib": "deep", "gprop": "deep", "quantity": "deep",
    "message": "periph", "numeric": "deep", "hull": "deep", "gcpnts": "deep",
    "geom": "deep", "geom2d": "deep", "math": "deep", "topo": "deep",
}


def now():
    return datetime.now(timezone.utc).isoformat()


def scan_files():
    """Return list of (crate_dir, abs_path, unit_id, mod, base, is_mod)."""
    files = []
    for crate in sorted(os.listdir(WORKSPACE)):
        src = os.path.join(WORKSPACE, crate, "src")
        if not os.path.isdir(src):
            continue
        prefix = CRATE_PREFIX.get(crate, crate.replace("occt-", ""))
        for dp, dn, fn in os.walk(src):
            if "target" in dp:
                continue
            for f in sorted(fn):
                if not f.endswith(".rs") or f == "lib.rs":
                    continue
                rel = os.path.relpath(dp, src)
                mod = rel.replace(os.sep, "/") if rel != "." else "."
                base = f[:-3]
                if base == "mod":
                    mid = "mod"
                elif mod == ".":
                    mid = base
                else:
                    mid = f"{mod}_{base}"
                uid = f"{prefix}_{mid}" if prefix else mid
                files.append((crate, os.path.join(dp, f), uid, mod, base, base == "mod"))
    return files


def first_symbol(path):
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


def extract_deps(path):
    """Extract `use X::Y::Z` and `use X::Y::{A,B}` tokens (crate:: or occt_*::).
    Returns set of (scope, "mod::[Symbol|file]")."""
    src = open(path, encoding="utf-8", errors="replace").read()
    deps = set()
    # use crate::a::b::{C, D} / use occt_x::a::{C,D}
    for m in re.finditer(r"use\s+(crate|occt_[a-z0-9_]+)::([A-Za-z0-9_:]+)::\{([A-Za-z0-9_]+)(?:\s*,\s*([A-Za-z0-9_]+))*\}", src):
        scope, path_part = m.group(1), m.group(2)
        for s in m.groups()[2:]:
            if s:
                deps.add((scope, f"{path_part}::{s}"))
    # use crate::a::b::C / use occt_x::a::C
    for m in re.finditer(r"use\s+(crate|occt_[a-z0-9_]+)::([A-Za-z0-9_:]+)", src):
        deps.add((m.group(1), m.group(2)))
    return deps


def main():
    write = "--write" in sys.argv
    files = scan_files()

    target_of = {}
    deps_of = {}
    syms_of = {}
    for crate, path, uid, mod, base, is_mod in files:
        target_of[uid] = first_symbol(path)
        deps_of[uid] = extract_deps(path)
        # pub symbol names (for module re-export resolution)
        src = open(path, encoding="utf-8", errors="replace").read()
        syms_of[uid] = set(re.findall(r"pub\s+(?:struct|enum|trait|type)\s+([A-Za-z_]\w*)", src))

    unit_ids = set(uid for _, _, uid, _, _, _ in files)

    # module symbol index: for intra-crate mod.rs re-exports and inline types
    # key: (crate_dir, "mod::Symbol") -> unit
    mod_symbol = {}
    for crate, path, uid, mod, base, is_mod in files:
        if not is_mod:
            continue
        prefix = CRATE_PREFIX.get(crate, crate.replace("occt-", ""))
        for sym in syms_of[uid]:
            key = (crate, f"{mod}::{sym}")
            mod_symbol[key] = uid
        src = open(path, encoding="utf-8", errors="replace").read()
        for pm in re.finditer(r"pub\s+use\s+(?:self::)?(\w+)::([A-Za-z_]\w*);", src):
            fname, sym = pm.group(1), pm.group(2)
            child = f"{prefix}_{mod}_{fname}" if prefix else f"{mod}_{fname}"
            if child in unit_ids:
                mod_symbol[(crate, f"{mod}::{sym}")] = child

    def resolve(crate, dep):
        """dep = (scope, rest). Returns unit id or None."""
        scope, rest = dep
        parts = rest.split("::")
        if scope == "crate":
            # intra-crate: same prefix as current crate
            prefix = CRATE_PREFIX.get(crate, crate.replace("occt-", ""))
            mod = parts[0] if parts else ""
            # try module::Symbol (mod.rs inline)
            if parts:
                key = (crate, f"{mod}::{parts[-1]}")
                if key in mod_symbol:
                    return mod_symbol[key]
            # try mod_file
            acc = mod
            for p in parts[1:] or [""]:
                cand = f"{prefix}_{acc}_{p}" if prefix else f"{acc}_{p}"
                if cand in unit_ids:
                    return cand
                if p:
                    acc = f"{acc}_{p}"
            # module fallback
            cand = f"{prefix}_{mod}_mod" if prefix else f"{mod}_mod"
            if cand in unit_ids:
                return cand
            # file at mod level (no submodule)
            cand = f"{prefix}_{mod}" if prefix else mod
            if cand in unit_ids:
                return cand
            return None
        else:
            # cross-crate: occt_geom::file::Sym -> geom_file
            prefix = NAME_TO_PREFIX.get(scope)
            crate_dir = NAME_TO_CRATE.get(scope)
            if prefix is None:
                return None
            if not parts:
                return None
            # try mod::Symbol via target crate's mod.rs re-export index
            if crate_dir and len(parts) >= 2:
                key = (crate_dir, f"{parts[0]}::{parts[-1]}")
                if key in mod_symbol:
                    return mod_symbol[key]
            # try geom_file
            cand = f"{prefix}_{parts[0]}"
            if cand in unit_ids:
                return cand
            # try geom_file_sym (nested) or geom_mod
            if len(parts) >= 2:
                cand2 = f"{prefix}_{parts[0]}_{parts[1]}"
                if cand2 in unit_ids:
                    return cand2
            cand3 = f"{prefix}_{parts[0]}_mod"
            if cand3 in unit_ids:
                return cand3
            return None

    now_ts = now()
    units = {}
    for crate, path, uid, mod, base, is_mod in files:
        deps = sorted({resolve(crate, d) for d in deps_of[uid] if resolve(crate, d)})
        deps = [d for d in deps if d != uid]
        top_mod = mod.split("/")[0] if mod != "." else "unit"
        # crate-level region: topo/math/geom are deep
        if top_mod == ".":
            top_mod = "unit"
        u = {
            "module": f"{crate}/{mod}".replace("/./", "/") if mod != "." else crate,
            "layer": 0,
            "status": "done",
            "region": REGION.get(top_mod, "unit"),
            "deps": deps,
            "target": target_of[uid],
            "target_path": os.path.relpath(path, os.path.join(WORKSPACE, "..")).replace(os.sep, "/"),
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
        print(f"  with target: {sum(1 for u in units.values() if u['target'])}/{len(units)}")
        print(f"  with deps: {sum(1 for u in units.values() if u['deps'])}/{len(units)}")
        c = Counter(u["module"].split("/")[0] for u in units.values())
        for m, n in sorted(c.items()):
            print(f"  {m:12s} {n} units")
        # broken deps
        bad = [(uid, dep) for uid, u in units.items() for dep in u["deps"] if dep not in units]
        print(f"  broken deps: {len(bad)}")
        for b in bad[:8]:
            print(f"    {b}")
        # samples
        for uid in ["gp_xyz", "geom_curve", "geom2d_curve", "math_bfgs", "topo_abs", "bvh_mod"]:
            if uid in units:
                print(f"  {uid}: target={units[uid]['target']} deps={units[uid]['deps'][:6]}")


if __name__ == "__main__":
    main()
