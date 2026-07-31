#!/usr/bin/env python3
"""OCCT CAD kernel harness — differential test via engine/differ.py."""
import json, math, os, random, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import differ

ROOT = Path(__file__).resolve().parent
OPS = json.loads((ROOT / "ops.json").read_text(encoding="utf-8"))

def _mig_root():
    v = os.environ.get("MIG_ROOT")
    if v: return Path(v)
    dotfile = ROOT / ".mig_root"
    if dotfile.exists(): return Path(dotfile.read_text(encoding="utf-8").strip())
    sys.exit("MIG_ROOT not set and .mig_root not found.")

OCCT = _mig_root()
DLL_DIRS = []

def vec(rng, scale=1.0):
    return [repr((rng.random() * 2 - 1) * scale) for _ in range(3)]

def mat9(rng):
    return " ".join(repr((rng.random()*2-1)*10) for _ in range(9))
def mat4(rng):
    return " ".join(repr((rng.random()*2-1)*10) for _ in range(4))
def vec2(rng, scale=1.0):
    return [repr((rng.random()*2-1)*scale) for _ in range(2)]

def make_corpus(n=2000):
    rng = random.Random(0xC0FFEE)
    lines = ["echo 42", "echo -1", "echo 0"]
    for _ in range(n):
        scale = rng.choice([1.0, 1e3, 1e6])
        a, b = vec(rng, scale), vec(rng, scale)
        lines.append("gp_XYZ::Dot " + " ".join(a + b))
        lines.append("gp_XYZ::Crossed " + " ".join(a + b))
        lines.append("gp_XYZ::Modulus " + " ".join(a))
        lines.append("gp_XYZ::Added " + " ".join(a + b))
        lines.append("gp_XYZ::MultipliedScalar " + " ".join(a + [repr(rng.random()*10-5)]))
        lines.append("gp_XYZ::CrossSquareMagnitude " + " ".join(a + b))
        c = vec(rng, scale)
        lines.append("gp_XYZ::DotCross " + " ".join(a + b + c))
        lines.append("gp_XYZ::MultipliedMat " + " ".join(a) + " " + mat9(rng))

        lines.append("gp_Mat::Multiply " + mat9(rng) + " " + mat9(rng))
        lines.append("gp_Mat::Determinant " + mat9(rng))
        lines.append("gp_Mat::Transposed " + mat9(rng))

        lines.append("gp_Vec::Dot " + " ".join(a + b))
        lines.append("gp_Vec::Crossed " + " ".join(a + b))
        lines.append("gp_Vec::Magnitude " + " ".join(a))

        # gp_Dir::Crossed
        da = [repr(rng.random()*2-1) for _ in range(3)]
        db = [repr(rng.random()*2-1) for _ in range(3)]
        lines.append("gp_Dir::Crossed " + " ".join(da + db))

        # ---- 2D ops ----
        a2, b2 = vec2(rng, scale), vec2(rng, scale)
        lines.append("gp_XY::Dot " + " ".join(a2 + b2))
        lines.append("gp_XY::Crossed " + " ".join(a2 + b2))
        lines.append("gp_XY::Modulus " + " ".join(a2))
        lines.append("gp_XY::Added " + " ".join(a2 + b2))
        lines.append("gp_XY::MultipliedScalar " + " ".join(a2 + [repr(rng.random()*10-5)]))
        if rng.random() > 0.1:
            lines.append("gp_XY::Normalized " + " ".join(a2))
        lines.append("gp_Vec2d::Dot " + " ".join(a2 + b2))
        lines.append("gp_Vec2d::Magnitude " + " ".join(a2))
        lines.append("gp_Vec2d::Crossed " + " ".join(a2 + b2))

        d2a = [repr(rng.random()*2-1) for _ in range(2)]
        d2b = [repr(rng.random()*2-1) for _ in range(2)]
        lines.append("gp_Dir2d::new " + " ".join(d2a))
        lines.append("gp_Dir2d::Dot " + " ".join(d2a + d2b))

    return "\n".join(lines) + "\n"

def main():
    print("occt: building driver (C++ OCCT) + runner (pure Rust)...")
    driver = differ.build(ROOT / "driver")
    runner = differ.build(ROOT / "runner")
    corpus = make_corpus()
    golden = differ.run(driver, corpus, DLL_DIRS)
    cand = differ.run(runner, corpus)

    mismatches, oks, errs, per_op = differ.diff(corpus.splitlines(), golden, cand, OPS)
    differ.write_results(ROOT / "last-run.json", per_op)
    print(f"records={len(corpus.splitlines())}  ok-matched={oks}  err-matched={errs}  mismatches={len(mismatches)}")
    for i, line, g, c in mismatches[:10]:
        print(f"  [{i}] {line}\n      occt: {g}\n      rust: {c}")
    if mismatches:
        sys.exit(1)
    print("PASS")

if __name__ == "__main__":
    main()
