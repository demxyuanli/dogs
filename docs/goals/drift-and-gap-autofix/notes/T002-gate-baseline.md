# T002 receipt — gate baseline frozen 2026-10-15 (HEAD e77be18e, clean tree)

## Commands and results

```
rtk cargo check --manifest-path crates/occt-math/Cargo.toml   --offline                  -> 0 errors,  27 warnings
rtk cargo check --manifest-path crates/occt-geom2d/Cargo.toml --offline                 -> 0 errors,  81 warnings
rtk cargo check --manifest-path crates/occt-geom/Cargo.toml   --offline                 -> 0 errors, 100 warnings
rtk cargo check --manifest-path crates/occt-core/Cargo.toml   --offline                 -> 0 errors,   8 warnings
rtk cargo check --manifest-path crates/occt-topo/Cargo.toml   --offline --all-targets   -> 0 errors, 375 warnings
rtk cargo test  --manifest-path crates/occt-topo/Cargo.toml   --offline --lib           -> 1251 passed, 0 failed (4.50s)
rtk cargo test  --manifest-path crates/occt-topo/Cargo.toml   --offline --test step_obj_gates -- --nocapture
                                                                                         -> 4 passed, 0 failed (182.12s)
rtk cargo run   --manifest-path crates/occt-topo/Cargo.toml   --offline --example export_data_obj
                                                                                         -> 15 ok, 0 err
```

## `export_data_obj` baseline (bare run, 15 models)

| model | v | f | exact bbox |
|---|---|---|---|
| Cone | 195 | 301 | [-4,-4,-0] ~ [4,4,10] |
| Cube | 24 | 12 | [-0,-0,-0] ~ [10,10,10] |
| Cylinder | 146 | 140 | [-2,-2,-0] ~ [2,2,10] |
| Extrusion | 24 | 12 | [-0,-0,-0] ~ [20,5,15] |
| HoledPlate | 180 | 128 | [-1.072,-1.072,-1.072] ~ [101.072,81.072,6.072] |
| Offset | 712 | 892 | [-2,-2,-2] ~ [12,12,12] |
| OffsetPlaneHoleEdge | 48 | 32 | [-2,-2,-2] ~ [12,12,4] |
| Shape-1 | 3335 | 4328 | [-25.772,-17.738,-0] ~ [25.772,17.738,77] |
| Shape-2 | 3081 | 4744 | [-63.958,-66.299,-0] ~ [63.958,66.299,150] |
| Shape | 6150 | 11372 | [-35,-35,-5] ~ [35,35,107.106] |
| Sphere | 642 | 1244 | [-5,-5,-5] ~ [5,5,5] |
| Torus | 1369 | 2592 | [-12.989,-12.989,-2] ~ [12.989,12.989,2] |
| linkrods | 3494 | 5078 | [3.063,2.438,-0] ~ [8.178,4.062,2] |
| rev | 104 | 92 | [-0,-0,-0] ~ [10,10,10] |
| screw | 601 | 792 | [-27.898,-10.827,-34.564] ~ [-7.898,9.174,7.751] |

## Corrections to the recorded baselines

1. The charter's blind spot "26 lib tests fail on `%TEMP%` write permission" did **not** reproduce:
   the run is 1251 / 0. A future lib failure must be treated as a real regression until proven
   environmental.
2. `export_data_obj` prints 15 models for the bare run; the 23/23 figure quoted in the drift canvas
   applies to the `-- data/occ data` argument form (15 `data/*.step` + 8 `data/occ/*`).
3. The drift canvas `GATES` array still carries `G-lib 1276/0`; the measured value is 1251/0. The
   canvas needs a pass (side finding).

## Artifacts

- `.target-gate/t000/lib.txt`
- `.target-gate/t000/gates.txt`
- `.target-gate/t000/export-obj.txt`
- `.target-gate/t000/check.txt`
