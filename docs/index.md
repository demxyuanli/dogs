# docs index — OCCT → Rust migration

> cross-session entry point. keyword → topic file.

## current status

| part | status |
|---|---|
| source scan / domain classification | ✅ done (trellis-init) |
| engine domain scaffold | ✅ done |
| echo op verified | ✅ PASS |
| migration index seeded | ✅ done (24 units, codegraph-verified deps) |
| code flow analysis | ✅ done (codegraph on d:/source/occt-src/.codegraph/) |
| rulebook (Layer 0) | ✅ AI survey done; grill-me pending |
| strategic plan | ○ pending (trellis-plan) |

## topic files

| file | content |
|---|---|
| [specs/_rules.md](../specs/_rules.md) | project rulebook — types, constants, exceptions, idioms, footguns |
| [specs/_tasks.md](../specs/_tasks.md) | atomic task matrix — 27 atoms across 2 layers (Layer 0-1) |
| [code-flow.json](../code-flow.json) | static call graph from codegraph — 5 dependency layers, 8 cross-cutting edges, blast radius |
| [migration-index.json](../migration-index.json) | machine-readable index (24 units, content-hashed) |
| [trellis/plan.md](trellis/plan.md) | strategic migration map (run trellis-plan) |

## key points

- original source: `d:/source/occt-src`
- domain kind: 3D CAD kernel (geometry + topology + algorithms) / default oracle: `numeric_tol`
- engine domain: `engine/occt/`
- build system: CMake (C++17+)
- scale: ~13,395 source files, ~73 toolkits, 7 module groups
- root layer: FoundationClasses (TKernel, TKMath) — 21 atoms identified
