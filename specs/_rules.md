# _rules — OCCT → Rust migration rulebook

> Grill-me completed: 2026-07-27 (AI survey phase; user grill-me pending).
> Rules below are derived from FoundationClasses source analysis.
> Higher-layer rules (ModelingData, ModelingAlgorithms, etc.) deferred to plan.

## Type map

> **Key insight:** OCCT 7.9+ deprecates ALL Standard_* typedefs. Port directly to Rust native types, not through intermediate OCCT type wrappers.

| OCCT type | Rust | Source | Notes |
|---|---|---|---|
| `Standard_Real` | `f64` | `Standard_TypeDef.hxx:74-76` | Deprecated typedef for `double` |
| `Standard_Integer` | `i32` | `Standard_TypeDef.hxx:67-68` | Deprecated typedef for `int` (32-bit on all platforms OCCT supports) |
| `Standard_Boolean` | `bool` | `Standard_TypeDef.hxx:79-80` | Deprecated typedef for `bool` |
| `Standard_Byte` | `u8` | `Standard_TypeDef.hxx:91-92` | Deprecated typedef for `uint8_t` |
| `Standard_Size` | `usize` | `Standard_TypeDef.hxx:99-100` | Deprecated typedef for `size_t` |
| `Standard_CString` | `&CStr` / `*const c_char` | `Standard_TypeDef.hxx:131-132` | Deprecated typedef for `const char*` |
| `Standard_ShortReal` | `f32` | `Standard_TypeDef.hxx:83-84` | Deprecated typedef for `float` |
| `Standard_Character` | `i8` / `c_char` | `Standard_TypeDef.hxx:87-88` | Deprecated typedef for `char` |
| `Standard_Utf8Char` | `c_char` | `Standard_TypeDef.hxx:107-108` | Deprecated typedef for `char` |
| `Standard_Utf16Char` | `u16` (`char16_t`) | `Standard_TypeDef.hxx:119-120` | |
| `Standard_Utf32Char` | `char` (`char32_t`) | `Standard_TypeDef.hxx:123-124` | |
| `gp_XYZ` | `struct GpXyz { x: f64, y: f64, z: f64 }` | `gp_XYZ.hxx:33-34` | 3D cartesian coordinate, all constexpr ops |
| `gp_Pnt` | `struct GpPnt(GpXyz)` | `gp_Pnt.hxx` | 3D point (newtype over XYZ) |
| `gp_Vec` | `struct GpVec(GpXyz)` | `gp_Vec.hxx` | 3D vector (newtype over XYZ) |
| `gp_Dir` | `struct GpDir(GpXyz)` | `gp_Dir.hxx` | Unit 3D direction (normalized XYZ) |
| `gp_Mat` | `struct GpMat([[f64; 3]; 3])` | `gp_Mat.hxx` | 3×3 matrix, row-major |
| `gp_Trsf` | `struct GpTrsf { matrix: GpMat, translation: GpXyz }` | `gp_Trsf.hxx` | Affine 3D transform (3×3 + translation) |
| `gp_Quaternion` | `struct GpQuaternion { x: f64, y: f64, z: f64, w: f64 }` | `gp_Quaternion.hxx` | Quaternion rotation |
| `Handle(T)` / `opencascade::handle<T>` | `Arc<T>` or `Rc<T>` | `Standard_Handle.hxx:53` | Intrusive refcounted smart pointer. T must extend Standard_Transient. |
| `Standard_Transient` | trait + `Arc<T>` (no custom base class needed) | `Standard_Transient.hxx:43` | Base for refcounted objects. Provides RTTI (DynamicType, IsKind). |
| `NCollection_Array1<T>` | `Vec<T>` | `NCollection_Array1.hxx` | 1-indexed contiguous array |
| `NCollection_Sequence<T>` | `Vec<T>` or `VecDeque<T>` | `NCollection_Sequence.hxx` | Indexed sequence |
| `NCollection_DataMap<K,V>` | `HashMap<K,V>` | `NCollection_DataMap.hxx` | Hash map |
| `NCollection_List<T>` | `LinkedList<T>` or `Vec<T>` | `NCollection_List.hxx` | Doubly-linked list |

## Constant map

> Source: `Precision.hxx` and `Standard_Real.hxx`. All values are exact — no platform variation.

| OCCT constant | Value | Rust | Source | Notes |
|---|---|---|---|---|
| `Precision::Confusion()` | `1e-7` | `const CONFUSION: f64 = 1e-7` | `Precision.hxx:165` | Default 3D point coincidence tolerance |
| `Precision::Angular()` | `1e-12` | `const ANGULAR: f64 = 1e-12` | `Precision.hxx:123` | Angle equality tolerance (radians) |
| `Precision::Intersection()` | `1e-9` | `const INTERSECTION: f64 = 1e-9` | `Precision.hxx:220` | Confusion() / 100 |
| `Precision::Approximation()` | `1e-6` | `const APPROXIMATION: f64 = 1e-6` | `Precision.hxx:235` | Confusion() * 10 |
| `Precision::Infinite()` | `2e+100` | `const INFINITE: f64 = 2e100` | `Precision.hxx:371` | Sentinel for "infinite" values |
| `Precision::Computational()` | `DBL_EPSILON` ≈ `2.22e-16` | `f64::EPSILON` | `Precision.hxx:192` | Machine epsilon for numerical guards |
| `Precision::PConfusion()` | `1e-9` | `const PCONFUSION: f64 = 1e-9` | `Precision.hxx:334` | Parametric confusion (Confusion/100) |
| `RealSmall()` | `DBL_MIN` ≈ `2.225e-308` | `f64::MIN_POSITIVE` | `Standard_Real.hxx:132-135` | Minimum positive double |
| `RealEpsilon()` | `DBL_EPSILON` ≈ `2.22e-16` | `f64::EPSILON` | `Standard_Real.hxx:161-164` | Machine epsilon |
| `gp::Resolution()` | `1e-12` (from gp.hxx) | `const RESOLUTION: f64 = 1e-12` | `gp.hxx` | Spatial tolerance for gp package |

## Exception map

> Source: `Standard_Failure.hxx`, `Standard_*Error.hxx`, macro definitions in `Standard_Macro.hxx`.

| OCCT exception | Condition | Rust equivalent | Source |
|---|---|---|---|
| `Standard_Failure` | Root of hierarchy | `Box<dyn Error>` | `Standard_Failure.hxx:28` |
| `Standard_ConstructionError` | Invalid construction args (e.g. zero-norm Dir) | `Result::Err` or `panic!` | `Standard_ConstructionError.hxx` |
| `Standard_DomainError` | Value outside function domain | `Result::Err` | `Standard_DomainError.hxx` |
| `Standard_OutOfRange` | Index out of 1..3 range for Coord/SetCoord | `Result::Err` | `Standard_OutOfRange.hxx` |
| `Standard_NullObject` | Null handle dereference | `Option::None` / `Result::Err` | `Standard_NullObject.hxx` |
| `Standard_NoSuchObject` | Key not found in collection | `Option::None` | `Standard_NoSuchObject.hxx` |
| `Standard_DivideByZero` | Division by zero in algorithm | `Result::Err` | `Standard_DivideByZero.hxx` |
| `Standard_NumericError` | Numeric computation failure | `Result::Err` | `Standard_NumericError.hxx` |
| `Standard_Overflow` / `Standard_Underflow` | Float overflow/underflow | `Result::Err` | respective headers |

**Exception mapping rule:**
- OCCT uses `*_Raise_if(cond, msg)` macros that throw on condition. Rust uses `Result<T, E>`.
- Functions that can fail → return `Result`. Pure math on value types that panic in OCCT → keep as `Result` in Rust.
- The `Standard_Failure::what()` message string → `Error::msg()` or display impl.

## Idiom map

| OCCT idiom | Rust translation | Source pattern |
|---|---|---|
| `Handle(T)` intrusive refcount | `Arc<T>` (stdlib) | `Standard_Handle.hxx` — custom intrusive_ptr, but Arc is equivalent for migration |
| `Standard_Transient` base class | No base needed; objects in `Arc` are automatically refcounted | `Standard_Transient.hxx:43` |
| `DynamicType()` / `IsKind()` RTTI | `Any` + `downcast_ref`, or enum-based dispatch | `Standard_Transient.hxx:79-95` |
| 1-indexed Coord/SetCoord | 0-indexed `[usize]` access or named fields `.x`, `.y`, `.z` | `gp_XYZ.hxx:64-83` |
| `IsEqual(tolerance)` methods | `approx_eq(a, b, tol)` free function or trait | `gp_XYZ.hxx:164-168` |
| `DEFINE_STANDARD_ALLOC` macro | Drop — Rust uses global allocator | `Standard_DefineAlloc.hxx` |
| `Standard_EXPORT` (DLL export) | Drop — Rust crate visibility | throughout |
| `*_Raise_if(cond, msg)` macros | `if cond { return Err(...) }` | `Standard_OutOfRange.hxx` |
| `DumpJson` serialization | `serde::Serialize` | `gp_XYZ.hxx:510` |
| `InitFromJson` deserialization | `serde::Deserialize` | `gp_XYZ.hxx:513` |
| `constexpr` everywhere in gp_* | Rust `const fn` where possible | throughout gp headers |
| Out-parameters `Coord(&x, &y, &z)` | Tuple return `(f64, f64, f64)` or destructure | `gp_XYZ.hxx:128-132` |
| `RealToInt()` with clamp | `(val as i32).clamp(i32::MIN, i32::MAX)` or `val as i32` (Rust saturates by default in debug, wraps in release — use explicit clamp) | `Standard_Real.hxx:319-328` |

## Footguns

| # | Issue | Source | Mitigation |
|---|---|---|---|
| F1 | **1-based indexing** — `SetCoord(1,x)` maps to X, `Coord(1)` returns X. Porting mindlessly to `[0]`/`[1]`/`[2]` produces off-by-one bugs. | `gp_XYZ.hxx:64-83` | Named fields `.x`, `.y`, `.z` — never index-based access |
| F2 | **DBL_MIN as small tolerance** — `RealSmall()` = DBL_MIN ≈ 2.2e-308. Using this for `IsEqual` means only exact equality passes. OCCT uses it as "effectively zero" check but DBL_MIN is subnormal. | `Standard_Real.hxx:132-135` | Use `Precision::Computational()` = f64::EPSILON for numerical guards, not DBL_MIN |
| F3 | **RealToInt clamps instead of UB** — OCCT's `RealToInt()` clamps out-of-range doubles to INT_MIN/INT_MAX because MSVC raises a trap on direct cast. Rust's `as i32` is defined differently (saturates in debug, wraps in release). | `Standard_Real.hxx:319-328` | Use `val.clamp(i32::MIN as f64, i32::MAX as f64) as i32` for equivalent behavior |
| F4 | **handle::operator->() is non-const** — `const handle<T>` still gives mutable `T*` access. No const correctness through handles. | `Standard_Handle.hxx:125` | Rust's `Arc<T>` gives `&T` from `Deref` — const-correct by default; use `Arc::get_mut` or interior mutability for mutation |
| F5 | **DumpJson precision** — `DumpJson` uses default ostream precision which may truncate doubles. The geometry reference harness notes `setprecision(17)` is required for round-trip fidelity. | `gp_XYZ.cxx` (impl) | Use `serde_json` with full precision (Rust default is lossless round-trip for f64) |
| F6 | **gp_XYZ::Cross() aliasing** — `Cross()` stores result in-place, must cache `this->x` and `this->y` before overwriting `this->z`. Port must match the same aliasing pattern. | `gp_XYZ.hxx:523-530` | `let (nx, ny, nz) = (self.y*rhs.z - self.z*rhs.y, ...);` then assign — Rust's let bindings naturally avoid aliasing |

## Deferred

- ModelingData layer (TKG2d, TKG3d, TKGeomBase, TKBRep) — geometry curves/surfaces, topology, B-rep data structures
- ModelingAlgorithms layer — boolean ops, fillets, offsets, mesh, hidden line removal
- DataExchange layer — STEP, IGES, STL, glTF, VRML parsers/writers
- Visualization layer — OpenGL, 3D viewer
- ApplicationFramework layer — OCAF document framework
- NCollection → Rust stdlib container mapping (detailed per-container analysis)
- Thread safety of handle<T> — atomic refcount but object access unprotected
- `Standard_MMgrOpt` custom memory manager — can be replaced with Rust's global allocator
