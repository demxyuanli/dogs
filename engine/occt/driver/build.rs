fn main() {
    let occt = std::env::var("MIG_ROOT").unwrap_or_else(|_| {
        std::fs::read_to_string("../.mig_root")
            .unwrap_or_else(|_| "d:/source/occt-src".into())
            .trim()
            .to_string()
    });
    let src = format!("{occt}/src/FoundationClasses");

    let fc_dirs = [
        "TKernel", "TKernel/Standard", "TKernel/Precision", "TKernel/NCollection",
        "TKernel/TCollection", "TKernel/TColStd", "TKernel/Message", "TKernel/OSD",
        "TKernel/Quantity", "TKernel/Resource", "TKernel/Storage", "TKernel/Units",
        "TKernel/UnitsAPI", "TKernel/UnitsMethods", "TKernel/FSD", "TKernel/StdFail",
        "TKMath", "TKMath/gp", "TKMath/math", "TKMath/Bnd", "TKMath/BSplCLib",
        "TKMath/BSplSLib", "TKMath/BVH", "TKMath/Convert", "TKMath/CSLib",
        "TKMath/ElCLib", "TKMath/ElSLib", "TKMath/PLib", "TKMath/Poly", "TKMath/TopLoc",
        "TKMath/GeomAbs", "TKMath/MathInteg", "TKMath/MathLin", "TKMath/MathOpt",
        "TKMath/MathPoly", "TKMath/MathRoot", "TKMath/MathSys", "TKMath/MathUtils",
    ];

    let mut build = cc::Build::new();
    build.cpp(true).std("c++17");
    build.define("Standard_EXPORT", "");

    for dir in &fc_dirs {
        build.include(format!("{}/{}", src, dir));
    }

    // Compile OCCT .cxx files needed for Standard_Failure + friends
    build.file(format!("{src}/TKernel/Standard/Standard_Failure.cxx"));
    build.file(format!("{src}/TKernel/Standard/Standard.cxx"));       // Standard::Free, Allocate
    build.file(format!("{src}/TKernel/Standard/Standard_OutOfMemory.cxx"));
    build.file(format!("{src}/TKernel/Standard/Standard_StackTrace.cxx"));

    build.file("cpp/driver.cpp").compile("occt_driver");
    println!("cargo:rerun-if-changed=cpp/driver.cpp");
}
