//! `GeomFill_Generator` — the ruled surface through the profiled sections.
//!
//! Source: `GeomFill_Generator.cxx:32-82`.

use std::sync::Arc;

use occt_core::precision::PCONFUSION;

use crate::curve::Curve;
use crate::geom_fill::profiler::Profiler;
use crate::surface::Surface;

/// `GeomFill_Generator::Perform` (`GeomFill_Generator.cxx:32-82`) after
/// `GeomFill_Profiler::AddCurve` of each section.
pub fn generator_surface(sections: &[Arc<dyn Curve>]) -> Result<Arc<dyn Surface>, String> {
    let mut profiler = Profiler::default();
    for curve in sections {
        profiler.add_curve(curve.as_ref())?;
    }
    profiler.perform(PCONFUSION)?;
    Ok(Arc::new(profiler.surface()?))
}
