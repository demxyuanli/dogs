//! `BRepFill_Draft::Init` (`BRepFill_Draft.cxx:412-460`): the location law of
//! the draft and the straight generating line.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::gp::{GpDir, GpLin, GpPnt, GpVec, GpXyz};
use occt_geom::curve::Curve;
use occt_geom::geom_fill::location_draft::LocationDraft;
use occt_geom::geom_fill::LocationLaw as GeomFillLocationLaw;
use occt_geom::line::GeomLine;

use crate::brep_fill_draft_law::DraftLaw;
use crate::brep_fill_section_law::SectionLaw;
use crate::brep_fill_shape_law::ShapeLaw;
use crate::brep_lib_make_edge;
use crate::brep_lib_make_wire::MakeWire;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Wire};

use super::geometry;
use super::BrepFillDraft;

impl BrepFillDraft {
    /// `Init(Surface, Length, Box)` (`cxx:412-460`).
    ///
    /// The surface argument is unused in the reference (its parameter has no
    /// name) and is therefore not taken here either.
    pub(super) fn init(&mut self, length: f64, b: &BndBox) -> Result<(), String> {
        // `cxx:416-417`: the positioning law.
        let loc = Arc::new(LocationDraft::new(&self.my_dir, self.my_angle));
        let handle: Arc<dyn GeomFillLocationLaw> = loc.clone();
        let mut law = DraftLaw::new(&self.my_wire, &handle);

        // `cxx:419-426`.
        let good = geometry::good_orientation(b, law.base(), &self.my_dir)?;
        if self.is_internal ^ (!good) {
            self.my_angle = -self.my_angle;
            loc.set_angle(self.my_angle)?;
            law = DraftLaw::new(&self.my_wire, &handle);
        }

        // `cxx:429`: "Clean small discontinuities."
        law.clean_law(self.angmin)?;

        // `cxx:434-448`: the generating direction starts along `+Y` and is
        // flipped when the law's binormal points away from `myDir`.
        let mut d = GpVec::from_xyz(&GpXyz::new(0.0, 1.0, 0.0));
        let law1 = law
            .base()
            .law(1)
            .ok_or_else(|| "Standard_NullObject: BRepFill_DraftLaw::Law(1)".to_string())?;
        let (f, l) = law1.get_domain();
        let (m, _bid) = law1.d0((f + l) / 2.0)?;
        let bn = GpDir::from_xyz(&m.column(2)).map_err(|e| format!("gp_Dir: {e}"))?;
        if self.my_dir.angle(&bn) > PI / 2.0 {
            d.reverse();
        }
        let dir = GpDir::from_vec(&d).map_err(|e| format!("gp_Dir: {e}"))?;

        // `cxx:449-458`: `Geom_Line(P, D)` trimmed to `[0, Length]`, one edge,
        // one wire.
        let line: Arc<dyn Curve> = Arc::new(GeomLine::new(GpLin::from_pnt_dir(
            GpPnt::new(0.0, 0.0, 0.0),
            dir,
        )));
        let eg = make_section_edge(&line, length)?;
        let mut mw = MakeWire::new();
        mw.add(&eg)
            .map_err(|_| "BRepLib_MakeWire: not done".to_string())?;
        if !mw.is_done() {
            return Err("BRepLib_MakeWire: not done".to_string());
        }
        let g: Wire = mw.wire();

        // `cxx:459`: `mySec = new BRepFill_ShapeLaw(G, true)`.
        let sec = ShapeLaw::with_wire(&g, true)?;
        self.my_sec = Some(Arc::new(sec) as Arc<dyn SectionLaw>);
        self.my_loc = Some(Arc::new(law));
        Ok(())
    }
}

/// `handle(Geom_TrimmedCurve) TC = new Geom_TrimmedCurve(L, 0, Length);` then
/// `BRepLib_MakeEdge ME(TC)` (`cxx:451-453`).
///
/// The reference stores the trimmed handle on the edge and gives it the range
/// `[0, Length]`; `Geom_TrimmedCurve` keeps basis parameters, so those two
/// numbers measure the same segment. The port's
/// [`GeomTrimmedCurve`](occt_geom::trimmed::GeomTrimmedCurve) reparameterises
/// its trim onto `[0, 1]`, which would make a stored range of `[0, Length]`
/// meaningless, so the edge is built straight on the basis line over
/// `[0, Length]` — the same curve, the same range, the same evaluated points.
fn make_section_edge(line: &Arc<dyn Curve>, length: f64) -> Result<Edge, String> {
    // `BRepLib_MakeEdge(TC)` = `Init(TC)` with the curve's own bounds
    // (`BRepLib_MakeEdge.cxx:606-610`, `:604`).
    let build = brep_lib_make_edge::init_3d(line, 0.0, length)
        .ok_or_else(|| "BRepLib_MakeEdge: not done".to_string())?;
    let b = TopoBuilder::new();
    let mut e = b.make_edge(line.clone(), build.first, build.last);
    brep_lib_make_edge::add_endpoint_vertices(&b, &mut e, &build);
    // `B.Range(E, p1, p2)` / `E.Closed(closed)` (`cxx:791-792`).
    if build.closed {
        e.0.set_closed(true);
    }
    Ok(e)
}
