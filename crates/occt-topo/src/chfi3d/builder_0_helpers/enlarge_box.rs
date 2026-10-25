//! `ChFi3d_EnlargeBox` — grow a `Bnd_Box` from a parameter of a curve, a
//! p-curve, an edge or a `ChFiDS_SurfData` end.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:335-472`
//! (OCCT 8.0.0):
//!
//! - `ChFi3d_EnlargeBox(Geom_Curve, wd, wf, box1, box2)` (`cxx:335-343`),
//! - `ChFi3d_EnlargeBox(Adaptor3d_Surface, Geom2d_Curve, wd, wf, box1, box2)`
//!   (`cxx:347-358`),
//! - `ChFi3d_EnlargeBox(Edge, FaceList, w, box)` (`cxx:363-382`),
//! - `ChFi3d_EnlargeBox(DataStructure, Stripe, SurfData, b1, b2, isfirst)`
//!   (`cxx:385-472`), whose file-static helper `Geometry(DStr, ind)`
//!   (`cxx:289-313`) picks the surface of an index of the DS.
//!
//! Translation notes:
//! - the four overloads become [`ch_fi3d_enlarge_box_curve`],
//!   [`ch_fi3d_enlarge_box_pcurve`], [`ch_fi3d_enlarge_box_edge`] and
//!   [`ch_fi3d_enlarge_box_surf_data`],
//! - `Bnd_Box&` out-parameters are `&mut BndBox`,
//! - `Geometry` answers the surface the DS index stands for: a positive index
//!   is a face (`BRepAdaptor_Surface(F, false)` in OCCT, the face's own
//!   surface here because only `Value(u, v)` is ever called on it), a negative
//!   one the `Geom_Surface` of `Surface(-ind)`, and a null handle for index
//!   `0`,
//! - the `BRepAdaptor_Curve::Initialize(E, F)` of the edge overload becomes
//!   the p-curve of `E` on `F` evaluated on `F`'s surface, which is what that
//!   adaptor evaluates (`BRepAdaptor_Curve.cxx`, the `CurveOnSurface` arm).

use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_geom::curve::Curve as GeomCurve;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_tool::BRepTool;
use crate::chfi_ds::SurfData;
use crate::chfi_ds::Stripe;
use crate::shape::{Edge, Face};
use crate::shape::TopoShape;
use crate::topopebrep::ds::DataStructure;

/// `ChFi3d_EnlargeBox(Geom_Curve, wd, wf, box1, box2)` (`cxx:335-343`).
pub fn ch_fi3d_enlarge_box_curve(
    c: &Arc<dyn GeomCurve>,
    wd: f64,
    wf: f64,
    box1: &mut BndBox,
    box2: &mut BndBox,
) {
    box1.add_point(&c.d0(wd));
    box2.add_point(&c.d0(wf));
}

/// `ChFi3d_EnlargeBox(Adaptor3d_Surface, Geom2d_Curve, wd, wf, box1, box2)`
/// (`cxx:347-358`).
pub fn ch_fi3d_enlarge_box_pcurve(
    s: &Arc<dyn Surface>,
    pc: &Arc<dyn Curve2d>,
    wd: f64,
    wf: f64,
    box1: &mut BndBox,
    box2: &mut BndBox,
) {
    let uv = pc.d0(wd);
    box1.add_point(&s.d0(uv.x(), uv.y()));
    let uv = pc.d0(wf);
    box2.add_point(&s.d0(uv.x(), uv.y()));
}

/// `Geometry(DStr, ind)` (`cxx:289-313`): the surface of a DS index, `None`
/// for index `0` and for an unbound entry.
fn geometry(ds: &DataStructure, ind: i32) -> Option<Arc<dyn Surface>> {
    if ind == 0 {
        return None;
    }
    if ind > 0 {
        // `TopoDS::Face(DStr.Shape(ind))` then
        // `BRepAdaptor_Surface(F, false)` (`cxx:296-302`).
        let shape = ds.shape(ind, false)?;
        return BRepTool::face_surface(&Face(shape));
    }
    // `DStr.Surface(-ind).Surface()` then `new GeomAdaptor_Surface(S)`
    // (`cxx:307-311`).
    ds.surface(-ind).surface().cloned()
}

/// `ChFi3d_EnlargeBox(Edge, FaceList, w, box)` (`cxx:363-382`).
///
/// The edge itself is always evaluated first; each face of the list then adds
/// the same parameter through the edge's p-curve on that face.
pub fn ch_fi3d_enlarge_box_edge(
    e: &TopoShape,
    lf: &[TopoShape],
    w: f64,
    bx: &mut BndBox,
) {
    let edge = Edge(e.clone());
    // `BRepAdaptor_Curve BC(E); box.Add(BC.Value(w));` (`cxx:369-370`).
    if let Some(p) = crate::topopebrep::tool::adaptor::CurveAdaptor::of_edge(&edge).d0(w) {
        bx.add_point(&p);
    }
    for shape in lf {
        // `TopoDS_Face F = TopoDS::Face(It.Value()); if (!F.IsNull()) ...`
        // (`cxx:375-380`).
        let face = Face(shape.clone());
        let Some(surface) = BRepTool::face_surface(&face) else {
            continue;
        };
        let Some((pc, _, _)) = curve_on_surface_range(&edge, &face) else {
            continue;
        };
        let cons = CurveOnSurface::new(pc, surface);
        bx.add_point(&cons.d0(w));
    }
}

/// `ChFi3d_EnlargeBox(DStr, st, sd, b1, b2, isfirst)` (`cxx:385-472`).
///
/// `st` is the `ChFiDS_Stripe` the surf data belongs to; OCCT tests it for
/// null and the port takes `Option`. `DStr.Surface(sd->Surf()).Surface()` is
/// read once at `cxx:399` and dereferenced by every p-curve branch; a null
/// handle there is the OCCT null dereference and answers `Err`.
pub fn ch_fi3d_enlarge_box_surf_data(
    ds: &DataStructure,
    st: Option<&Stripe>,
    sd: &SurfData,
    b1: &mut BndBox,
    b2: &mut BndBox,
    isfirst: bool,
) -> Result<(), String> {
    // `cxx:396-397`.
    let cp1 = sd.vertex(isfirst, 1);
    let cp2 = sd.vertex(isfirst, 2);
    b1.add_point(cp1.point());
    b2.add_point(cp2.point());

    // `cxx:398-404`.
    let fi1 = sd.interference_on_s1();
    let fi2 = sd.interference_on_s2();
    let surface = ds.surface(sd.surf()).surface().cloned();
    let null_surface = || {
        "ChFi3d_EnlargeBox: DStr.Surface(Surf()).Surface() is null".to_string()
    };
    let pcs1 = fi1.p_curve_on_surf();
    let pcs2 = fi2.p_curve_on_surf();
    let c3d1 = ds.curve(fi1.line_index()).curve();
    let c3d2 = ds.curve(fi2.line_index()).curve();
    let f1 = geometry(ds, sd.index_of_s1());
    let f2 = geometry(ds, sd.index_of_s2());

    // `cxx:407-414`.
    let mut p1 = fi1.parameter(isfirst);
    if let Some(c3d1) = c3d1 {
        b1.add_point(&c3d1.d0(p1));
    }
    if let Some(pcs1) = pcs1 {
        let s = surface.as_ref().ok_or_else(null_surface)?;
        let uv = pcs1.d0(p1);
        b1.add_point(&s.d0(uv.x(), uv.y()));
    }
    if let Some(f1) = &f1 {
        if let Some(pcf1) = fi1.p_curve_on_face() {
            let uv = pcf1.d0(p1);
            b1.add_point(&f1.d0(uv.x(), uv.y()));
        }
    }

    // `cxx:425-444`.
    let mut p2 = fi2.parameter(isfirst);
    if let Some(c3d2) = c3d2 {
        b2.add_point(&c3d2.d0(p2));
    }
    if let Some(pcs2) = pcs2 {
        let s = surface.as_ref().ok_or_else(null_surface)?;
        let uv = pcs2.d0(p2);
        b2.add_point(&s.d0(uv.x(), uv.y()));
    }
    if let Some(f2) = &f2 {
        if let Some(pcf2) = fi2.p_curve_on_face() {
            let uv = pcf2.d0(p2);
            b2.add_point(&f2.d0(uv.x(), uv.y()));
        }
    }

    // `cxx:445-468`.
    let Some(st) = st else {
        return Ok(());
    };
    let c3d = ds.curve(st.curve(isfirst)).curve();
    let c2d = st.p_curve(isfirst);
    // `if (st->Orientation(isfirst) == TopAbs_FORWARD) st->Parameters(isfirst,
    // p1, p2); else st->Parameters(isfirst, p2, p1);`
    let (pdeb, pfin) = st.parameters(isfirst);
    if st.orientation(isfirst) == Orientation::Forward {
        p1 = pdeb;
        p2 = pfin;
    } else {
        p1 = pfin;
        p2 = pdeb;
    }
    if let Some(c3d) = c3d {
        b1.add_point(&c3d.d0(p1));
        b2.add_point(&c3d.d0(p2));
    }
    if let Some(c2d) = c2d {
        let s = surface.as_ref().ok_or_else(null_surface)?;
        let uv = c2d.d0(p1);
        b1.add_point(&s.d0(uv.x(), uv.y()));
        let uv = c2d.d0(p2);
        b2.add_point(&s.d0(uv.x(), uv.y()));
    }
    Ok(())
}
