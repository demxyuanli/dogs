//! Port of `BRepClass_Intersector` (`BRepClass_Intersector.cxx`). The OCCT
//! class derives from `IntRes2d_Intersection`; here the result is held as
//! the `base` field and exposed through the same accessors.

use occt_core::bnd::BndBox2d;
use occt_core::gp::{GpDir2d, GpLin2d};
use occt_core::intres2d::{
    IntRes2dDomain, IntRes2dIntersection, IntRes2dIntersectionPoint,
    IntRes2dIntersectionSegment,
};
use occt_core::precision::{Precision, CONFUSION, PCONFUSION, PINTERSECTION};
use occt_geom2d::cl_props2d::GeomLPropCLProps2d;
use occt_geom2d::geom2d_int::Geom2dIntGInter;
use occt_geom2d::line::Geom2dLine;

use crate::bnd_lib_add2d::add_geom2d_range;
use crate::boptools_2d::curve_on_surface_range;
use crate::topo_tools_full::edge_vertices;

use super::edge::BRepClassEdge;
use super::intersector_checks::{
    adaptor_curve, check_on, check_skip, get_tangent_as_chord, is_inter, lin_value,
    max_tol_2d_cur_edge,
};

/// `BRepClass_Intersector`.
pub struct BRepClassIntersector {
    base: IntRes2dIntersection,
}

impl Default for BRepClassIntersector {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepClassIntersector {
    /// `BRepClass_Intersector()`.
    pub fn new() -> Self {
        Self {
            base: IntRes2dIntersection::new(),
        }
    }

    /// `IsDone()`.
    pub fn is_done(&self) -> bool {
        self.base.is_done()
    }

    /// `NbPoints()`.
    pub fn nb_points(&self) -> usize {
        self.base.nb_points()
    }

    /// `Point(N)`, 1-based as in OCCT.
    pub fn point(&self, n: usize) -> &IntRes2dIntersectionPoint {
        self.base.point(n)
    }

    /// `NbSegments()`.
    pub fn nb_segments(&self) -> usize {
        self.base.nb_segments()
    }

    /// `Segment(N)`, 1-based as in OCCT.
    pub fn segment(&self, n: usize) -> &IntRes2dIntersectionSegment {
        self.base.segment(n)
    }

    /// `Append(P)` (`IntRes2d_Intersection::Append`): appends without
    /// clearing the list.
    pub fn append(&mut self, pnt: &IntRes2dIntersectionPoint) {
        self.base.append_point(pnt);
    }

    /// `Perform(L, P, Tol, E)` (`BRepClass_Intersector.cxx:330-441`).
    /// Intersects the line `l` (its parameter bound `p`, `RealLast` for an
    /// unbounded ray) with edge `e`, using tolerance `tol`.
    pub fn perform(&mut self, l: &GpLin2d, p: f64, tol: f64, e: &BRepClassEdge) {
        let mut tol_z = tol;
        let (Some(edge), Some(face)) = (e.edge(), e.face()) else {
            self.base.done = false;
            return;
        };
        let Some((c2d, deb, fin)) = curve_on_surface_range(edge, face) else {
            self.base.done = false;
            return;
        };

        let mut bond = BndBox2d::new();
        let pnt_f = l.location();
        if e.use_bnd_box() {
            add_geom2d_range(c2d.as_ref(), deb, fin, 0.0, &mut bond);
            bond.set_gap(tol_z);
        }

        let c = adaptor_curve(&c2d, deb, fin);

        // "ON" case: direct check of belonging to the edge, with tolerance.
        if !e.use_bnd_box() || !bond.is_out(&pnt_f) {
            let (mut deb_tol, mut fin_tol) = (deb, fin);
            if tol_z > CONFUSION {
                deb_tol = deb - tol_z;
                fin_tol = fin + tol_z;
            }
            let cur = adaptor_curve(&c2d, deb_tol, fin_tol);
            if let Some(pnt) = check_on(face, l, cur.as_ref(), &mut tol_z, fin, deb) {
                self.base.append_point(&pnt);
                self.base.done = true;
                return;
            }
        }

        if e.use_bnd_box() {
            // `TopExp::Vertices(EE, aVF, aVL)`: plain orientation.
            let (vf, vl) = edge_vertices(edge);
            tol_z = max_tol_2d_cur_edge(vf.as_ref(), vl.as_ref(), face, tol);
            bond.set_gap(tol_z);
            if !is_inter(&bond, l, p) {
                self.base.done = false;
                return;
            }
        }

        let pdeb = c.d0(deb);
        let pfin = c.d0(fin);
        let toldeb = 1.0e-5;
        let tolfin = 1.0e-5;

        let mut dl = IntRes2dDomain::new();
        if p != f64::MAX {
            dl.set_bounded(
                &l.location(),
                0.0,
                PCONFUSION,
                &lin_value(l, p),
                p,
                PCONFUSION,
            );
        } else {
            dl.set_semi_infinite(&l.location(), 0.0, PCONFUSION, true);
        }

        let mut de = IntRes2dDomain::bounded(&pdeb, deb, toldeb, &pfin, fin, tolfin);
        if c.is_periodic() {
            let first = c.first_parameter();
            de.set_equivalent_parameters(first, first + c.period());
        }

        // `Geom2d_Line(L)` with `Geom2dAdaptor_Curve(GL)`.
        let gl = Geom2dLine::new(*l.position());

        let mut inter = Geom2dIntGInter::new();
        inter.perform(&gl, &dl, c.as_ref(), &de, PCONFUSION, PINTERSECTION);

        // A miss may be a hit through a high-tolerance vertex.
        if inter.result().is_empty() {
            if let Some(skip) = check_skip(l, &gl, e, c2d.as_ref(), &dl, deb, fin, e.max_tolerance())
            {
                inter = skip;
            }
        }

        self.base.set_values(inter.result());
    }

    /// `LocalGeometry(E, U, Tang, Norm, C)` (`BRepClass_Intersector.cxx:445-476`).
    /// Returns the tangent, the normal and the curvature of the edge's
    /// pcurve at `u`.
    pub fn local_geometry(&self, e: &BRepClassEdge, u: f64) -> (GpDir2d, GpDir2d, f64) {
        let (Some(edge), Some(face)) = (e.edge(), e.face()) else {
            panic!("BRepClass_Intersector::LocalGeometry: null edge or face");
        };
        let (pcurve, fpar, lpar) = curve_on_surface_range(edge, face)
            .expect("BRepClass_Intersector::LocalGeometry: no CurveOnSurface");
        let mut prop = GeomLPropCLProps2d::new(pcurve.as_ref(), u, 2, PCONFUSION);

        let mut curv = 0.0;
        // OCCT leaves `Tang` uninitialised when the chord is degenerate; the
        // X axis is used here as a neutral placeholder.
        let mut tang = GpDir2d::new(1.0, 0.0).expect("unit X direction");
        if prop.is_tangent_defined() {
            tang = prop.tangent();
            curv = prop.curvature();
        } else if let Some(chord) = get_tangent_as_chord(pcurve.as_ref(), u, fpar, lpar) {
            tang = chord;
        }

        let norm = if curv > PCONFUSION && !Precision::is_infinite(curv) {
            prop.normal()
        } else {
            GpDir2d::new(tang.y(), -tang.x()).expect("unit normal")
        };
        (tang, norm, curv)
    }
}
