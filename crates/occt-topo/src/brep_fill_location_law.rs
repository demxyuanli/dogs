//! `BRepFill_LocationLaw` — the law that places a section frame along a path.
//!
//! Source: `BRepFill_LocationLaw.cxx:48-744`, `.hxx:37-142`.
//!
//! OCCT declares this as a class with protected data that the `BRepFill_*Law`
//! subclasses (`BRepFill_Edge3DLaw`, `BRepFill_EdgeOnSurfLaw`,
//! `BRepFill_ACRLaw`) fill in through `Init` + `myEdges->SetValue` +
//! `myLaws->SetValue`. The port keeps that data in [`LocationLawBase`] (held by
//! each concrete law as a field) and exposes the two setters the subclasses
//! need ([`LocationLawBase::set_edge`], [`LocationLawBase::set_law`]).
//!
//! Translation notes:
//! - `NCollection_HArray1<T>` with OCCT's 1-based bounds becomes a `Vec` whose
//!   element `i - 1` is OCCT's `Value(i)`; the accessors keep the 1-based
//!   indices of the reference so call sites read the same.
//! - `occ::handle<GeomFill_LocationLaw>` is `Arc<dyn GeomFillLocationLaw>`; a
//!   null handle is `None`, and every place that would dereference a null
//!   handle in OCCT answers the `Err` of the function instead of crashing.
//! - `TopExp::FirstVertex(E)` / `LastVertex(E)` are called with OCCT's default
//!   `CumOri = false` (`TopExp.hxx:116`, `:121`), i.e. the vertices *stored*
//!   `FORWARD` / `REVERSED` inside the edge (`crate::topo_tools::edge_vertices`).
//! - `BRep_Tool::Pnt` / `Tolerance` are [`BRepTool::vertex_point`] /
//!   [`BRepTool::vertex_tolerance`]; `BRep_Builder::MakeVertex(P, Tol)` is
//!   [`crate::builder::TopoBuilder::make_vertex`].
//! - `BRepBuilderAPI_Transform(W, fila, true)` (`cxx:715`) is
//!   [`crate::shape_ops::transformed_copy`] — the copying variant.
//! - `Trsf.SetRotation(OZ, Angle)` (`cxx:194`) builds a rotation about the world
//!   Z axis with no translation; `SetTrsf` consumes only the 3x3 part
//!   (`GeomFill_CurveAndTrihedron::SetTrsf` stores the matrix), so the port
//!   builds that matrix directly with `gp_Mat::SetRotation`
//!   (`gp_Mat.cxx:122-158`).
//! - `myType` (`hxx:142`) is written by `TangentIsMain` / `NormalIsMain` /
//!   `BiNormalIsMain` and never read by `BRepFill_LocationLaw.cxx`; the field is
//!   kept for the member layout.
//! - `Parameter` leaves its `U` out-parameter untouched when no interval
//!   matches (`cxx:676-678`); the port answers `0.` there and `Index == 0`, the
//!   value OCCT does write.

use std::sync::Arc;

use occt_core::gp::{GpAx1, GpDir, GpMat, GpPnt, GpTrsf, GpVec, GpXyz};
use occt_geom::gcpnts::{abscissa_point_with_tolerance, curve_length_between};
use occt_geom::geom_fill::{GeomFillPipeError, LocationLaw as GeomFillLocationLaw};

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Vertex, Wire};
use crate::shape_ops::transformed_copy;
use crate::topo_tools::{edge_vertices, is_same, wire_vertices};

/// `EpsNul` (`cxx:500`) — a local constant of `IsG1`, not `gp::Resolution()`.
const EPS_NUL: f64 = 1.0e-12;

/// Angular tolerance default of `IsG1` (`hxx:103`).
const IS_G1_ANGULAR_TOLERANCE: f64 = 1.0e-4;

/// `Norm(const gp_Mat&)` (`cxx:46-67`): the max over rows of the sum of the
/// absolute values of the row.
fn norm(m: &GpMat) -> f64 {
    let mut max = 0.0_f64;
    for row in 1..=3 {
        let coord = m.row(row);
        let r = coord.x().abs() + coord.y().abs() + coord.z().abs();
        if r > max {
            max = r;
        }
    }
    max
}

/// `ToG0(M1, M2, T)` (`cxx:73-78`): `T = M2^-1 * M1`, the transformation that
/// takes the frame `M2` onto `M1`.
///
/// `BRepFill_DraftLaw.cxx:34-38` defines the same static function; the port
/// shares this one between both classes.
pub(crate) fn to_g0(m1: &GpMat, m2: &GpMat) -> Result<GpMat, String> {
    let inverted = m2
        .inverted()
        .map_err(|e| format!("gp_Mat::Inverted: {e}"))?;
    Ok(inverted.multiply(m1))
}

/// `BRepFill_LocationLaw` (`BRepFill_LocationLaw.hxx:37-142`).
pub struct LocationLawBase {
    /// `myPath` (`hxx:134`).
    path: Wire,
    /// `myTol` (`hxx:135`).
    tol: f64,
    /// `myLaws` (`hxx:136`).
    laws: Vec<Option<Arc<dyn GeomFillLocationLaw>>>,
    /// `myLength` (`hxx:137`). `RefCell` keeps `CurvilinearBounds` callable
    /// through the `Arc<dyn LocationLaw>` handle `BRepFill_Sweep` holds, which
    /// OCCT does with a plain mutable member.
    length: std::cell::RefCell<Vec<f64>>,
    /// `myEdges` (`hxx:138`).
    edges: Vec<Option<Edge>>,
    /// `myDisc` (`hxx:139`), `None` for a null handle. `RefCell` keeps
    /// `NbHoles` / `Holes` callable through the `Arc<dyn LocationLaw>` handle
    /// `BRepFill_Sweep` holds, which OCCT does with a plain mutable member.
    disc: std::cell::RefCell<Option<Vec<i32>>>,
    /// `myType` (`hxx:142`).
    law_type: i32,
}

impl Default for LocationLawBase {
    /// The compiler-generated default OCCT relies on before `Init` runs: a
    /// null path, `myTol == 0`, null arrays.
    fn default() -> Self {
        Self {
            path: Wire::new(),
            tol: 0.0,
            laws: Vec::new(),
            length: std::cell::RefCell::new(Vec::new()),
            edges: Vec::new(),
            disc: std::cell::RefCell::new(None),
            law_type: 0,
        }
    }
}

/// `BRepFill_LocationLaw` (`BRepFill_LocationLaw.hxx:37-142`) as a consumer
/// sees it.
///
/// OCCT's `BRepFill_Sweep` holds a `handle(BRepFill_LocationLaw)` and calls
/// only the base-class queries (`NbLaw`, `Law`, `Wire`, `IsClosed`, `IsG1`).
/// The port keeps the data in [`LocationLawBase`], which the concrete laws
/// embed as a field, so the trait exposes that base and leaves the concrete
/// types free to add their own data.
pub trait LocationLaw {
    /// The embedded `BRepFill_LocationLaw` data.
    fn base(&self) -> &LocationLawBase;

    /// The embedded `BRepFill_LocationLaw` data, mutably: the base-class
    /// mutators (`DeleteTransform`, `TransformInG0Law`,
    /// `TransformInCompatibleLaw`) reach the owning law through it.
    fn base_mut(&mut self) -> &mut LocationLawBase;
}

impl LocationLaw for LocationLawBase {
    fn base(&self) -> &LocationLawBase {
        self
    }

    fn base_mut(&mut self) -> &mut LocationLawBase {
        self
    }
}

impl LocationLawBase {

    pub fn new() -> Self {
        Self::default()
    }


    /// `Init(Path)` (`cxx:82-110`): count the non-degenerated edges of the
    /// path, size `myLaws` / `myLength` / `myEdges` from that count, and start
    /// in the tangent-main mode.
    pub fn init(&mut self, path: &Wire) {
        self.path = path.clone();
        self.tol = 1.0e-4;

        let mut nb_edge = 0usize;
        let mut wexp = WireExplorer::new();
        wexp.init(path, None);
        while wexp.more() {
            // `if (!BRep_Tool::Degenerated(wexp.Current())) NbEdge++;`
            match wexp.current_edge() {
                Some(e) if !BRepTool::is_degenerated(&e) => nb_edge += 1,
                _ => {}
            }
            wexp.next();
        }

        self.laws = vec![None; nb_edge];
        self.length = std::cell::RefCell::new(vec![-1.0; nb_edge + 1]);
        // `myLength->SetValue(1, 0.)` (`cxx:106`).
        self.length.borrow_mut()[0] = 0.0;
        self.edges = vec![None; nb_edge];
        self.disc = std::cell::RefCell::new(None);
        self.tangent_is_main();
    }

    /// `GetStatus` (`cxx:114-123`): the first non-`PipeOk` status of the
    /// elementary laws.
    pub fn get_status(&self) -> GeomFillPipeError {
        let mut status = GeomFillPipeError::PipeOk;
        for slot in &self.laws {
            if status != GeomFillPipeError::PipeOk {
                break;
            }
            if let Some(law) = slot {
                status = law.error_status();
            }
        }
        status
    }

    /// `TangentIsMain` (`cxx:127-130`).
    pub fn tangent_is_main(&mut self) {
        self.law_type = 1;
    }

    /// `NormalIsMain` (`cxx:132-135`).
    pub fn normal_is_main(&mut self) {
        self.law_type = 2;
    }

    /// `BiNormalIsMain` (`cxx:137-140`).
    pub fn bi_normal_is_main(&mut self) {
        self.law_type = 3;
    }

    /// `myType` (`hxx:142`) — written by the three modes above, never read by
    /// `BRepFill_LocationLaw.cxx`.
    pub fn law_type(&self) -> i32 {
        self.law_type
    }

    /// `TransformInCompatibleLaw(AngularTolerance)` (`cxx:149-198`): rotate
    /// each elementary law so the frames are continuous across the path's
    /// vertices — a G0 correction when the tangents match, otherwise the
    /// normal is carried over and only the rotation about the world Z axis is
    /// applied.
    pub fn transform_in_compatible_law(&mut self, tol_angular: f64) -> Result<(), String> {
        let nb = self.laws.len();
        let mut domain = self.law_at(1)?.get_domain();

        for ipath in 2..=nb {
            let (m1, _v) = self.law_at(ipath - 1)?.d0(domain.1)?;
            domain = self.law_at(ipath)?.get_domain();
            let (m2, _v) = self.law_at(ipath)?.d0(domain.0)?;

            let t1 = GpVec::from_xyz(&m1.column(3));
            let t2 = GpVec::from_xyz(&m2.column(3));
            let n1 = GpVec::from_xyz(&m1.column(1));
            let mut n2 = GpVec::from_xyz(&m2.column(1));

            let trsf = if t1.is_parallel_ang(&t2, tol_angular) && !t1.is_opposite_ang(&t2, tol_angular) {
                // Correction G0 (`cxx:169-171`).
                to_g0(&m1, &m2)?
            } else {
                let cross = t1.crossed(&t2);
                let alpha = t2.angle_with_ref(&t1, &cross);
                let dir = GpDir::new(cross.x(), cross.y(), cross.z())
                    .map_err(|e| format!("gp_Dir: {e}"))?;
                // `gp_Ax1 axe(gp::Origin(), cross.XYZ())` (`cxx:178`).
                let axe = GpAx1::new(GpPnt::zero(), dir);
                n2 = n2.rotated(&axe, alpha);
                let angle = n2.angle_with_ref(&n1, &t1);
                let mut rot = GpMat::zero();
                // `Trsf.SetRotation(gp_XYZ(0, 0, 1), Angle)` (`cxx:194`).
                rot.set_rotation(&GpXyz::new(0.0, 0.0, 1.0), angle)
                    .map_err(|e| format!("gp_Mat::SetRotation: {e}"))?;
                rot
            };
            self.law_at(ipath)?.set_trsf(&trsf);
        }
        Ok(())
    }

    /// `TransformInG0Law` (`cxx:202-227`): make every elementary law start on
    /// the previous law's end frame.
    pub fn transform_in_g0_law(&mut self) -> Result<(), String> {
        let nb = self.laws.len();
        let mut domain = self.law_at(1)?.get_domain();
        for ipath in 2..=nb {
            let (m1, _) = self.law_at(ipath - 1)?.d0(domain.1)?;
            domain = self.law_at(ipath)?.get_domain();
            let (m2, _) = self.law_at(ipath)?.d0(domain.0)?;
            let aux = to_g0(&m1, &m2)?;
            self.law_at(ipath)?.set_trsf(&aux);
        }

        // Is the law periodical ? (`cxx:219-226`): OCCT evaluates the two end
        // frames of a closed path and discards them — no `SetTrsf` follows — so
        // the calls are kept only for the control-flow shape of the reference.
        if self.path.is_closed() {
            let (m1, _) = self.law_at(nb)?.d0(domain.1)?;
            let first_domain = self.law_at(1)?.get_domain();
            let (m2, _) = self.law_at(1)?.d0(first_domain.0)?;
            let _ = (m1, m2);
        }
        Ok(())
    }

    /// `DeleteTransform` (`cxx:232-241`): reset every law's transformation and
    /// drop the cached discontinuities.
    pub fn delete_transform(&mut self) -> Result<(), String> {
        let identity = GpMat::identity();
        for ii in 1..=self.edges.len() {
            self.law_at(ii)?.set_trsf(&identity);
        }
        self.disc = std::cell::RefCell::new(None);
        Ok(())
    }

    /// `NbHoles(Tol)` (`cxx:247-277`): the number of G0 discontinuities of the
    /// law, computed once and cached in `myDisc`.
    pub fn nb_holes(&self, tol: f64) -> usize {
        if self.disc.borrow().is_none() {
            let mut seq: Vec<i32> = Vec::new();
            let nb = self.laws.len();
            let mut ii = 2usize;
            while ii <= nb + 1 {
                if self.is_g1(ii - 1, tol, EPS_NUL) == -1 {
                    seq.push(ii as i32);
                }
                ii += 1;
            }
            let nb_disc = seq.len();
            if nb_disc > 0 {
                *self.disc.borrow_mut() = Some(seq);
            }
        }
        match &*self.disc.borrow() {
            Some(d) => d.len(),
            None => 0,
        }
    }

    /// `Holes(Disc)` (`cxx:281-289`): copy the cached discontinuity indices.
    /// OCCT's `Disc` is a 1-based `NCollection_Array1`; the port writes
    /// `disc[i - 1]` for `Disc(i)`.
    pub fn holes(&self, disc: &mut [i32]) {
        if let Some(values) = &*self.disc.borrow() {
            for (i, v) in values.iter().enumerate() {
                if i < disc.len() {
                    disc[i] = *v;
                }
            }
        }
    }

    /// `NbLaw` (`cxx:292-295`).
    pub fn nb_law(&self) -> usize {
        self.laws.len()
    }

    /// `Law(Index)` (`cxx:299-301`), 1-based as in OCCT; `None` is a null
    /// handle.
    pub fn law(&self, index: usize) -> Option<Arc<dyn GeomFillLocationLaw>> {
        index
            .checked_sub(1)
            .and_then(|i| self.laws.get(i))
            .and_then(|slot| slot.clone())
    }

    /// `Wire()` (`cxx:307-309`).
    pub fn wire(&self) -> &Wire {
        &self.path
    }

    /// `Edge(Index)` (`cxx:313-316`).
    pub fn edge(&self, index: usize) -> Option<Edge> {
        index
            .checked_sub(1)
            .and_then(|i| self.edges.get(i))
            .and_then(|slot| slot.clone())
    }

    /// `Vertex(Index)` (`cxx:322-349`): the path vertex of rank `Index`, i.e.
    /// the start vertex of edge `Index`, the end vertex of the last edge when
    /// `Index == NbLaw() + 1`.
    ///
    /// `CumOri = false` on both `TopExp` calls, so the `REVERSED` branch picks
    /// the vertex stored `REVERSED` inside the edge. `Index == 0` would index
    /// `myEdges->Value(0)`, outside OCCT's `1..NbEdge` bounds; the port answers
    /// `None`.
    pub fn vertex(&self, index: usize) -> Option<Vertex> {
        let nb = self.edges.len();
        let edge = if index <= nb {
            self.edge(index)?
        } else if index == nb + 1 {
            self.edge(index - 1)?
        } else {
            return None;
        };
        let (first, last) = edge_vertices(&edge);
        if index <= nb {
            if edge.0.orientation() == Orientation::Reversed {
                last
            } else {
                first
            }
        } else if edge.0.orientation() == Orientation::Reversed {
            first
        } else {
            last
        }
    }

    /// `PerformVertex(Index, Input, TolMin, Output, Location)` (`cxx:358-453`):
    /// the swept image of a section vertex.
    ///
    /// `iloc == -1` uses the law before the vertex, `1` the law after it and
    /// `0` the barycentre of both.
    pub fn perform_vertex(
        &self,
        index: usize,
        input: &Vertex,
        tol_min: f64,
        iloc: i32,
    ) -> Result<Vertex, String> {
        let nb = self.laws.len();
        let is_bary = iloc == 0;
        let mut m1 = GpMat::zero();
        let mut v1 = GpVec::zero();
        let mut m2 = GpMat::zero();
        let mut v2 = GpVec::zero();

        if index > 0 && index < nb {
            if iloc <= 0 {
                let domain = self.law_at(index)?.get_domain();
                let (m, v) = self.law_at(index)?.d0(domain.1)?;
                m1 = m;
                v1 = v;
            }
            if iloc >= 0 {
                let domain = self.law_at(index + 1)?.get_domain();
                let (m, v) = self.law_at(index + 1)?.d0(domain.0)?;
                if iloc == 0 {
                    m2 = m;
                    v2 = v;
                } else {
                    m1 = m;
                    v1 = v;
                }
            }
        }

        if index == 0 || index == nb {
            if !self.path.is_closed() || self.is_g1(index, tol_min, IS_G1_ANGULAR_TOLERANCE) != 1 {
                let (m, v) = if index == 0 {
                    let domain = self.law_at(1)?.get_domain();
                    self.law_at(1)?.d0(domain.0)?
                } else {
                    let domain = self.law_at(nb)?.get_domain();
                    self.law_at(nb)?.d0(domain.1)?
                };
                m1 = m;
                v1 = v;
            } else {
                if iloc <= 0 {
                    let domain = self.law_at(nb)?.get_domain();
                    let (m, v) = self.law_at(nb)?.d0(domain.1)?;
                    m1 = m;
                    v1 = v;
                }
                if iloc >= 0 {
                    let domain = self.law_at(1)?.get_domain();
                    let (m, v) = self.law_at(1)?.d0(domain.0)?;
                    if iloc == 0 {
                        m2 = m;
                        v2 = v;
                    } else {
                        m1 = m;
                        v1 = v;
                    }
                }
            }
        }

        let p = BRepTool::vertex_point(input);
        if is_bary {
            let p1 = p.xyz().multiplied_mat(&m1).added(v1.xyz());
            let p2 = p.xyz().multiplied_mat(&m2).added(v2.xyz());
            let mut mid = GpXyz::new(0.0, 0.0, 0.0);
            mid.set_linear_form(0.5, &p1, 0.5, &p2);
            let tol = p1.subtracted(&p2).modulus() / 2.0 + tol_min;
            Ok(TopoBuilder::new().make_vertex(GpPnt::from_xyz(&mid), tol))
        } else {
            let moved = p.xyz().multiplied_mat(&m1).added(v1.xyz());
            Ok(TopoBuilder::new().make_vertex(GpPnt::from_xyz(&moved), tol_min))
        }
    }

    /// `CurvilinearBounds(Index, First, Last)` (`cxx:455-475`): the arc-length
    /// window of law `Index`, computing the whole table on first use.
    pub fn curvilinear_bounds(&self, index: usize) -> Result<(f64, f64), String> {
        let mut first = self.length_at(index)?;
        let mut last = self.length_at(index + 1)?;
        if last < 0.0 {
            let nb_e = self.edges.len();
            let mut length = 0.0_f64;
            for ii in 1..=nb_e {
                let law = self.law_at(ii)?;
                // `GetDomain(f, l)` (`cxx:468`) is evaluated for its side
                // effects only: the `Length` call below spans the whole curve.
                let _ = law.get_domain();
                let curve = law
                    .get_curve()
                    .ok_or_else(|| NULL_CURVE.to_string())?;
                length += curve_length_between(
                    &*curve,
                    curve.first_parameter(),
                    curve.last_parameter(),
                    Some(self.tol),
                );
                // `myLength->SetValue(ii + 1, Length)` (`cxx:470`).
                if let Some(slot) = self.length.borrow_mut().get_mut(ii) {
                    *slot = length;
                }
            }
            first = self.length_at(index)?;
            last = self.length_at(index + 1)?;
        }
        Ok((first, last))
    }

    /// `IsClosed()` (`cxx:477-488`): the path flag, or the two boundary
    /// vertices being the same.
    pub fn is_closed(&self) -> bool {
        if self.path.is_closed() {
            return true;
        }
        let (v1, v2) = wire_vertices(&self.path);
        match (v1, v2) {
            (Some(a), Some(b)) => is_same(&a.0, &b.0),
            // Two null `TopoDS_Shape`s are `IsSame` in OCCT.
            (None, None) => true,
            _ => false,
        }
    }

    /// `IsG1(Index, SpatialTolerance, AngularTolerance)` (`cxx:493-623`):
    /// `-1` not connected, `0` G0, `1` G1.
    pub fn is_g1(&self, index: usize, spatial_tolerance: f64, angular_tolerance: f64) -> i32 {
        let nb = self.laws.len();
        let mut tol_eps = spatial_tolerance;
        let mut ok_d1 = false;
        let mut m1 = GpMat::zero();
        let mut v1 = GpVec::zero();
        let mut dm1 = GpMat::zero();
        let mut dv1 = GpVec::zero();
        let mut m2 = GpMat::zero();
        let mut v2 = GpVec::zero();
        let mut dm2 = GpMat::zero();
        let mut dv2 = GpVec::zero();
        let mut edge: Option<Edge> = None;

        if index > 0 && index < nb {
            let Ok(law1) = self.law_at(index) else {
                return -1;
            };
            let domain = law1.get_domain();
            match law1.d1(domain.1) {
                Ok(r) => {
                    ok_d1 = true;
                    m1 = r.m;
                    v1 = r.v;
                    dm1 = r.dm;
                    dv1 = r.dv;
                }
                Err(_) => ok_d1 = false,
            }
            if !ok_d1 {
                let Ok((m, v)) = law1.d0(domain.1) else {
                    return -1;
                };
                m1 = m;
                v1 = v;
            }

            let Ok(law2) = self.law_at(index + 1) else {
                return -1;
            };
            let domain2 = law2.get_domain();
            if ok_d1 {
                match law2.d1(domain2.0) {
                    Ok(r) => {
                        m2 = r.m;
                        v2 = r.v;
                        dm2 = r.dm;
                        dv2 = r.dv;
                    }
                    Err(_) => ok_d1 = false,
                }
            }
            if !ok_d1 {
                let Ok((m, v)) = law2.d0(domain2.0) else {
                    return -1;
                };
                m2 = m;
                v2 = v;
            }

            // `E = TopoDS::Edge(myEdges->Value(Index + 1))` (`cxx:527`).
            edge = self.edge(index + 1);
        }

        if index == 0 || index == nb {
            if !self.path.is_closed() {
                return -1;
            }
            let Ok(law_last) = self.law_at(nb) else {
                return -1;
            };
            let domain = law_last.get_domain();
            match law_last.d1(domain.1) {
                Ok(r) => {
                    ok_d1 = true;
                    m1 = r.m;
                    v1 = r.v;
                    dm1 = r.dm;
                    dv1 = r.dv;
                }
                Err(_) => ok_d1 = false,
            }
            if !ok_d1 {
                let Ok((m, v)) = law_last.d0(domain.1) else {
                    return -1;
                };
                m1 = m;
                v1 = v;
            }

            let Ok(law_first) = self.law_at(1) else {
                return -1;
            };
            let domain_first = law_first.get_domain();
            if ok_d1 {
                match law_first.d1(domain_first.0) {
                    Ok(r) => {
                        m2 = r.m;
                        v2 = r.v;
                        dm2 = r.dm;
                        dv2 = r.dv;
                    }
                    Err(_) => ok_d1 = false,
                }
            }
            if !ok_d1 {
                let Ok((m, v)) = law_first.d0(domain_first.0) else {
                    return -1;
                };
                m2 = m;
                v2 = v;
            }

            // `E = TopoDS::Edge(myEdges->Value(1))` (`cxx:552`).
            edge = self.edge(1);
        }

        let Some(edge) = edge else {
            return -1;
        };
        let (first, last) = edge_vertices(&edge);
        let v = if edge.0.orientation() == Orientation::Reversed {
            last
        } else {
            first
        };
        let Some(v) = v else {
            return -1;
        };

        tol_eps += 2.0 * BRepTool::vertex_tolerance(&v);

        let mut is_g0 = true;
        let mut is_g1 = true;

        if v1.subtracted(&v2).magnitude() > tol_eps {
            is_g0 = false;
        }
        if norm(&m1.subtract(&m2)) > spatial_tolerance {
            is_g0 = false;
        }
        if !is_g0 {
            return -1;
        }
        if !ok_d1 {
            return 0; // No control of the derivative
        }

        if dv1.magnitude() > EPS_NUL
            && dv2.magnitude() > EPS_NUL
            && dv1.angle(&dv2) > angular_tolerance
        {
            is_g1 = false;
        }

        let norm1 = norm(&dm1);
        let norm2 = norm(&dm2);
        if norm1 > EPS_NUL || norm2 > EPS_NUL {
            if norm1 > EPS_NUL && norm2 > EPS_NUL {
                let ndm1 = dm1.multiply_scalar(1.0 / norm1);
                let ndm2 = dm2.multiply_scalar(1.0 / norm2);
                if norm(&ndm1.subtract(&ndm2)) > angular_tolerance {
                    is_g1 = false;
                }
            } else {
                is_g1 = false; // 1 Null the other is not
            }
        }

        if is_g1 {
            1
        } else {
            0
        }
    }

    /// `Parameter(Abcissa, Index, U)` (`cxx:625-680`): the elementary law and
    /// the parameter holding `Abcissa`.
    ///
    /// `Index == 0` is OCCT's "not found"; `U` is `0.` there because OCCT
    /// leaves its out-parameter untouched. `Err` is the `StdFail_NotDone` of
    /// `GCPnts_AbscissaPoint::Parameter()` (`cxx:672`).
    pub fn parameter(&self, abcissa: f64) -> Result<(usize, f64), String> {
        let nb_e = self.edges.len();

        // Control that the lengths are calculated (`cxx:631-635`).
        if self.length_at(nb_e + 1)? < 0.0 {
            self.curvilinear_bounds(nb_e)?;
        }

        let mut trouve = false;
        let mut iedge = 1usize;
        while iedge <= nb_e && !trouve {
            if self.length_at(iedge + 1)? >= abcissa {
                trouve = true;
            } else {
                iedge += 1;
            }
        }

        if !trouve {
            return Ok((0, 0.0));
        }

        let law = self.law_at(iedge)?;
        let (f, l) = law.get_domain();
        let length_hi = self.length_at(iedge + 1)?;
        let length_lo = self.length_at(iedge)?;
        if abcissa == length_hi {
            Ok((iedge, l))
        } else if abcissa == length_lo {
            Ok((iedge, f))
        } else {
            let curve = law.get_curve().ok_or_else(|| NULL_CURVE.to_string())?;
            let u = abscissa_point_with_tolerance(&*curve, abcissa - length_lo, f, self.tol)?;
            Ok((iedge, u))
        }
    }

    /// `D0(Abcissa, W)` (`cxx:688-724`): `W` moved onto the frame the law has
    /// at `Abcissa`. `None` is OCCT's `W.Nullify()` for an abscissa outside the
    /// path.
    pub fn d0(&self, abcissa: f64, w: &TopoShape) -> Result<Option<TopoShape>, String> {
        let (ind, u) = self.parameter(abcissa)?;
        if ind == 0 {
            return Ok(None);
        }
        let (m, v) = self.law_at(ind)?.d0(u)?;
        // `fila.SetValues(M(1,1) ... M(3,3), V.X() ... V.Z())` (`cxx:699-711`):
        // three columns of the frame, then the translation.
        let mut fila = GpTrsf::identity();
        fila.set_values(
            m.value(1, 1),
            m.value(1, 2),
            m.value(1, 3),
            v.x(),
            m.value(2, 1),
            m.value(2, 2),
            m.value(2, 3),
            v.y(),
            m.value(3, 1),
            m.value(3, 2),
            m.value(3, 3),
            v.z(),
        )
        .map_err(|e| format!("gp_Trsf::SetValues: {e}"))?;
        // `W = BRepBuilderAPI_Transform(W, fila, true)` (`cxx:715`).
        Ok(Some(transformed_copy(w, &fila)?))
    }

    /// `Abscissa(Index, Param)` (`cxx:729-744`): the arc length of `Param` on
    /// edge `Index`, measured from the start of the path.
    ///
    /// Read-only: the only state it can fill is the `myLength` table, which is
    /// already behind a `RefCell` through `CurvilinearBounds`.
    pub fn abscissa(&self, index: usize, param: f64) -> Result<f64, String> {
        let mut length = self.length_at(index)?;
        if length < 0.0 {
            let (_first, last) = self.curvilinear_bounds(index)?;
            length = last;
        }
        let law = self.law_at(index)?;
        let curve = law.get_curve().ok_or_else(|| NULL_CURVE.to_string())?;
        length += curve_length_between(
            &*curve,
            curve.first_parameter(),
            param,
            Some(self.tol),
        );
        Ok(length)
    }

    // ---- state the `BRepFill_*Law` subclasses write ----

    /// `myEdges->SetValue(Index, E)` as the subclasses do while walking the
    /// path (`BRepFill_Edge3DLaw.cxx:56`).
    pub fn set_edge(&mut self, index: usize, edge: &Edge) {
        if let Some(slot) = index.checked_sub(1).and_then(|i| self.edges.get_mut(i)) {
            *slot = Some(edge.clone());
        }
    }

    /// `myLaws->SetValue(Index, Law)` (`BRepFill_Edge3DLaw.cxx:62`).
    pub fn set_law(&mut self, index: usize, law: Arc<dyn GeomFillLocationLaw>) {
        if let Some(slot) = self.laws.get_mut(index.saturating_sub(1)) {
            *slot = Some(law);
        }
    }

    /// `myLaws.Nullify()` (`BRepFill_EdgeOnSurfLaw.cxx:76`): drop every
    /// elementary law.
    pub fn nullify_laws(&mut self) {
        self.laws = Vec::new();
    }

    /// Number of path edges, i.e. the OCCT `myEdges->Length()`.
    pub fn nb_edge(&self) -> usize {
        self.edges.len()
    }

    // ---- internals ----

    /// `myLaws->Value(Index)` with OCCT's 1-based index; `Err` where OCCT
    /// would dereference a null handle.
    fn law_at(&self, index: usize) -> Result<Arc<dyn GeomFillLocationLaw>, String> {
        self.law(index)
            .ok_or_else(|| "Standard_NullObject: BRepFill_LocationLaw::myLaws".to_string())
    }

    /// `myLength->Value(Index)` with OCCT's 1-based index.
    fn length_at(&self, index: usize) -> Result<f64, String> {
        let i = index
            .checked_sub(1)
            .ok_or_else(|| "Standard_OutOfRange: BRepFill_LocationLaw::myLength".to_string())?;
        self.length
            .borrow()
            .get(i)
            .copied()
            .ok_or_else(|| "Standard_OutOfRange: BRepFill_LocationLaw::myLength".to_string())
    }
}

const NULL_CURVE: &str = "Standard_NullObject: GeomFill_LocationLaw::GetCurve";
