//! `TopOpeBRepDS_GapFiller` — closes the small gaps left between the vertices
//! of a section and the section itself.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_GapFiller.hxx:30-106`,
//! `TopOpeBRepDS_GapFiller.cxx:43-537` (OCCT 8.0.0).
//!
//! `Perform` (`cxx:52-84`) has two halves:
//!
//! 1. **Pairing** (`cxx:58-79`): for every `POINT` interference of every curve,
//!    once per point index (`View`), `FindAssociatedPoints` collects the
//!    interferences on the same edge that face the same way, and
//!    `myAsso->Associate(I, ALI)` records the couple.
//! 2. **Rebuilding** (`cxx:81-84` → `BuildNewGeometries` `cxx:517-537`): every
//!    point index reachable from a curve is rebuilt once by `ReBuildGeom`
//!    (`cxx:465-514`) into the mid-parameter point of all the associated
//!    vertices, and every associated interference is moved onto it.
//!
//! OCCT reaches the edge/face connexity tables through the `FDSCNX_*` free
//! functions with their `thread_local` state; the port threads the prepared
//! [`Connexity`] through `perform` (`FDSCNX_Prepare` is called by
//! `TopOpeBRep_DSFiller` at `DSFiller.cxx:705`, before `Complete`).
//!
//! Failure modes that are `Option`/`return` here but raise in OCCT:
//!
//! * `myHDS->Point(G)` (`cxx:406`, `:429`, `:486`) is a `DataMap::Find` that
//!   throws on an unbound index; the port skips the interference.
//! * `BRepAdaptor_Curve C(E)` (`cxx:499`) on a null `E` (`EdgeSupport` failed
//!   at `cxx:480`) dereferences a null shape; [`GapFiller::rebuild_geom`]
//!   returns instead.

use occt_core::gp::{GpDir, GpPnt, GpVec};
use occt_core::precision::{INFINITE, RESOLUTION};

use crate::abs::ShapeType;
use crate::shape::{Face, TopoShape};
use crate::shape_upgrade::GEOMABS_C1;
use crate::topexp::Explorer;
use crate::topo_tools::is_same;
use crate::topopebrep::ds::association::Association;
use crate::topopebrep::ds::connex::Connexity;
use crate::topopebrep::ds::gap_tool::GapTool;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference::InterferenceHandle;
use crate::topopebrep::ds::interference_iterator::InterferenceList;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::point::Point;
use crate::topopebrep::tool::{CurveAdaptor, SurfaceAdaptor};
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;

/// `TopOpeBRepDS_GapFiller` (`hxx:30-106`).
pub struct GapFiller {
    /// `myGapTool` (`hxx:103`).
    gap_tool: GapTool,
    /// `myAsso` (`hxx:104`).
    asso: Association,
}

impl GapFiller {
    /// `TopOpeBRepDS_GapFiller(HDS)` (`cxx:44-50`): the constructor builds the
    /// tool (which runs `Init`) and an empty association.
    pub fn new(hds: &HDataStructure) -> Self {
        Self {
            gap_tool: GapTool::with_hds(hds),
            asso: Association::new(),
        }
    }

    /// `GapTool()` accessor, the port has no handle indirection (`hxx:47`).
    pub fn gap_tool(&self) -> &GapTool {
        &self.gap_tool
    }

    /// `Perform()` (`cxx:52-84`).
    pub fn perform(&mut self, hds: &mut HDataStructure, cx: &Connexity) {
        // `cxx:54`
        self.gap_tool.init(hds);
        // `cxx:55`: `NCollection_Map<int> View`.
        let mut view: Vec<i32> = Vec::new();

        // `cxx:58-79`: APPAIRAGE.
        let nb_curves = hds.ds().nb_curves();
        for i in 1..=nb_curves {
            let li = hds.ds().curve_interferences(i).to_vec();
            for handle in li {
                let g = match handle.read() {
                    Ok(itf) => itf.geometry_type() == Kind::Point,
                    Err(_) => false,
                };
                if !g {
                    continue;
                }
                let ip = match handle.read() {
                    Ok(itf) => itf.geometry(),
                    Err(_) => continue,
                };
                // `View.Add` answers false when the index is already in.
                if view.contains(&ip) {
                    continue;
                }
                view.push(ip);
                let ali = self.find_associated_points(&handle, hds, cx);
                self.asso.associate_list(&handle, &ali);
            }
        }

        // `cxx:81-84`
        self.build_new_geometries(hds);
    }

    /// `Contains(F, E)` (`cxx:87-100`): is `E` an occurrence of `F`?
    fn contains(f: &TopoShape, e: &TopoShape) -> bool {
        let mut ex = Explorer::new(f, e.shape_type());
        while ex.more() {
            if is_same(ex.current(), e) {
                return true;
            }
            ex.next();
        }
        false
    }

    /// `FindAssociatedPoints(I, LI)` (`cxx:103-199`): fills and answers `LI`.
    fn find_associated_points(
        &mut self,
        i: &InterferenceHandle,
        hds: &HDataStructure,
        cx: &Connexity,
    ) -> InterferenceList {
        let mut li: InterferenceList = Vec::new();

        // `cxx:108-115`
        for handle in self.gap_tool.same_interferences(i) {
            if self.asso.has_association(handle) {
                return li;
            }
        }

        // `cxx:122-131`: LI = the interferences on the edge.
        let Some(e) = self.gap_tool.edge_support(i) else {
            return li;
        };
        self.add_points_on_shape(&e, &mut li, hds);

        // `cxx:133-139`
        let Some((f1, f2)) = self.gap_tool.faces_support(i, hds) else {
            li.clear();
            return li;
        };
        // `cxx:141-145`
        let (f1, f2) = if !Self::contains(&f1, &e) {
            (f2.clone(), f1.clone())
        } else {
            (f1.clone(), f2.clone())
        };

        // `cxx:148`: `FDSCNX_EdgeConnexitySameShape(E, HDS)`.
        let lf = cx.edge_connexity_same_shape(&e, hds).to_vec();

        // `cxx:156-164`: a point lying on two faces connected through E means
        // the sections are already connected, so nothing is paired.
        for f in &lf {
            if !is_same(f, &f1) && self.is_on_face(i, &Face(f.clone()), hds) {
                li.clear();
                return li;
            }
        }

        // `cxx:171-177`
        for f in &lf {
            if !is_same(f, &f1) {
                self.filter_by_face(&Face(f.clone()), &mut li, hds);
            }
        }
        // `cxx:181-185`
        if !li.is_empty() {
            self.filter_by_face(&Face(f2.clone()), &mut li, hds);
        }
        // `cxx:189-192`
        if !li.is_empty() {
            self.filter_by_incident_distance(&Face(f2), i, &mut li, hds);
        }
        // `cxx:194-197`
        if !li.is_empty() {
            li.push(i.clone());
        }
        li
    }

    /// `CheckConnexity(LI)` (`cxx:205-210`): always `true` in OCCT.
    pub fn check_connexity(&self, _li: &mut InterferenceList) -> bool {
        true
    }

    /// `AddPointsOnShape(S, LI)` (`cxx:213-225`): append every interference the
    /// DS stores on `S`.
    fn add_points_on_shape(&self, s: &TopoShape, li: &mut InterferenceList, hds: &HDataStructure) {
        for handle in hds.ds().shape_interferences(s, true) {
            li.push(handle.clone());
        }
    }

    /// `AddPointsOnConnexShape(F, LI)` (`cxx:231-235`): empty in OCCT, the body
    /// is commented out in the source.
    fn add_points_on_connex_shape(&self, _f: &TopoShape, _li: &InterferenceList) {}

    /// `FilterByFace(F, LI)` (`cxx:239-258`): keep only the interferences that
    /// have a representation on `F`.
    fn filter_by_face(&self, f: &Face, li: &mut InterferenceList, hds: &HDataStructure) {
        li.retain(|handle| self.is_on_face(handle, f, hds));
    }

    /// `IsOnFace(I, F)` (`cxx:262-282`).
    fn is_on_face(&self, i: &InterferenceHandle, f: &Face, hds: &HDataStructure) -> bool {
        let Some(c) = self.gap_tool.curve(i, hds) else {
            return false;
        };
        let (s1, s2) = c.get_shapes();
        // `S1.IsSame(F)` / `S2.IsSame(F)`: a null shape is never same as a face.
        if let Some(s1) = s1 {
            if is_same(s1, &f.0) {
                return true;
            }
        }
        if let Some(s2) = s2 {
            if is_same(s2, &f.0) {
                return true;
            }
        }
        false
    }

    /// `FilterByEdge(E, LI)` (`cxx:287-303`): keep only the interferences whose
    /// support is `E`.
    fn filter_by_edge(&self, e: &TopoShape, li: &mut InterferenceList, hds: &HDataStructure) {
        li.retain(|handle| self.is_on_edge(handle, e, hds));
    }

    /// `IsOnEdge(I, E)` (`cxx:306-326`).
    fn is_on_edge(&self, i: &InterferenceHandle, e: &TopoShape, hds: &HDataStructure) -> bool {
        for handle in self.gap_tool.same_interferences(i) {
            let Ok(itf) = handle.read() else {
                continue;
            };
            if itf.support_type() != Kind::Edge {
                continue;
            }
            if let Some(s1) = hds.ds().shape(itf.support(), true) {
                if is_same(&s1, e) {
                    return true;
                }
            }
        }
        false
    }

    /// `FilterByIncidentDistance(F, I, LI)` (`cxx:398-462`): keep the single
    /// closest candidate within `20 * tolerance`.
    fn filter_by_incident_distance(
        &self,
        f: &Face,
        i: &InterferenceHandle,
        li: &mut InterferenceList,
        hds: &HDataStructure,
    ) {
        let mut dist_min = INFINITE;
        let mut i_sol: Option<InterferenceHandle> = None;

        // `cxx:406-408`
        let Some(pi1) = point_of(i, hds) else {
            li.clear();
            return;
        };
        let gpi: &GpPnt = pi1.point();

        // `cxx:413-416`
        let tol_def = 0.94; // cos(20 degrees)
        let tol_dist = 20.0 * pi1.tolerance();
        let n1 = normal(&self.gap_tool, hds, i, f);
        let ok1 = n1.is_some();

        // `cxx:418-455`
        for handle in li.iter() {
            // `cxx:422-425`
            if let (Ok(ci), Ok(ii)) = (handle.read(), i.read()) {
                if ci.has_same_geometry(&ii) {
                    continue;
                }
            }
            let n2 = normal(&self.gap_tool, hds, handle, f);
            let ok2 = n2.is_some();
            let Some(p) = point_of(handle, hds) else {
                continue;
            };
            let gp: &GpPnt = p.point();
            let dist = gp.distance(gpi);

            // `cxx:437-448`: both arms `continue`, so the normal test of
            // `cxx:440-446` only decides which `continue` is taken. The
            // structure is kept because the OCCT code keeps it.
            if dist > tol_dist {
                if ok1 && ok2 {
                    if let (Some(n1), Some(n2)) = (n1.as_ref(), n2.as_ref()) {
                        if n1.dot(n2) < tol_def {
                            continue; // Reject
                        }
                    }
                }
                continue; // Reject
            }

            // `cxx:450-454`
            if dist < dist_min {
                dist_min = dist;
                i_sol = Some(handle.clone());
            }
        }

        // `cxx:457-461`
        li.clear();
        if let Some(sol) = i_sol {
            li.push(sol);
        }
    }

    /// `ReBuildGeom(I, View)` (`cxx:465-514`).
    fn rebuild_geom(
        &mut self,
        i: &InterferenceHandle,
        view: &mut Vec<i32>,
        hds: &mut HDataStructure,
    ) {
        if !self.asso.has_association(i) {
            return;
        }

        // `cxx:475`: `LI = myAsso->Associated(I)` (a writable reference in
        // OCCT; the list is snapshotted because `SetPoint` below appends to the
        // gap tool's index, not to this list).
        let li = match self.asso.associated(i) {
            Some(l) => l.clone(),
            None => return,
        };

        // `cxx:478-480`
        let mut tol_max = 0.0_f64;
        let mut u_min = INFINITE;
        let mut u_max = -u_min;
        let Some(e) = self.gap_tool.edge_support(i) else {
            // UNPORTED-adjacent guard: OCCT carries a null `E` into
            // `BRepAdaptor_Curve C(E)` (`cxx:499`), which raises.
            return;
        };

        // `cxx:482-497`: the new point.
        for handle in &li {
            let Some(pp) = point_of(handle, hds) else {
                continue;
            };
            tol_max = tol_max.max(pp.tolerance());
            if let Some(u) = self.gap_tool.parameter_on_edge(handle, &e, hds) {
                u_min = u_min.min(u);
                u_max = u_max.max(u);
            }
            // `cxx:493-496`
            match self.gap_tool.edge_support(handle) {
                Some(ce) => {
                    if !is_same(&ce, &e) {
                        return;
                    }
                }
                None => return,
            }
        }
        let u = (u_max + u_min) * 0.5;
        // `cxx:499-501`
        let adaptor = CurveAdaptor::of_edge(&crate::shape::Edge(e.clone()));
        let Some(gp) = adaptor.d0(u) else {
            return;
        };
        let p = Point::with_tolerance(gp, tol_max);

        // `cxx:504-512`
        let ip = hds.change_ds().add_point(&p);
        for handle in &li {
            if let Ok(itf) = handle.read() {
                view.push(itf.geometry());
            }
            self.gap_tool.set_parameter_on_edge(handle, &e, u);
            self.gap_tool.set_point(handle, ip);
        }
        self.gap_tool.set_parameter_on_edge(i, &e, u);
        self.gap_tool.set_point(i, ip);
        view.push(ip);
    }

    /// `BuildNewGeometries()` (`cxx:517-537`).
    fn build_new_geometries(&mut self, hds: &mut HDataStructure) {
        let mut view: Vec<i32> = Vec::new();
        let nb_curves = hds.ds().nb_curves();
        let nb_points = hds.ds().nb_points();
        for i in 1..=nb_curves {
            let li = hds.ds().curve_interferences(i).to_vec();
            for handle in li {
                let ip = match handle.read() {
                    Ok(itf) => itf.geometry(),
                    Err(_) => continue,
                };
                // `cxx:531-534`
                if !view.contains(&ip) && ip <= nb_points {
                    view.push(ip);
                    self.rebuild_geom(&handle, &mut view, hds);
                }
            }
        }
    }
}

/// `Normal(A, HDS, I, F, D)` (`cxx:330-394`): the normal to `F` at the point of
/// `I`. `None` is OCCT's `false`.
fn normal(
    a: &GapTool,
    hds: &HDataStructure,
    i: &InterferenceHandle,
    f: &Face,
) -> Option<GpDir> {
    // `cxx:343-348`
    let s = SurfaceAdaptor::of_face(f);
    if s.get_type() == GeomAbsSurfaceType::Plane {
        let pln = s.surface()?.gp_pln()?;
        return Some(*pln.axis().direction());
    }

    // `cxx:350-392`
    for handle in a.same_interferences(i) {
        let Ok(itf) = handle.read() else {
            continue;
        };
        if itf.support_type() != Kind::Curve {
            continue;
        }
        // `cxx:358-359`
        if !itf.is_curve_point() {
            continue;
        }
        let Some(par) = itf.parameter() else {
            continue;
        };
        let c = hds.ds().curve(itf.support());
        let (s1, s2) = c.get_shapes();
        // `cxx:362-378`
        let p2d = match (s1, s2) {
            (Some(s1), _) if is_same(s1, &f.0) => {
                let c1 = c.curve1()?;
                c1.value(par)
            }
            _ => {
                let c2 = c.curve2()?;
                c2.value(par)
            }
        };
        // `cxx:380-391`: `UContinuity() >= C1 && VContinuity() >= C1`. The
        // port's `Surface::continuity` is the surface-wide value, i.e. the
        // minimum over both directions, so the conjunction is one test here.
        let surf = s.surface()?;
        if surf.continuity() >= GEOMABS_C1 {
            let (_, tu, tv) = surf.d1(p2d.x(), p2d.y());
            let n: GpVec = tu.crossed(&tv);
            if n.square_magnitude() >= RESOLUTION {
                return GpDir::from_vec(&n).ok();
            }
        }
        return None;
    }
    None
}

/// `myHDS->Point(G)` for the geometry of `I` (`cxx:406`, `:429`, `:486`).
fn point_of<'a>(i: &InterferenceHandle, hds: &'a HDataStructure) -> Option<&'a Point> {
    let g = i.read().ok()?.geometry();
    hds.ds().point(g)
}
