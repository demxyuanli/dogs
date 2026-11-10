//! `TopOpeBRepDS_HDataStructure` — the `Standard_Transient` wrapper around
//! `TopOpeBRepDS_DataStructure` that the boolean engine passes around.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_HDataStructure.hxx:30-160`,
//! `TopOpeBRepDS_HDataStructure.cxx:1-813`. Split in two files: this one holds
//! the queries, `AddAncestors`, `MakeCurve` and the parameter sorts;
//! [`store`] holds the `StoreInterference*` / `ClearStoreInterferences` family.

mod store;


use std::sync::Arc;

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::topexp::Explorer;

use super::curve::Curve;
use super::curve_iterator::CurveIterator;
use super::data_structure::DataStructure;
use super::interference::Interference;
use super::interference_iterator::InterferenceList;
use super::kind::{is_geometry, is_topology, kind_to_shape, Config, Kind};
use super::point::Point;
use super::point_iterator::PointIterator;
use super::surface::Surface;
use super::surface_iterator::SurfaceIterator;

/// `TopOpeBRepDS_HDataStructure` (`hxx:30-160`). The C++ class holds
/// `TopOpeBRepDS_DataStructure myDS` by value and is reference-counted as a
/// `handle`; the port keeps the same value semantics.
#[derive(Clone, Default)]
pub struct HDataStructure {
    /// `myDS`
    ds: DataStructure,
}

/// `FUN_HDS_data` (`cxx:40-54`): the `(GeometryType, Geometry, SupportType,
/// Support)` quadruple of an interference. A null handle leaves the outputs
/// untouched in OCCT (they were initialised by the caller to
/// `UNKNOWN`/`0`), so the port returns `None`.
pub fn fun_hds_data(i: Option<&Arc<std::sync::RwLock<Interference>>>) -> Option<(Kind, i32, Kind, i32)> {
    let i = i?.read().expect("interference lock");
    Some((i.geometry_type(), i.geometry(), i.support_type(), i.support()))
}

/// `FUN_HDS_FACESINTERFER(F1, F2, HDS)` (`cxx:56-90`): true when `F1` has an
/// edge/face interference whose support is `F2`.
pub fn fun_hds_faces_interfer(f1: &TopoShape, f2: &TopoShape, hds: &HDataStructure) -> bool {
    let ds = hds.ds();
    if !ds.has_shape(f1, true) {
        return false;
    }

    let i_f2 = ds.shape_index(f2, true);
    for h in ds.shape_interferences(f1, true) {
        let i = h.read().expect("interference lock");
        // interference face1/edge/face2
        let mut fef = true;
        fef = fef && i.geometry_type() == Kind::Edge;
        fef = fef && i.support_type() == Kind::Face;
        fef = fef && i.support() == i_f2;
        if fef {
            return true;
        }
    }
    false
}

impl HDataStructure {
    /// `TopOpeBRepDS_HDataStructure()` (`cxx:96`): the only member is the
    /// by-value `myDS`, default-constructed.
    pub fn new() -> Self {
        Self {
            ds: DataStructure::new(),
        }
    }

    /// `DS()` (`cxx:143-146`).
    pub fn ds(&self) -> &DataStructure {
        &self.ds
    }

    /// `ChangeDS()` (`cxx:150-153`).
    pub fn change_ds(&mut self) -> &mut DataStructure {
        &mut self.ds
    }

    /// `AddAncestors(S)` (`cxx:99-107`): the four `(T1, T2)` passes.
    pub fn add_ancestors(&mut self, s: &TopoShape) {
        self.add_ancestors_of(s, ShapeType::Wire, ShapeType::Edge);
        self.add_ancestors_of(s, ShapeType::Face, ShapeType::Wire);
        self.add_ancestors_of(s, ShapeType::Shell, ShapeType::Face);
        self.add_ancestors_of(s, ShapeType::Solid, ShapeType::Shell);
    }

    /// `AddAncestors(S, T1, T2)` (`cxx:109-129`): every sub-shape `ST1` of
    /// type `T1` that contains a sub-shape `ST2` of type `T2` already in the
    /// DS is added with `rankS = AncestorRank(S)`. The inner loop `break`s on
    /// the first hit, so a `ST1` is added at most once.
    pub fn add_ancestors_of(&mut self, s: &TopoShape, t1: ShapeType, t2: ShapeType) {
        let rank_s = self.ds.ancestor_rank_of(s);
        let mut to_add: Vec<TopoShape> = Vec::new();
        let mut e_t1 = Explorer::new(s, t1);
        while e_t1.more() {
            let st1 = e_t1.current().clone();
            let mut e_t2 = Explorer::new(&st1, t2);
            while e_t2.more() {
                let st2 = e_t2.current().clone();
                if self.ds.has_shape(&st2, true) {
                    to_add.push(st1.clone());
                    break;
                }
                e_t2.next();
            }
            e_t1.next();
        }
        for st1 in to_add {
            self.ds.add_shape_with_ancestor(&st1, rank_s);
        }
    }

    /// `ChkIntg()` (`cxx:133-139`): builds a `TopOpeBRepDS_Check`, runs
    /// `ChkIntg()` on it and discards the returned status, exactly as OCCT
    /// does. The status is returned here so a caller can still see it.
    pub fn chk_intg(&self) -> bool {
        super::check::Check::new(self).chk_intg()
    }

    /// `NbSurfaces()` (`cxx:157-160`).
    pub fn nb_surfaces(&self) -> i32 {
        self.ds.nb_surfaces()
    }

    /// `NbCurves()` (`cxx:164-167`).
    pub fn nb_curves(&self) -> i32 {
        self.ds.nb_curves()
    }

    /// `NbPoints()` (`cxx:171-174`).
    pub fn nb_points(&self) -> i32 {
        self.ds.nb_points()
    }

    /// `Surface(I)` (`cxx:178-181`).
    pub fn surface(&self, i: i32) -> &Surface {
        self.ds.surface(i)
    }

    /// `SurfaceCurves(I)` (`cxx:185-189`).
    pub fn surface_curves(&self, i: i32) -> CurveIterator {
        CurveIterator::new(self.ds.surface_interferences(i).to_vec())
    }

    /// `Curve(I)` (`cxx:193-196`).
    pub fn curve(&self, i: i32) -> &Curve {
        self.ds.curve(i)
    }

    /// `ChangeCurve(I)` (`cxx:200-203`).
    pub fn change_curve(&mut self, i: i32) -> &mut Curve {
        self.ds.change_curve(i)
    }

    /// `CurvePoints(I)` (`cxx:207-211`).
    pub fn curve_points(&self, i: i32) -> PointIterator {
        PointIterator::new(self.ds.curve_interferences(i).to_vec())
    }

    /// `Point(I)` (`cxx:215-218`).
    pub fn point(&self, i: i32) -> Option<&Point> {
        self.ds.point(i)
    }

    /// `HasGeometry(S)` (`cxx:222-226`).
    pub fn has_geometry(&self, s: &TopoShape) -> bool {
        self.ds.has_geometry(s)
    }

    /// `HasShape(S, FindKeep)` (`cxx:230-233`).
    pub fn has_shape(&self, s: &TopoShape, find_keep: bool) -> bool {
        self.ds.has_shape(s, find_keep)
    }

    /// `HasSameDomain(S, FindKeep)` (`cxx:238-248`).
    pub fn has_same_domain(&self, s: &TopoShape, find_keep: bool) -> bool {
        if !self.has_shape(s, find_keep) {
            return false;
        }
        !self.ds.shape_same_domain(s).is_empty()
    }

    /// `SameDomain(S)` (`cxx:251-257`): the same-domain list.
    pub fn same_domain(&self, s: &TopoShape) -> &[TopoShape] {
        self.ds.shape_same_domain(s)
    }

    /// `SameDomainOrientation(S)` (`cxx:260-267`).
    pub fn same_domain_orientation(&self, s: &TopoShape) -> Config {
        if !self.has_shape(s, true) {
            return Config::UnshGeometry;
        }
        self.ds.same_domain_ori_of(s)
    }

    /// `SameDomainReference(S)` (`cxx:271-279`).
    pub fn same_domain_reference(&self, s: &TopoShape) -> i32 {
        if !self.has_shape(s, true) {
            return 0;
        }
        self.ds.same_domain_ref_of(s)
    }

    /// `NbShapes()` (`cxx:282-285`).
    pub fn nb_shapes(&self) -> i32 {
        self.ds.nb_shapes()
    }

    /// `Shape(I, FindKeep)` (`cxx:289-292`).
    pub fn shape(&self, i: i32, find_keep: bool) -> Option<TopoShape> {
        self.ds.shape(i, find_keep)
    }

    /// `Shape(S, FindKeep)` (`cxx:296-299`).
    pub fn shape_index(&self, s: &TopoShape, find_keep: bool) -> i32 {
        self.ds.shape_index(s, find_keep)
    }

    /// `SolidSurfaces(S)` (`cxx:303-307`).
    pub fn solid_surfaces_of_shape(&self, s: &TopoShape) -> SurfaceIterator {
        SurfaceIterator::new(self.ds.shape_interferences(s, true).to_vec())
    }

    /// `SolidSurfaces(I)` (`cxx:311-315`).
    pub fn solid_surfaces(&self, i: i32) -> SurfaceIterator {
        SurfaceIterator::new(self.ds.shape_interferences_at(i, true).to_vec())
    }

    /// `FaceCurves(F)` (`cxx:319-323`).
    pub fn face_curves_of_shape(&self, f: &TopoShape) -> CurveIterator {
        CurveIterator::new(self.ds.shape_interferences(f, true).to_vec())
    }

    /// `FaceCurves(I)` (`cxx:327-331`).
    pub fn face_curves(&self, i: i32) -> CurveIterator {
        CurveIterator::new(self.ds.shape_interferences_at(i, true).to_vec())
    }

    /// `EdgePoints(E)` (`cxx:335-339`).
    pub fn edge_points(&self, e: &TopoShape) -> PointIterator {
        PointIterator::new(self.ds.shape_interferences(e, true).to_vec())
    }

    /// `MakeCurve(curC, newC)` (`cxx:343-395`): clone the surface/curve
    /// interferences of `curC` onto the same shapes, register the new curve and
    /// point both SCIs at it.
    ///
    /// `SCI1` / `SCI2` are the *same* handles that are appended to the shape
    /// interference lists, so the geometry index written afterwards is visible
    /// through both (`cxx:380-393`); the port keeps the appended handle and
    /// updates it in place.
    pub fn make_curve(&mut self, cur_c: &Curve, new_c: &mut Curve) -> Result<i32, String> {
        let sci1 = cur_c
            .get_sci1()
            .map(|i| Arc::new(std::sync::RwLock::new(Interference::surface_curve_from_interference(
                &i.read().expect("interference lock"),
            ))));
        let sci2 = cur_c
            .get_sci2()
            .map(|i| Arc::new(std::sync::RwLock::new(Interference::surface_curve_from_interference(
                &i.read().expect("interference lock"),
            ))));

        let s1 = cur_c.shape1().cloned();
        let s2 = cur_c.shape2().cloned();

        // add the new SCIs in the lists of SCIs connected to the shapes
        if let (Some(sh), Some(i)) = (&s1, &sci1) {
            self.ds.add_shape_interference(sh, i.clone());
        }
        if let (Some(sh), Some(i)) = (&s2, &sci2) {
            self.ds.add_shape_interference(sh, i.clone());
        }

        // the shapes of the new curve are the shapes of curve curC
        new_c.set_shapes(s1, s2);

        // set surface/curve interferences of the new curve
        new_c.set_sci(sci1.clone(), sci2.clone());

        // newC is a daughter of curC
        new_c.change_mother(cur_c.ds_index());

        // add the new curve in the DS
        let inew_c = self.ds.add_curve(new_c);

        // the geometry of the new surface/curve interf is new curve inewC
        if let Some(i) = &sci1 {
            i.write().expect("interference lock").set_geometry(inew_c);
        }
        if let Some(i) = &sci2 {
            i.write().expect("interference lock").set_geometry(inew_c);
        }

        Ok(inew_c)
    }

    /// `RemoveCurve(icurC)` (`cxx:399-402`).
    pub fn remove_curve(&mut self, icur_c: i32) {
        self.ds.remove_curve(icur_c);
    }

    /// `NbGeometry(K)` (`cxx:406-431`): 0 for anything that is not a geometry
    /// kind (`TopOpeBRepDS::IsGeometry`).
    pub fn nb_geometry(&self, k: Kind) -> i32 {
        if !is_geometry(k) {
            return 0;
        }
        match k {
            Kind::Point => self.nb_points(),
            Kind::Curve => self.nb_curves(),
            Kind::Surface => self.nb_surfaces(),
            _ => 0,
        }
    }

    /// `NbTopology()` (`cxx:435-439`).
    pub fn nb_topology(&self) -> i32 {
        self.ds.nb_shapes()
    }

    /// `NbTopology(K)` (`cxx:443-464`): `FindKeep = false` here, unlike
    /// [`HDataStructure::nb_geometry`].
    pub fn nb_topology_of(&self, k: Kind) -> i32 {
        if !is_topology(k) {
            return 0;
        }
        let find_keep = false;
        let t = kind_to_shape(k);
        let mut res = 0;
        for i in 1..=self.nb_topology() {
            if self.ds.keep_shape(i, find_keep) {
                if let Some(sh) = self.ds.shape(i, find_keep) {
                    if sh.shape_type() == t {
                        res += 1;
                    }
                }
            }
        }
        res
    }

    /// `EdgesSameParameter()` (`cxx:467-483`).
    pub fn edges_same_parameter(&self) -> bool {
        for i in 1..=self.ds.nb_shapes() {
            let Some(s) = self.ds.shape(i, true) else {
                continue;
            };
            if s.shape_type() == ShapeType::Edge {
                let e = crate::shape::Edge(s);
                if !crate::brep_tool::BRepTool::same_parameter(&e) {
                    return false;
                }
            }
        }
        true
    }

    /// `SortOnParameter(List, SList)` (`cxx:533-573`).
    ///
    /// The sort itself is `FUN_TopOpeBRepDS_SortOnParameter` (`cxx:488-527`),
    /// an O(n^2) selection over `TopOpeBRepDS_PointIterator::Parameter()`;
    /// afterwards the list is split at the first interference whose
    /// `Transition().Orientation(IN)` is FORWARD, that element and everything
    /// after it going to the second list.
    pub fn sort_on_parameter(
        &self,
        list: &InterferenceList,
        s_list: &mut InterferenceList,
    ) -> Result<(), String> {
        fun_topopebrepds_sort_on_parameter(list, s_list)?;

        // tete = interf FORWARD
        let mut found = false;
        let mut l1: InterferenceList = Vec::new();
        let mut l2: InterferenceList = Vec::new();
        for i in s_list.iter() {
            if !found {
                let o = i
                    .read()
                    .expect("interference lock")
                    .transition()
                    .orientation(crate::abs::State::In, crate::abs::ShapeType::Face);
                if o == crate::abs::Orientation::Forward {
                    found = true;
                    l1.push(i.clone());
                } else {
                    l2.push(i.clone());
                }
            } else {
                l1.push(i.clone());
            }
        }

        s_list.clear();
        s_list.extend(l1);
        s_list.extend(l2);
        Ok(())
    }

    /// `SortOnParameter(List)` (`cxx:577-585`): only sorts when the list has
    /// at least one point interference (`PointIterator::More()`).
    pub fn sort_on_parameter_in_place(&self, list: &mut InterferenceList) -> Result<(), String> {
        let it = PointIterator::new(list.clone());
        if it.more() {
            let mut s_list = InterferenceList::new();
            self.sort_on_parameter(list, &mut s_list)?;
            *list = s_list;
        }
        Ok(())
    }

    /// `MinMaxOnParameter(List, parmin, parmax)` (`cxx:589-607`).
    ///
    /// An empty list leaves both outputs untouched in OCCT; the port returns
    /// `None` so the caller keeps its own values.
    pub fn min_max_on_parameter(&self, list: &InterferenceList) -> Result<Option<(f64, f64)>, String> {
        if list.is_empty() {
            return Ok(None);
        }
        let mut parmin = f64::MAX;
        let mut parmax = -f64::MAX;
        let mut it = PointIterator::new(list.clone());
        while it.more() {
            let parline = it.parameter()?;
            parmin = parmin.min(parline);
            parmax = parmax.max(parline);
            it.next();
        }
        Ok(Some((parmin, parmax)))
    }

    /// `ScanInterfList(IT, PDS)` (`cxx:618-645`): advances `it` to the first
    /// interference whose geometry is a point equal to `PDS` and reports
    /// whether it found one. The caller's iterator is left on that element, so
    /// the port takes `&mut`.
    pub fn scan_interf_list(&self, it: &mut PointIterator, pds: &Point) -> Result<bool, String> {
        while it.more() {
            let (gt, g) = {
                let i = it.current();
                (i.geometry_type(), i.geometry())
            };
            if gt == Kind::Point {
                if let Some(oopds) = self.ds.point(g) {
                    if pds.is_equal(oopds) {
                        return Ok(true);
                    }
                }
            } else if gt == Kind::Vertex {
                if let Some(sh) = self.ds.shape(g, true) {
                    let oopds = Point::of_shape(&sh)?;
                    if pds.is_equal(&oopds) {
                        return Ok(true);
                    }
                }
            }
            it.next();
        }
        Ok(false)
    }

    /// `GetGeometry(IT, PDS, G, K)` (`cxx:650-663`).
    pub fn get_geometry(
        &self,
        it: &mut PointIterator,
        pds: &Point,
    ) -> Result<Option<(i32, Kind)>, String> {
        let found = self.scan_interf_list(it, pds)?;
        if found {
            let i = it.current();
            return Ok(Some((i.geometry(), i.geometry_type())));
        }
        Ok(None)
    }
}

/// `FUN_TopOpeBRepDS_SortOnParameter` (`cxx:488-527`).
///
/// The `T` array is a `1..=n` `bool` mask; the outer loop runs `n` times and
/// each pass scans the whole list for the smallest still-unused parameter,
/// comparing with `<` (so the first of equal parameters wins).
///
/// OCCT keeps `Intf` across passes and appends the stale value when a pass
/// finds no candidate (a `NaN` parameter, or a parameter `>= RealLast()`);
/// `Err` is that degenerate case, which cannot be represented with
/// non-nullable handles here.
pub fn fun_topopebrepds_sort_on_parameter(
    list: &InterferenceList,
    s_list: &mut InterferenceList,
) -> Result<(), String> {
    let n_intf = list.len();
    if n_intf == 0 {
        return Ok(());
    }

    let mut used = vec![false; n_intf];

    for _ in 1..=n_intf {
        let mut parmin = f64::MAX;
        let mut i_intf = 0usize;
        let mut found = false;

        let mut it = PointIterator::new(list.clone());
        let mut itest = 1usize;
        while it.more() {
            if !used[itest - 1] {
                let par = it.parameter()?;
                if par < parmin {
                    parmin = par;
                    i_intf = itest;
                    found = true;
                }
            }
            it.next();
            itest += 1;
        }

        if !found {
            return Err("FUN_TopOpeBRepDS_SortOnParameter: no candidate".into());
        }
        s_list.push(list[i_intf - 1].clone());
        used[i_intf - 1] = true;
    }
    Ok(())
}
