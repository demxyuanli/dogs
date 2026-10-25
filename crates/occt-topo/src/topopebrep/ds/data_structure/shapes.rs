//! `TopOpeBRepDS_DataStructure` — the shape half: the `myShapes` indexed
//! data map, the section edges, and every accessor of the same-domain,
//! interference and ancestor-rank records.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_DataStructure.cxx:246-1001`.
//!
//! `TopOpeBRepDS_ShapeData`'s fields are private in OCCT and written here
//! through the `friend class` declared at `TopOpeBRepDS_DataStructure.hxx:311`;
//! the Rust fields are `pub(super)` to `topopebrep::ds` and reached through
//! the accessors of [`crate::topopebrep::ds::ShapeData`].

use crate::shape::TopoShape;

use super::super::interference::InterferenceHandle;
use super::super::kind::Config;
use super::super::shape_data::ShapeData;
use super::DataStructure;

impl DataStructure {
    // --- shapes -----------------------------------------------------------

    /// `AddShape(S)` (`cxx:246-258`).
    pub fn add_shape(&mut self, s: &TopoShape) -> i32 {
        let mut is = self.shapes.find_index(s);
        if is == 0 {
            is = self.shapes.add(s, ShapeData::new());
            // a shape is its own reference, oriented as itself
            self.set_same_domain_ref(is, is);
            self.set_same_domain_ori(is, Config::SameOriented);
        }
        is
    }

    /// `AddShape(S, Ianc)` (`cxx:262-275`).
    pub fn add_shape_with_ancestor(&mut self, s: &TopoShape, ianc: i32) -> i32 {
        let mut is = self.shapes.find_index(s);
        if is == 0 {
            is = self.shapes.add(s, ShapeData::new());
            // a shape is its own reference, oriented as itself
            self.set_same_domain_ref(is, is);
            self.set_same_domain_ori(is, Config::SameOriented);
            self.set_ancestor_rank(is, ianc);
        }
        is
    }

    /// `ChangeShapes()` (`cxx:360-364`).
    pub fn change_shapes(&mut self) -> &mut super::IndexedShapeData<ShapeData> {
        &mut self.shapes
    }

    /// `NbShapes()` (`cxx:1036-1041`).
    pub fn nb_shapes(&self) -> i32 {
        self.shapes.extent()
    }

    /// `KeepShape(I, FindKeep)` (`cxx:279-293`).
    pub fn keep_shape(&self, i: i32, find_keep: bool) -> bool {
        if i >= 1 && i <= self.shapes.extent() {
            if find_keep {
                return self
                    .shapes
                    .find_from_index(i)
                    .map(|sd| sd.keep())
                    .unwrap_or(false);
            }
            return true;
        }
        false
    }

    /// `KeepShape(S, FindKeep)` (`cxx:299-313`).
    pub fn keep_shape_of(&self, s: &TopoShape, find_keep: bool) -> bool {
        if self.shapes.is_empty() || s.is_null() || !self.shapes.contains(s) {
            return false;
        }
        if find_keep {
            return self
                .shapes
                .find_from_key(s)
                .map(|sd| sd.keep())
                .unwrap_or(false);
        }
        true
    }

    /// `ChangeKeepShape(I, FindKeep)` (`cxx:319-326`).
    pub fn change_keep_shape(&mut self, i: i32, find_keep: bool) {
        if i >= 1 && i <= self.shapes.extent() {
            if let Some(s) = self.shapes.find_key(i).cloned() {
                self.change_keep_shape_of(&s, find_keep);
            }
        }
    }

    /// `ChangeKeepShape(S, FindKeep)` (`cxx:330-337`).
    pub fn change_keep_shape_of(&mut self, s: &TopoShape, find_keep: bool) {
        if !s.is_null() && self.shapes.contains(s) {
            if let Some(sd) = self.shapes.change_from_key(s) {
                sd.change_keep(find_keep);
            }
        }
    }

    /// `Shape(I, FindKeep)` (`cxx:1147-1154`): `myEmptyShape` (a null shape)
    /// is `None`.
    pub fn shape(&self, i: i32, find_keep: bool) -> Option<TopoShape> {
        if self.keep_shape(i, find_keep) {
            return self.shapes.find_key(i).cloned();
        }
        None
    }

    /// `Shape(S, FindKeep)` (`cxx:1159-1167`): the index, 0 when absent.
    pub fn shape_index(&self, s: &TopoShape, find_keep: bool) -> i32 {
        if self.has_shape(s, find_keep) {
            return self.shapes.find_index(s);
        }
        0
    }

    /// `HasShape(S, FindKeep)` (`cxx:1220-1224`).
    pub fn has_shape(&self, s: &TopoShape, find_keep: bool) -> bool {
        self.keep_shape_of(s, find_keep)
    }

    /// `HasGeometry(S)` (`cxx:1208-1216`).
    pub fn has_geometry(&self, s: &TopoShape) -> bool {
        self.has_shape(s, true) && !self.shape_interferences(s, true).is_empty()
    }

    // --- section edges ----------------------------------------------------

    /// `InitSectionEdges()` (`cxx:341-345`).
    pub fn init_section_edges(&mut self) {
        self.section_edges.clear();
    }

    /// `AddSectionEdge(E)` (`cxx:348-355`).
    pub fn add_section_edge(&mut self, e: &TopoShape) -> i32 {
        let mut ie = self.section_edges.find_index(e);
        if ie == 0 {
            ie = self.section_edges.add(e);
        }
        ie
    }

    /// `NbSectionEdges()` (`cxx:1044-1049`).
    pub fn nb_section_edges(&self) -> i32 {
        self.section_edges.extent()
    }

    /// `SectionEdge(I, FindKeep)` (`cxx:1172-1180`).
    ///
    /// OCCT returns `TopoDS::Edge(myEmptyShape)` when the shape is absent;
    /// that is `None` here.
    pub fn section_edge(&self, i: i32, find_keep: bool) -> Option<TopoShape> {
        let s = self.section_edges.find_key(i)?;
        if self.has_shape(s, find_keep) {
            return Some(s.clone());
        }
        None
    }

    /// `SectionEdge(E, FindKeep)` (`cxx:1184-1193`): the index, 0 when absent.
    pub fn section_edge_index(&self, e: &TopoShape, find_keep: bool) -> i32 {
        if self.keep_shape_of(e, find_keep) {
            return self.section_edges.find_index(e);
        }
        0
    }

    /// `IsSectionEdge(E, FindKeep)` (`cxx:1196-1204`).
    pub fn is_section_edge(&self, e: &TopoShape, find_keep: bool) -> bool {
        if self.keep_shape_of(e, find_keep) {
            return self.section_edges.contains(e);
        }
        false
    }

    // --- interference lists ----------------------------------------------

    /// `SurfaceInterferences(I)` (`cxx:368-379`).
    pub fn surface_interferences(&self, i: i32) -> &[InterferenceHandle] {
        match self.surfaces.get(&i) {
            Some(sd) => sd.geometry().interferences(),
            None => &self.empty_interferences,
        }
    }

    /// `ChangeSurfaceInterferences(I)` (`cxx:383-395`).
    pub fn change_surface_interferences(&mut self, i: i32) -> &mut Vec<InterferenceHandle> {
        match self.surfaces.get_mut(&i) {
            Some(sd) => sd.geometry_mut().change_interferences(),
            None => &mut self.empty_interferences,
        }
    }

    /// `CurveInterferences(I)` (`cxx:398-409`).
    pub fn curve_interferences(&self, i: i32) -> &[InterferenceHandle] {
        match self.curves.get(&i) {
            Some(cd) => cd.geometry().interferences(),
            None => &self.empty_interferences,
        }
    }

    /// `ChangeCurveInterferences(I)` (`cxx:413-424`).
    pub fn change_curve_interferences(&mut self, i: i32) -> &mut Vec<InterferenceHandle> {
        match self.curves.get_mut(&i) {
            Some(cd) => cd.geometry_mut().change_interferences(),
            None => &mut self.empty_interferences,
        }
    }

    /// `PointInterferences(I)` (`cxx:428-439`).
    pub fn point_interferences(&self, i: i32) -> &[InterferenceHandle] {
        match self.points.get(&i) {
            Some(pd) => pd.geometry().interferences(),
            None => &self.empty_interferences,
        }
    }

    /// `ChangePointInterferences(I)` (`cxx:443-453`).
    pub fn change_point_interferences(&mut self, i: i32) -> &mut Vec<InterferenceHandle> {
        match self.points.get_mut(&i) {
            Some(pd) => pd.geometry_mut().change_interferences(),
            None => &mut self.empty_interferences,
        }
    }

    /// `ShapeInterferences(S, FindKeep)` (`cxx:458-466`).
    pub fn shape_interferences(&self, s: &TopoShape, find_keep: bool) -> &[InterferenceHandle] {
        if self.has_shape(s, find_keep) {
            if let Some(sd) = self.shapes.find_from_key(s) {
                return sd.interferences();
            }
        }
        &self.empty_interferences
    }

    /// `ChangeShapeInterferences(S)` (`cxx:470-478`): appends to a
    /// throw-away list when `S` is not in `myShapes`.
    pub fn change_shape_interferences(&mut self, s: &TopoShape) -> &mut Vec<InterferenceHandle> {
        if !self.has_shape(s, true) {
            return &mut self.empty_interferences;
        }
        match self.shapes.change_from_key(s) {
            Some(sd) => sd.change_interferences(),
            None => &mut self.empty_interferences,
        }
    }

    /// `ShapeInterferences(I, FindKeep)` (`cxx:483-490`).
    pub fn shape_interferences_at(&self, i: i32, find_keep: bool) -> &[InterferenceHandle] {
        if find_keep && !self.keep_shape(i, true) {
            return &self.empty_interferences;
        }
        match self.shapes.find_from_index(i) {
            Some(sd) => sd.interferences(),
            None => &self.empty_interferences,
        }
    }

    /// `ChangeShapeInterferences(I)` (`cxx:495-500`).
    pub fn change_shape_interferences_at(&mut self, i: i32) -> &mut Vec<InterferenceHandle> {
        match self.shapes.change_from_index(i) {
            Some(sd) => sd.change_interferences(),
            None => &mut self.empty_interferences,
        }
    }

    /// `AddShapeInterference(S, I)` (`cxx:804-809`).
    pub fn add_shape_interference(&mut self, s: &TopoShape, i: InterferenceHandle) {
        self.change_shape_interferences(s).push(i);
    }

    /// `RemoveShapeInterference(S, I)` (`cxx:813-823`).
    ///
    /// `FindInterference(IT, I)` (`cxx:1255-1266`) advances the list iterator
    /// to the first handle equal to `I`; the found entry is then removed.
    pub fn remove_shape_interference(&mut self, s: &TopoShape, i: &InterferenceHandle) {
        let list = self.change_shape_interferences(s);
        if let Some(pos) = list.iter().position(|it| std::sync::Arc::ptr_eq(it, i)) {
            list.remove(pos);
        }
    }

    // --- same domain ------------------------------------------------------

    /// `ShapeSameDomain(S)` (`cxx:504-516`).
    pub fn shape_same_domain(&self, s: &TopoShape) -> &[TopoShape] {
        if !s.is_null() {
            if let Some(sd) = self.shapes.find_from_key(s) {
                return sd.same_domain();
            }
        }
        &[]
    }

    /// `ChangeShapeSameDomain(S)` (`cxx:521-526`).
    ///
    /// OCCT calls `ChangeFromKey`, which throws `Standard_NoSuchObject` when
    /// `S` is absent; that is `None` here.
    pub fn change_shape_same_domain(&mut self, s: &TopoShape) -> Option<&mut Vec<TopoShape>> {
        self.shapes.change_from_key(s).map(|sd| sd.change_same_domain())
    }

    /// `ShapeSameDomain(I)` (`cxx:530-542`).
    pub fn shape_same_domain_at(&self, i: i32) -> &[TopoShape] {
        if i >= 1 && i <= self.shapes.extent() {
            if let Some(sd) = self.shapes.find_from_index(i) {
                return sd.same_domain();
            }
        }
        &[]
    }

    /// `ChangeShapeSameDomain(I)` (`cxx:546-550`).
    pub fn change_shape_same_domain_at(&mut self, i: i32) -> Option<&mut Vec<TopoShape>> {
        self.shapes
            .change_from_index(i)
            .map(|sd| sd.change_same_domain())
    }

    /// `AddShapeSameDomain(S, SSD)` (`cxx:554-572`).
    pub fn add_shape_same_domain(&mut self, s: &TopoShape, ssd: &TopoShape) {
        let append = !self
            .shape_same_domain(s)
            .iter()
            .any(|it_s| crate::topo_tools::is_same(it_s, ssd));
        if append {
            if let Some(list) = self.change_shape_same_domain(s) {
                list.push(ssd.clone());
            }
        }
    }

    /// `RemoveShapeSameDomain(S, SSD)` (`cxx:577-595`).
    pub fn remove_shape_same_domain(&mut self, s: &TopoShape, ssd: &TopoShape) {
        if let Some(list) = self.change_shape_same_domain(s) {
            list.retain(|it_s| !crate::topo_tools::is_same(it_s, ssd));
        }
    }

    /// `UnfillShapesSameDomain(S1, S2)` (`cxx:996-1001`).
    pub fn unfill_shapes_same_domain(&mut self, s1: &TopoShape, s2: &TopoShape) {
        self.remove_shape_same_domain(s1, s2);
        self.remove_shape_same_domain(s2, s1);
    }

    /// `FillShapesSameDomain(S1, S2, refFirst)` (`cxx:828-951`).
    ///
    /// `Err` covers the two `throw Standard_ProgramError` sites: `r == 0`
    /// (`cxx:908-911`) and the "S1,S2 have SameDomain reference" mismatch,
    /// whose `throw` is commented out in OCCT (`cxx:899-903`) and therefore
    /// never fires here either.
    pub fn fill_shapes_same_domain(
        &mut self,
        s1: &TopoShape,
        s2: &TopoShape,
        ref_first: bool,
    ) -> Result<(), String> {
        let o1 = s1.orientation();
        let o2 = s2.orientation();

        let is1 = self.add_shape_with_ancestor(s1, 1);
        let sd1 = self
            .shapes
            .change_from_index(is1)
            .ok_or_else(|| "FacesFiller::Insert SD 1".to_string())?;
        let isdef1 = sd1.orientation_def();
        let mut todef1 = !isdef1;
        if isdef1 && sd1.same_domain_ori() == Config::UnshGeometry {
            todef1 = true;
        }
        if todef1 {
            sd1.set_orientation(o1);
            sd1.set_orientation_def(true);
        }

        let is2 = self.add_shape_with_ancestor(s2, 2);
        let sd2 = self
            .shapes
            .change_from_index(is2)
            .ok_or_else(|| "FacesFiller::Insert SD 1".to_string())?;
        let isdef2 = sd2.orientation_def();
        let mut todef2 = !isdef2;
        if isdef2 && sd2.same_domain_ori() == Config::UnshGeometry {
            todef2 = true;
        }
        if todef2 {
            sd2.set_orientation(o2);
            sd2.set_orientation_def(true);
        }

        let n1 = self.shape_same_domain(s1).len() as i32;
        let n2 = self.shape_same_domain(s2).len() as i32;

        self.add_shape_same_domain(s1, s2);
        self.add_shape_same_domain(s2, s1);

        let n11 = self.shape_same_domain(s1).len() as i32;
        let n22 = self.shape_same_domain(s2).len() as i32;

        let cond = (n11 == n1) && (n22 == n2) && !todef1 && !todef2;
        if cond {
            // nothing changed in SameDomain data of S1 and S2 : return
            return Ok(());
        }

        let r1 = self.same_domain_ref_of(s1);
        let r2 = self.same_domain_ref_of(s2);
        let r;
        // r1 == i1  r2 == i2 : shapes have no SameDomain ref : take S1 as reference
        // r1 == i1  r2 != i2 : S2 has a SameDomain reference : give it to S1
        // r1 != i1  r2 == i2 : S1 has a SameDomain reference : give it to S2
        // r1 != i1  r2 != i2 : S1,S2 have SameDomain reference : check equality
        if r1 == is1 && r2 == is2 {
            r = if ref_first { is1 } else { is2 };
        } else if r1 == is1 && r2 != is2 {
            r = r2;
        } else if r1 != is1 && r2 == is2 {
            r = r1;
        } else {
            // `if (r1 != r2) { // throw Standard_ProgramError(...) }` (`cxx:899-903`)
            r = if ref_first { r1 } else { r2 };
        }

        if r == 0 {
            return Err("FacesFiller::Insert SD 2".into());
        }

        let mut sr = self
            .shape(r, true)
            .ok_or_else(|| "FacesFiller::Insert SD 2".to_string())?;
        let o_sr = self
            .shapes
            .find_from_key(&sr)
            .map(|sd| sd.orientation())
            .ok_or_else(|| "FacesFiller::Insert SD 2".to_string())?;
        sr.set_orientation(o_sr);

        if r != r1 || todef1 {
            // S1 gets a new reference r
            let mut o = Config::SameOriented;
            if r != is1 || todef1 {
                let sso = crate::topopebrep::tool::shapes_same_oriented(s1, &sr);
                if !sso {
                    o = Config::DiffOriented;
                }
            }
            self.set_same_domain_ref(is1, r);
            self.set_same_domain_ori(is1, o);
        }

        if r != r2 || todef2 {
            // S2 gets a new reference r
            let mut o = Config::SameOriented;
            if r != is2 || todef2 {
                let sso = crate::topopebrep::tool::shapes_same_oriented(s2, &sr);
                if !sso {
                    o = Config::DiffOriented;
                }
            }
            self.set_same_domain_ref(is2, r);
            self.set_same_domain_ori(is2, o);
        }

        // index
        self.set_same_domain_ind_of(s1, 1);
        self.set_same_domain_ind_of(s2, 2);
        Ok(())
    }

    /// `FillShapesSameDomain(S1, S2, c1, c2, refFirst)` (`cxx:955-991`).
    pub fn fill_shapes_same_domain_with_config(
        &mut self,
        s1: &TopoShape,
        s2: &TopoShape,
        c1: Config,
        c2: Config,
        ref_first: bool,
    ) -> Result<(), String> {
        let is1 = self.add_shape_with_ancestor(s1, 1);
        let isdef1 = self
            .shapes
            .change_from_index(is1)
            .map(|sd| sd.orientation_def())
            .ok_or_else(|| "FacesFiller::Insert SD 1".to_string())?;
        let mut todef1 = true;
        if c1 == Config::UnshGeometry && isdef1 {
            todef1 = false;
        }

        let is2 = self.add_shape_with_ancestor(s2, 2);
        let isdef2 = self
            .shapes
            .change_from_index(is2)
            .map(|sd| sd.orientation_def())
            .ok_or_else(|| "FacesFiller::Insert SD 1".to_string())?;
        let mut todef2 = true;
        if c2 == Config::UnshGeometry && isdef2 {
            todef2 = false;
        }

        if todef1 || todef2 {
            self.fill_shapes_same_domain(s1, s2, ref_first)?;
        }

        if todef1 && c1 == Config::UnshGeometry {
            self.set_same_domain_ori_of(s1, Config::UnshGeometry);
        }
        if todef2 && c2 == Config::UnshGeometry {
            self.set_same_domain_ori_of(s2, Config::UnshGeometry);
        }
        Ok(())
    }


    // --- same domain / ancestor records ----------------------------------

    /// `SameDomainRef(I)` (`cxx:599-607`).
    pub fn same_domain_ref(&self, i: i32) -> i32 {
        if i >= 1 && i <= self.shapes.extent() {
            if let Some(sd) = self.shapes.find_from_index(i) {
                return sd.same_domain_ref();
            }
        }
        0
    }

    /// `SameDomainRef(S)` (`cxx:610-620`).
    pub fn same_domain_ref_of(&self, s: &TopoShape) -> i32 {
        if s.is_null() {
            return 0;
        }
        match self.shapes.find_from_key(s) {
            Some(sd) => sd.same_domain_ref(),
            None => 0,
        }
    }

    /// `SameDomainRef(I, Ref)` (`cxx:625-632`).
    pub fn set_same_domain_ref(&mut self, i: i32, r: i32) {
        if i >= 1 && i <= self.shapes.extent() {
            if let Some(sd) = self.shapes.change_from_index(i) {
                sd.set_same_domain_ref(r);
            }
        }
    }

    /// `SameDomainRef(S, Ref)` (`cxx:636-646`).
    pub fn set_same_domain_ref_of(&mut self, s: &TopoShape, r: i32) {
        if s.is_null() {
            return;
        }
        if let Some(sd) = self.shapes.change_from_key(s) {
            sd.set_same_domain_ref(r);
        }
    }

    /// `SameDomainOri(I)` (`cxx:651-658`).
    pub fn same_domain_ori(&self, i: i32) -> Config {
        if i >= 1 && i <= self.shapes.extent() {
            if let Some(sd) = self.shapes.find_from_index(i) {
                return sd.same_domain_ori();
            }
        }
        Config::UnshGeometry
    }

    /// `SameDomainOri(S)` (`cxx:662-671`).
    pub fn same_domain_ori_of(&self, s: &TopoShape) -> Config {
        if !s.is_null() {
            if let Some(sd) = self.shapes.find_from_key(s) {
                return sd.same_domain_ori();
            }
        }
        Config::UnshGeometry
    }

    /// `SameDomainOri(I, Ori)` (`cxx:676-682`).
    pub fn set_same_domain_ori(&mut self, i: i32, ori: Config) {
        if i >= 1 && i <= self.shapes.extent() {
            if let Some(sd) = self.shapes.change_from_index(i) {
                sd.set_same_domain_ori(ori);
            }
        }
    }

    /// `SameDomainOri(S, Ori)` (`cxx:687-696`).
    pub fn set_same_domain_ori_of(&mut self, s: &TopoShape, ori: Config) {
        if !s.is_null() {
            if let Some(sd) = self.shapes.change_from_key(s) {
                sd.set_same_domain_ori(ori);
            }
        }
    }

    /// `SameDomainInd(I)` (`cxx:701-709`).
    pub fn same_domain_ind(&self, i: i32) -> i32 {
        if i >= 1 && i <= self.shapes.extent() {
            if let Some(sd) = self.shapes.find_from_index(i) {
                return sd.same_domain_ind();
            }
        }
        0
    }

    /// `SameDomainInd(S)` (`cxx:712-721`).
    pub fn same_domain_ind_of(&self, s: &TopoShape) -> i32 {
        if !s.is_null() {
            if let Some(sd) = self.shapes.find_from_key(s) {
                return sd.same_domain_ind();
            }
        }
        0
    }

    /// `SameDomainInd(I, Ind)` (`cxx:726-733`).
    pub fn set_same_domain_ind(&mut self, i: i32, ind: i32) {
        if i >= 1 && i <= self.shapes.extent() {
            if let Some(sd) = self.shapes.change_from_index(i) {
                sd.set_same_domain_ind(ind);
            }
        }
    }

    /// `SameDomainInd(S, Ind)` (`cxx:737-748`).
    pub fn set_same_domain_ind_of(&mut self, s: &TopoShape, ind: i32) {
        if s.is_null() {
            return;
        }
        if let Some(sd) = self.shapes.change_from_key(s) {
            sd.set_same_domain_ind(ind);
        }
    }

    /// `AncestorRank(I)` (`cxx:752-760`).
    pub fn ancestor_rank(&self, i: i32) -> i32 {
        if i >= 1 && i <= self.shapes.extent() {
            if let Some(sd) = self.shapes.find_from_index(i) {
                return sd.ancestor_rank();
            }
        }
        0
    }

    /// `AncestorRank(S)` (`cxx:763-773`).
    pub fn ancestor_rank_of(&self, s: &TopoShape) -> i32 {
        if s.is_null() {
            return 0;
        }
        match self.shapes.find_from_key(s) {
            Some(sd) => sd.ancestor_rank(),
            None => 0,
        }
    }

    /// `AncestorRank(I, Ianc)` (`cxx:778-784`).
    pub fn set_ancestor_rank(&mut self, i: i32, ianc: i32) {
        if i >= 1 && i <= self.shapes.extent() {
            if let Some(sd) = self.shapes.change_from_index(i) {
                sd.set_ancestor_rank(ianc);
            }
        }
    }

    /// `AncestorRank(S, Ianc)` (`cxx:789-800`).
    pub fn set_ancestor_rank_of(&mut self, s: &TopoShape, ianc: i32) {
        if s.is_null() {
            return;
        }
        if let Some(sd) = self.shapes.change_from_key(s) {
            sd.set_ancestor_rank(ianc);
        }
    }
}
