//! `StepToTopoDS_TranslateMappedItem` with `StepToGeom::MakeTransformation3d` (`StepToGeom.cxx:2153-2214`): assembly instances and the axis-to-axis transform.
use super::*;

impl<'a> Resolver<'a> {
    /// `STEPControl_ActorRead::TransferEntity(StepRepr_MappedItem)`
    /// (`STEPControl_ActorRead.cxx:1971-2045`): an assembly instance. The shape
    /// of `mapping_source.mapped_representation` is placed by the transform
    /// from `mapping_source.mapping_origin` to `mapping_target`, or, when the
    /// mapping target is a `CARTESIAN_TRANSFORMATION_OPERATOR_3D`, by the
    /// transform carried by that operator (`cxx:2013-2032`). `MAPPED_ITEM(name,
    /// mapping_source, mapping_target)` (`RWStepRepr_RWMappedItem.cxx:37-63`);
    /// `REPRESENTATION_MAP(mapping_origin, mapped_representation)`
    /// (`RWStepRepr_RWRepresentationMap.cxx:34-57`).
    pub(in crate::step) fn resolve_mapped_item(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let source_ref = parse_ref(&rec.args[1]).ok_or("MAPPED_ITEM: bad mapping_source ref")?;
        let target_ref = parse_ref(&rec.args[2]).ok_or("MAPPED_ITEM: bad mapping_target ref")?;
        let source = self.record(source_ref)?;
        if source.type_name != "REPRESENTATION_MAP" {
            return Err(format!(
                "MAPPED_ITEM: mapping_source #{source_ref} is {}",
                source.type_name
            ));
        }
        let origin_ref =
            parse_ref(&source.args[0]).ok_or("REPRESENTATION_MAP: bad mapping_origin ref")?;
        let maprep_ref = parse_ref(&source.args[1])
            .ok_or("REPRESENTATION_MAP: bad mapped_representation ref")?;
        // `cxx:1989-2001`: transfer the mapped representation, then warn when it
        // produced no shape.
        let (_, shapes) = self.resolve_representation(maprep_ref)?;
        if shapes.is_empty() {
            return Err(format!(
                "MAPPED_ITEM: mapped representation #{maprep_ref} produced no shape"
            ));
        }
        let mut mapped = if shapes.len() == 1 {
            shapes.into_iter().next().unwrap()
        } else {
            TopoBuilder::new().make_compound_of(&shapes).0
        };
        // `cxx:2013-2032`: two placement formulas.
        let target = self.record(target_ref)?;
        let trsf = if target.type_name == "CARTESIAN_TRANSFORMATION_OPERATOR_3D" {
            self.make_transformation3d(target_ref).ok()
        } else {
            match (self.resolve_axis2(origin_ref), self.resolve_axis2(target_ref)) {
                (Ok(ax_orig), Ok(ax_targ)) => Some(self.compute_axis_transform(&ax_orig, &ax_targ)),
                _ => None,
            }
        };
        match trsf {
            // `cxx:2038 ApplyTransformation`: `shape.Move(TopLoc_Location(Trsf))`.
            // This port's geometry registry is location-blind (the mesh path
            // reads the registered geometry, not `TopoShape::location`), so the
            // placement is baked into a copy, the same mechanism `Assembly`
            // export uses (`shape_ops::transformed_copy`).
            Some(t) if t.form() != occt_core::gp::TrsfForm::Identity => {
                match crate::shape_ops::transformed_copy(&mapped, &t) {
                    Ok(copy) => mapped = copy,
                    Err(e) => self.warn(format!("MAPPED_ITEM: transform failed: {e}")),
                }
            }
            Some(_) => {}
            // `cxx:2042-2044`: "Mapped Item, case not recognized, location
            // ignored" - keep the untransformed shape.
            None => self.warn(format!(
                "MAPPED_ITEM #{target_ref}: case not recognized, location ignored"
            )),
        }
        Ok(mapped)
    }

    /// `StepToGeom::MakeTransformation3d` (`StepToGeom.cxx:2153-2214`) for a
    /// `CARTESIAN_TRANSFORMATION_OPERATOR_3D`. The reader starts at parameter 3
    /// (two `functionally_defined_transformation` names are skipped,
    /// `RWStepGeom_RWCartesianTransformationOperator.cxx:36-49`), so axis1 is
    /// param 4, axis2 param 5, local_origin param 6, scale param 7 and axis3
    /// param 8 (`RWStepGeom_RWCartesianTransformationOperator3d.cxx:82-128`).
    pub(in crate::step) fn make_transformation3d(&self, id: usize) -> Result<occt_core::gp::GpTrsf, String> {
        let rec = self.record(id)?;
        let origin_ref = parse_ref(&rec.args[5]).ok_or("CTO3D: bad local_origin ref")?;
        // `cxx:2158`: `MakeCartesianPoint(LocalOrigin)` - scaled by LengthFactor.
        let p = self.resolve_point(origin_ref)?;
        // `cxx:2163-2198`: default X / Y, overridden by axis1 / axis2.
        let mut d1 = dir_x();
        if let Some(r) = parse_ref(&rec.args[3]) {
            if let Ok(d) = self.resolve_direction(r) {
                d1 = d;
            }
        }
        let mut d2 = dir_y();
        if let Some(r) = parse_ref(&rec.args[4]) {
            if let Ok(d) = self.resolve_direction(r) {
                d2 = d;
            }
        }
        // `cxx:2186-2202`: axis3, or `D1.Crossed(D2)` when absent.
        let d3 = match parse_ref(&rec.args[7]) {
            Some(r) => self
                .resolve_direction(r)
                .unwrap_or_else(|_| d1.crossed(&d2).unwrap_or_else(|_| dir_z())),
            None => d1.crossed(&d2).unwrap_or_else(|_| dir_z()),
        };
        // `gp_Ax3(P, D3, D1)`: Vy = D3 ^ D1, X re-orthogonalized as Vy ^ D3.
        let ydir = d3.crossed(&d1).map_err(|e| format!("CTO3D: {e}"))?;
        let xdir = ydir.crossed(&d3).map_err(|e| format!("CTO3D: {e}"))?;
        let ax3 = GpAx3 {
            axis: GpAx1::new(p, d3),
            vxdir: xdir,
            vydir: ydir,
        };
        let mut t = occt_core::gp::GpTrsf::identity();
        t.set_transformation(&ax3);
        // `cxx:2205-2208`: `if (HasScale) CT.SetScaleFactor(Scale())`.
        if let Some(s) = rec.args.get(6).and_then(|a| parse_f64(a).ok()) {
            t.scale = s;
        }
        // `cxx:2210`: `CT = CT.Inverted()`.
        t.invert().map_err(|e| format!("CTO3D: {e}"))?;
        Ok(t)
    }

    /// `STEPControl_ActorRead::ComputeTransformation` axis-pair arm
    /// (`STEPControl_ActorRead.cxx:2471-2489`):
    /// `Trsf.SetTransformation(ax3Targ, ax3Orig)`
    /// (`gp_Trsf.cxx:172-192`) maps `ax3Orig`-frame coordinates into the
    /// `ax3Targ` frame, i.e. `world->ax3Orig` composed with the inverse of
    /// `world->ax3Targ`.
    pub(in crate::step) fn compute_axis_transform(&self, orig: &GpAx2, targ: &GpAx2) -> occt_core::gp::GpTrsf {
        let ax3_orig = GpAx3 {
            axis: orig.axis,
            vxdir: orig.vxdir,
            vydir: orig.vydir,
        };
        let ax3_targ = GpAx3 {
            axis: targ.axis,
            vxdir: targ.vxdir,
            vydir: targ.vydir,
        };
        let mut t_orig = occt_core::gp::GpTrsf::identity();
        t_orig.set_transformation(&ax3_orig);
        let mut t_targ = occt_core::gp::GpTrsf::identity();
        t_targ.set_transformation(&ax3_targ);
        let inv_targ = t_targ
            .inverted()
            .unwrap_or_else(|_| occt_core::gp::GpTrsf::identity());
        t_orig.multiplied(&inv_targ)
    }

}
