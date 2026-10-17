//! `StepToTopoDS_TranslateShell` / `TranslateSolid` with `StepToTopoDS_Builder`: `CLOSED_SHELL`, `OPEN_SHELL`, `SHELL_BASED_SURFACE_MODEL`, `MANIFOLD_SOLID_BREP`, `BREP_WITH_VOIDS` and the oriented closed shell.
use super::*;

impl<'a> Resolver<'a> {
    pub(in crate::step) fn resolve_shell(&self, rec: &'a Record) -> Result<TopoShape, String> {
        // `FromSTEP.exec.op` = `FixShape`: one `ShapeFix_Shape` pass over the
        // root shape with one `ShapeBuild_ReShape` (`ShapeFix_Shape.cxx:75`),
        // closed by `Context()->Apply(S)` (`cxx:257`).
        let owned = self.open_heal_context();
        let out = self.resolve_shell_faces(rec);
        match out {
            Ok(shell) => {
                // `ShapeFix_Shell::Perform` runs `FixFaceOrientation` after the
                // per-face `ShapeFix_Face` pass (`ShapeFix_Shell.cxx:139-142`),
                // on `Context()->Apply(myShell)`.
                let applied = self.close_heal_context(&owned, shell);
                Ok(crate::shhealing::fix_shell_face_orientation(applied))
            }
            Err(e) => {
                self.pop_heal_context(&owned);
                Err(e)
            }
        }
    }

    pub(super) fn resolve_shell_faces(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let faces = parse_ref_list(&rec.args[1]);
        let mut fs = Vec::with_capacity(faces.len());
        for &it in &faces {
            let s = self.resolve_shape(it)?;
            if s.is_face() {
                fs.push(Face(s));
            } else if s.is_shell() || s.is_compound() {
                // `ShapeBuild_ReShape::applyImpl` (`cxx:282-299`): when the
                // replacement has a different type than the original, its own
                // children of the original's type are added instead. A
                // `FixMissingSeam` shell and a `FixSplitFace` compound
                // (`ShapeFix_Face.cxx:2988-2999`) therefore contribute their
                // faces here, which is how OCCT's `Context()->Apply(S)`
                // (`ShapeFix_Shape.cxx:257`) turns 1772 imported faces into
                // 1778 (`FromSTEP.exec.op = FixShape`).
                let faces = crate::topo_tools_full::faces_of(&s);
                if faces.is_empty() {
                    return Err(format!("#{it}: expected FACE in CLOSED_SHELL"));
                }
                fs.extend(faces);
            } else {
                return Err(format!("#{it}: expected FACE in CLOSED_SHELL"));
            }
        }
        Ok(self.b.make_shell(&fs).0)
    }

    pub(in crate::step) fn resolve_solid(&self, rec: &'a Record) -> Result<TopoShape, String> {
        // One `FixShape` pass per transferred root: every shell of the solid
        // shares the context (`ShapeFix_Shape.cxx:150-171` -> `ShapeFix_Solid`
        // -> `ShapeFix_Shell`, all with `Context()`).
        let owned = self.open_heal_context();
        let out = self.resolve_solid_shells(rec);
        match out {
            Ok(solid) => Ok(self.close_heal_context(&owned, solid)),
            Err(e) => {
                self.pop_heal_context(&owned);
                Err(e)
            }
        }
    }

    pub(super) fn resolve_solid_shells(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let outer = parse_ref(&rec.args[1]).ok_or("MANIFOLD_SOLID_BREP: bad shell ref")?;
        let s = self.resolve_shape(outer)?;
        if !s.is_shell() {
            return Err("MANIFOLD_SOLID_BREP: outer is not a shell".into());
        }
        // `ShapeFix_Solid::Perform` single-shell arm (`ShapeFix_Solid.cxx:520-565`)
        // calls `SolidFromShell` (`:655-702`) after `ShapeFix_Shell`.
        let mut s = s;
        crate::shhealing::orient_single_shell_solid(&mut s);
        Ok(self.b.make_solid(&[Shell(s)]).0)
    }

    /// `StepToTopoDS_Builder::Init(BrepWithVoids)` (`StepToTopoDS_Builder.cxx:
    /// 179-251`), reached from `STEPControl_ActorRead::TransferEntity`
    /// (`cxx:1827-1831`, checked before `MANIFOLD_SOLID_BREP` because
    /// `BREP_WITH_VOIDS` is a subtype). The outer `CLOSED_SHELL` becomes the
    /// solid's first shell; every `ORIENTED_CLOSED_SHELL` of the `voids` set
    /// becomes an inner shell, reversed when its `orientation` attribute is
    /// `.F.` (`cxx:237-241`).
    pub(in crate::step) fn resolve_brep_with_voids(&self, rec: &'a Record) -> Result<TopoShape, String> {
        // One `FixShape` pass for the whole solid, outer shell and voids alike
        // (`ShapeFix_Shape.cxx:158-160` -> `ShapeFix_Solid`).
        let owned = self.open_heal_context();
        let out = self.resolve_brep_with_voids_shells(rec);
        match out {
            Ok(solid) => Ok(self.close_heal_context(&owned, solid)),
            Err(e) => {
                self.pop_heal_context(&owned);
                Err(e)
            }
        }
    }

    pub(super) fn resolve_brep_with_voids_shells(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let outer_ref = parse_ref(&rec.args[1]).ok_or("BREP_WITH_VOIDS: bad outer ref")?;
        let outer = self.resolve_shape(outer_ref)?;
        if !outer.is_shell() {
            // `cxx:210-214`: "OuterShell from BrepWithVoids not mapped to TopoDS".
            return Err("BREP_WITH_VOIDS: outer is not a shell".into());
        }
        let mut shells = vec![Shell(outer)];
        for void in parse_ref_list(&rec.args[2]) {
            match self.resolve_oriented_closed_shell(void) {
                Ok(s) => shells.push(Shell(s)),
                // `cxx:244-247`: "A Void from BrepWithVoids not mapped to
                // TopoDS" is a warning; the solid keeps the shells that mapped.
                Err(e) => self.warn(format!("void #{void}: {e}")),
            }
        }
        Ok(self.b.make_solid(&shells).0)
    }

    /// `StepToTopoDS_Builder::Init(ShellBasedSurfaceModel)`
    /// (`StepToTopoDS_Builder.cxx:384-459`), reached from
    /// `STEPControl_ActorRead::TransferEntity` (`cxx:1837-1843`). Every entry of
    /// `sbsm_boundary` (a `shell` select) goes through
    /// `StepToTopoDS_TranslateShell::Init` (`StepToTopoDS_TranslateShell.cxx:
    /// 53-118`) and is added to a compound; an open shell is flagged
    /// `Closed(false)` (`StepToTopoDS_Builder.cxx:425`) and a closed one
    /// `Closed(true)` (`cxx:442`). The result is that compound when the boundary
    /// set holds more than one entry, else the last shell (`cxx:452-458`). A
    /// shell that does not map is a warning, not a failure (`cxx:431-433`,
    /// `:447-449`).
    pub(in crate::step) fn resolve_shell_based_surface_model(
        &self,
        rec: &'a Record,
    ) -> Result<TopoShape, String> {
        let boundaries = parse_ref_list(&rec.args[1]);
        let mut shells: Vec<TopoShape> = Vec::with_capacity(boundaries.len());
        let mut last: Option<TopoShape> = None;
        for &b in &boundaries {
            // `aShell.OpenShell()` / `aShell.ClosedShell()`: a `shell` select
            // that is neither (`VERTEX_SHELL` / `WIRE_SHELL`) takes no branch.
            let kind = match self.records.get(&b) {
                Some(r) => r.type_name.clone(),
                None => continue,
            };
            let closed = match kind.as_str() {
                "OPEN_SHELL" => false,
                "CLOSED_SHELL" => true,
                _ => continue,
            };
            match self.resolve_shape(b) {
                Ok(shell) if shell.is_shell() => {
                    shell.set_closed(closed);
                    last = Some(shell.clone());
                    shells.push(shell);
                }
                Ok(_) => self.warn(format!(
                    "#{b}: shell from SHELL_BASED_SURFACE_MODEL is not a shell"
                )),
                Err(e) => self.warn(format!(
                    "#{b}: shell from SHELL_BASED_SURFACE_MODEL not mapped to TopoDS: {e}"
                )),
            }
        }
        if boundaries.len() > 1 {
            // `cxx:452-456`: the compound, even when no boundary mapped.
            return Ok(self.b.make_compound_of(&shells).0);
        }
        last.ok_or_else(|| "SHELL_BASED_SURFACE_MODEL: no shell mapped".into())
    }

    /// One void of a `BREP_WITH_VOIDS`: `ORIENTED_CLOSED_SHELL(name, *,
    /// closed_shell, orientation)`. Its faces are the derived
    /// `SELF\connected_face_set.cfs_faces := closed_shell.cfs_faces`
    /// (`StepShape_OrientedClosedShell.cxx:81-87`), so the shell is translated
    /// through the referenced `CLOSED_SHELL`.
    pub(super) fn resolve_oriented_closed_shell(&self, id: usize) -> Result<TopoShape, String> {
        let rec = self.record(id)?;
        if rec.type_name != "ORIENTED_CLOSED_SHELL" {
            return Err(format!("expected ORIENTED_CLOSED_SHELL at #{id}"));
        }
        let shell_ref = parse_ref(&rec.args[2]).ok_or("ORIENTED_CLOSED_SHELL: bad shell ref")?;
        let mut s = self.resolve_shape(shell_ref)?;
        if !s.is_shell() {
            return Err(format!("#{shell_ref}: void is not a shell"));
        }
        // `StepToTopoDS_Builder.cxx:239-241`:
        // `if (!anOCShell->Orientation()) aShape.Reverse();`
        if !parse_logical(rec.args.get(3).map(String::as_str), true) {
            s.set_orientation(Orientation::Reversed);
        }
        Ok(s)
    }

}
