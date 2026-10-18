//! The incremental writer and its emit plumbing, plus shell/solid assembly and the sectioned file layout. Source: `IGESControl_Writer` -> `Interface_InterfaceModel::AddWithRefs`, `BRepToIGES_BRShell` / `BRepToIGES_BRSolid`.
use super::prelude::*;
use super::*;

/// Incremental IGES writer holding the entity list.
pub(super) struct IgesWriter {
    pub(super) entities: Vec<Ent>,
    /// Top-level entity of every shape handed to [`IgesWriter::emit_shape`], in
    /// that order. `IGESControl_Writer::AddEntity` feeds each one to
    /// `Interface_InterfaceModel::AddWithRefs` (`IGESControl_Writer.cxx:243-252`),
    /// so these are the DFS roots of the written model (T-85 step 2).
    pub(super) roots: Vec<usize>,
    pub(super) point_entities: HashMap<usize, usize>,
    /// Largest absolute model coordinate (`IGESControl_Writer.cxx:180-189`,
    /// `MaxMaxCoords` of the shape box). `0` leaves global field 20 empty,
    /// which is `hasMaxCoord == false` (`IGESData_GlobalSection.cxx:519-527`).
    pub(super) max_coord: f64,
}

impl IgesWriter {
    pub(super) fn new() -> Self {
        Self {
            entities: Vec::new(),
            roots: Vec::new(),
            point_entities: HashMap::new(),
            max_coord: 0.0,
        }
    }

    /// Register one entity. The DE label, line weight, color and status fields are
    /// not parameters here: OCCT's writer takes them from the entity's own
    /// `IGESData_IGESEntity` state, which is at its defaults for generated
    /// geometry (`IGESData_IGESWriter.cxx:276-365`), see
    /// [`Ent::directory_lines`].
    pub(super) fn emit(&mut self, ty: i32, form: i32, params: String) -> usize {
        self.entities.push(Ent {
            ty,
            form,
            trsf: None,
            refs: Vec::new(),
            params,
            status: [0; 4],
        });
        self.entities.len()
    }

    /// Register one entity that references others, recording those pointers in
    /// record order (T-85's prerequisite for OCCT's reachability-based writing).
    ///
    /// The parameter text carries a `#k` placeholder for the `k`-th entry of
    /// `refs` instead of a DE number; [`IgesWriter::final_entities`] substitutes
    /// the final number after the reachability renumbering. OCCT can write the
    /// numbers directly because its model is numbered while it is built
    /// (`Interface_InterfaceModel::AddWithRefs`, `Interface_InterfaceModel.cxx:652-692`);
    /// the port numbers at write time, so the reference *index* is what the
    /// parameter text can hold.
    pub(super) fn emit_refs(&mut self, ty: i32, form: i32, params: String, refs: &[usize]) -> usize {
        debug_assert!(
            pointers_all_known(&params, refs.len()),
            "entity type {ty}: parameter text has placeholders outside 0..{}",
            refs.len()
        );
        let de = self.emit(ty, form, params);
        self.entities[de - 1].refs = refs.to_vec();
        de
    }

    /// Every recorded pointer must address an entity that exists and precede the
    /// referencing entity (`IGESData_IGESWriter::Send` resolves pointers through
    /// the model, so a dangling number would be written verbatim). Debug-only:
    /// the emitters below build the pointers themselves.
    pub(super) fn check_refs(&self) {
        for (i, e) in self.entities.iter().enumerate() {
            for r in &e.refs {
                debug_assert!(
                    *r >= 1 && *r <= self.entities.len() && *r != i + 1,
                    "entity {} (type {}) references {} - out of range or self",
                    i + 1,
                    e.ty,
                    r
                );
            }
        }
    }

    /// The entities OCCT's `AddWithRefs` would put in the model, in the order it
    /// adds them, renumbered and with the `#k` placeholders resolved.
    ///
    /// `Interface_InterfaceModel::AddWithRefs` (`Interface_InterfaceModel.cxx:652-692`)
    /// adds an entity and then recurses into its shared entities **in order**,
    /// skipping anything already added: a pre-order DFS from the roots
    /// (`IGESControl_Writer::AddEntity`, `IGESControl_Writer.cxx:243-252`). The
    /// shared entities are `IGESData_GeneralModule::FillSharedCase`'s list —
    /// the DE-part entities first (only field 7, the transformation matrix, is
    /// ever set by this writer) and then the own-parameter references, which is
    /// the record order of [`Ent::refs`].
    pub(super) fn final_entities(&self) -> Vec<Ent> {
        let mut order: Vec<usize> = Vec::new();
        let mut seen: HashSet<usize> = HashSet::new();
        let mut stack: Vec<usize> = self.roots.iter().rev().copied().collect();
        while let Some(old) = stack.pop() {
            if old < 1 || old > self.entities.len() || !seen.insert(old) {
                continue;
            }
            order.push(old);
            let e = &self.entities[old - 1];
            let mut kids: Vec<usize> = Vec::new();
            if let Some(t) = e.trsf {
                kids.push(t);
            }
            kids.extend(e.refs.iter().copied());
            for k in kids.into_iter().rev() {
                stack.push(k);
            }
        }

        let new_of: HashMap<usize, usize> =
            order.iter().enumerate().map(|(i, old)| (*old, i + 1)).collect();
        order
            .iter()
            .map(|old| {
                let e = &self.entities[old - 1];
                Ent {
                    ty: e.ty,
                    form: e.form,
                    // `v[6] = themodel->DNum(trsf)` (`IGESData_IGESWriter.cxx:324-331`):
                    // the transformation matrix field is an entity pointer like
                    // any other, so it carries the first Directory card number.
                    // A pointer that is not in `new_of` cannot happen (every
                    // child was pushed onto the DFS stack); keep the old value
                    // rather than emitting a hole.
                    trsf: e.trsf.map(|t| {
                        new_of
                            .get(&t)
                            .map(|n| 2 * n - 1)
                            .unwrap_or(t)
                    }),
                    refs: e
                        .refs
                        .iter()
                        .map(|r| new_of.get(r).copied().unwrap_or(*r))
                        .collect(),
                    params: resolve_placeholders(&e.params, &e.refs, &new_of),
                    status: [0; 4],
                }
            })
            .collect()
    }

    /// `IGESData_IGESEntity::InitTransf` — record the entity's transformation
    /// matrix pointer, written to DE card 1 field 7 (`IGESData_IGESWriter.cxx:324-331`).
    pub(super) fn set_trsf(&mut self, de: usize, trsf: usize) {
        self.entities[de - 1].trsf = Some(trsf);
    }


    /// `Group` (402) when there is more than one item, the item itself otherwise -
    /// `IGESBasic_ToolGroup::WriteOwnParams` (`IGESBasic_ToolGroup.cxx:104-115`)
    /// writes `402, n, entity...;` and both `BRepToIGES_BRShell::TransferShell`
    /// (`BRepToIGES_BRShell.cxx:463-471`) and
    /// `BRepToIGES_BRSolid::TransferSolid` (`BRepToIGES_BRSolid.cxx:154-163`) use
    /// exactly this rule. The form is 1: `IGESBasic_Group::IGESBasic_Group`
    /// calls `InitTypeAndForm(402, 1)` (`IGESBasic_Group.cxx:29-31`), and
    /// `DirChecker` requires that form (`IGESBasic_ToolGroup.cxx:161`).
    pub(super) fn group_or_single(&mut self, items: Vec<usize>) -> Option<usize> {
        match items.len() {
            0 => None,
            1 => Some(items[0]),
            n => {
                Some(self.emit_refs(402, 1, format!("402,{n},{};", ph_list(n)), &items))
            }
        }
    }

    /// `BRepToIGES_BRShell::TransferShell` (`BRepToIGES_BRShell.cxx:411-476`): the
    /// shell's faces are transferred as trimmed surfaces (144) and grouped.
    pub(super) fn emit_shell(&mut self, sh: &Shell) -> Option<usize> {
        if std::env::var("IGES_TRACE_FACE").is_ok() {
            eprintln!("TRACESHELL ori={:?}", sh.0.orientation());
        }
        // `TopExp_Explorer(shell, FACE)` composes the shell orientation
        // (`BRepToIGES_BRShell.cxx:433`). A shell `SolidFromShell` reversed
        // (`ShapeFix_Solid.cxx:686`) therefore flips every face.
        let faces: Vec<TopoShape> = crate::iterator::cumulated_children(&sh.0)
            .into_iter()
            .filter(|h| h.shape_type() == ShapeType::Face)
            .collect();
        let refs: Vec<usize> = faces
            .iter()
            .filter_map(|f| self.emit_face(&Face(f.clone())))
            .collect();
        self.group_or_single(refs)
    }

    /// `BRepToIGES_BRSolid::TransferSolid` (`BRepToIGES_BRSolid.cxx:100-168`): the
    /// solid's shells are transferred and grouped the same way.
    pub(super) fn emit_solid(&mut self, s: &Solid) -> Option<usize> {
        // `TopExp_Explorer(solid, SHELL)` (`BRepToIGES_BRSolid.cxx:124`).
        let shells: Vec<TopoShape> = crate::iterator::cumulated_children(&s.0)
            .into_iter()
            .filter(|h| h.shape_type() == ShapeType::Shell)
            .collect();
        let refs: Vec<usize> = shells
            .iter()
            .filter_map(|sh| self.emit_shell(&Shell(sh.clone())))
            .collect();
        self.group_or_single(refs)
    }

    pub(super) fn note_max_coord(&mut self, shape: &TopoShape) {
        // `BRepBndLib::Add` then skip a void or open box
        // (`IGESControl_Writer.cxx:181-189`).
        let box_ = crate::bopalgo_tools_wires::wires_shape_box(shape);
        if box_.is_void()
            || box_.is_open_xmin()
            || box_.is_open_xmax()
            || box_.is_open_ymin()
            || box_.is_open_ymax()
            || box_.is_open_zmin()
            || box_.is_open_zmax()
        {
            return;
        }
        let Some((xmin, xmax, ymin, ymax, zmin, zmax)) = box_.get() else {
            return;
        };
        for v in [xmin, xmax, ymin, ymax, zmin, zmax] {
            let a = v.abs();
            if a.is_finite() && a > self.max_coord {
                self.max_coord = a;
            }
        }
    }

    pub(super) fn emit_shape(&mut self, shape: &TopoShape) {
        self.note_max_coord(shape);
        match shape.shape_type() {
            ShapeType::Compound => {
                let kids: Vec<TopoShape> = shape
                    .tshape
                    .read()
                    .unwrap()
                    .children
                    .clone();
                for k in kids {
                    self.emit_shape(&k);
                }
            }
            ShapeType::Solid => {
                if let Some(r) = self.emit_solid(&Solid(shape.clone())) {
                    self.roots.push(r);
                }
            }
            ShapeType::Shell => {
                if let Some(r) = self.emit_shell(&Shell(shape.clone())) {
                    self.roots.push(r);
                }
            }
            ShapeType::Face => {
                // `BRepToIGES_BRShell.cxx:257-261`: a face with a null basis
                // surface entity contributes nothing to the model.
                if let Some(r) = self.emit_face(&Face(shape.clone())) {
                    self.roots.push(r);
                }
            }
            _ => {}
        }
    }

    // ---- file assembly ----

    pub(super) fn global_lines(&self) -> Vec<String> {
        // `IGESData::Init` (`IGESData.cxx:249-279`) plus the per-shape update
        // in `IGESControl_Writer::AddShape` (`IGESControl_Writer.cxx:177-189`).
        // Default comma and semicolon are written as two empty fields
        // (`IGESData_GlobalSection::Params`, `cxx:445-462`), not as `1H,`.
        // Single and double precision are both `RealLast10Exp()` / `RealDigits()`
        // (`cxx:258-261`), which is `308` / `15`, and the version is
        // `SetIGESVersion(11)` (`cxx:276`). The date is `YYYYMMDD.HHMMSS`
        // (`NewDateString`, `cxx:247`); its Hollerith count has to be the
        // character length or every later field, including the version, shifts.
        let date = iges_date_now();
        let product = "Open CASCADE IGES processor 8.0";
        let system = "Open CASCADE 8.0";
        let max_coord = if self.max_coord > 0.0 {
            num(self.max_coord)
        } else {
            String::new()
        };
        let data = format!(
            ",,{product_h},{file_h},{system_h},{product_h},32,308,15,308,15,,{scale},2,{unit},1,{weight},{date_h},{resolution},{max_coord},,,11,0,{date_h},;",
            product_h = hollerith(product),
            file_h = hollerith("Filename.iges"),
            system_h = hollerith(system),
            scale = num(1.0),
            unit = hollerith("MM"),
            weight = num(0.01),
            date_h = hollerith(&date),
            resolution = num(occt_core::CONFUSION),
        );
        let mut out = Vec::new();
        let mut s = data.as_str();
        while !s.is_empty() {
            let end = s.len().min(MAXCARS_G);
            out.push(s[..end].to_string());
            s = &s[end..];
        }
        out
    }

    pub(super) fn finish(self) -> String {
        self.check_refs();
        // T-85 step 2: OCCT's model holds only the entities reachable from the
        // shapes handed to `AddEntity` (`AddWithRefs`), numbered in the order it
        // adds them; the placeholders in the parameter text become the final
        // numbers here.
        let mut entities = self.final_entities();
        // `IGESControl_Writer::ComputeModel` (`IGESControl_Writer.cxx:256-262`).
        compute_status(&mut entities);
        let mut out = String::new();
        let mut s_seq = 1usize;
        let mut g_seq = 1usize;
        let mut d_seq = 1usize;
        let mut p_seq = 1usize;

        out.push_str(&sec_line(
            'S',
            s_seq,
            "IGES B-REP MODEL GENERATED BY RUST OCCT PORT",
            MAXCARS_G,
        ));
        out.push('\n');
        s_seq += 1;

        for gl in self.global_lines() {
            out.push_str(&sec_line('G', g_seq, &gl, MAXCARS_G));
            out.push('\n');
            g_seq += 1;
        }

        // Directory entries: entity `i` owns the two DE pointers `2i-1` and `2i`.
        let mut p_start = 1usize;
        for (i, ent) in entities.iter().enumerate() {
            let (l1, l2) = ent.directory_lines(p_start);
            p_start += ent.param_line_count();
            out.push_str(&sec_line('D', 2 * i + 1, &l1, MAXCARS_G));
            out.push('\n');
            out.push_str(&sec_line('D', 2 * i + 2, &l2, MAXCARS_G));
            out.push('\n');
            d_seq += 2;
        }

        for (i, ent) in entities.iter().enumerate() {
            for chunk in ent.param_chunks() {
                out.push_str(&param_line(&chunk, 2 * i + 1, p_seq));
                out.push('\n');
                p_seq += 1;
            }
        }

        // Terminate card (`IGESData_IGESWriter.cxx:942-947`): the **last sequence
        // number** of each section - `nbs`, `nbg`, `nbd * 2` and
        // `thepnum.Value(thepnum.Length()) - 1` - blank-filled, then 40 blanks and
        // `T0000001`. The Directory section therefore reports `2 * entities`, not
        // `2 * entities + 1` (as of batch 70).
        let t = format!(
            "S{:>7}G{:>7}D{:>7}P{:>7}{}T0000001",
            s_seq - 1,
            g_seq - 1,
            d_seq - 1,
            p_seq - 1,
            " ".repeat(40)
        );
        out.push_str(&t);
        out.push('\n');
        out
    }
}
