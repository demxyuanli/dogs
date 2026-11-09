//! Entity bookkeeping: the directory-entry metadata and status of one written entity. Source: `IGESData_IGESEntity` / `IGESData_DirChecker`, filled by `IGESControl_Writer::ComputeModel`.
use super::*;

// ---------------------------------------------------------------------------
// Entity bookkeeping
// ---------------------------------------------------------------------------

/// One IGES entity: directory-entry metadata + parameter-data string.
pub(super) struct Ent {
    pub(super) ty: i32,
    pub(super) form: i32,
    /// Pointer to the entity's transformation matrix (`#124`), `HasTransf()` in
    /// `IGESData_IGESEntity`; written to DE card 1 field 7
    /// (`IGESData_IGESWriter.cxx:324-331`, `v[6] = themodel->DNum(...)`).
    pub(super) trsf: Option<usize>,
    /// The other entities this one references, in record order - `OwnShared` in
    /// OCCT terms (`IGESData_GeneralModule::OwnShared` feeds the writer's pointer
    /// resolution). Recorded so the pointers can be checked (and, for T-85, later
    /// remapped) without parsing the parameter text back into numbers.
    pub(super) refs: Vec<usize>,
    /// Parameter data, always ending with the record delimiter `;`.
    pub(super) params: String,
    /// The entity's status number, DE card 1 columns 65-72: BlankStatus,
    /// SubordinateStatus, UseFlag and HierarchyStatus, two columns each.
    /// Filled by [`compute_status`] from the finished model;
    /// `IGESData_IGESEntity` starts at `InitStatus(0, 0, 0, 0)`
    /// (`IGESData_IGESEntity.cxx:53-63`).
    pub(super) status: [i32; 4],
}

impl Ent {
    /// Number of Parameter cards this entity occupies, which card 2 of its
    /// directory entry reports (`IGESData_IGESWriter.cxx:835`:
    /// `v[15] = thepnum.Value(i + 1) - thepnum.Value(i)`).
    pub(super) fn param_line_count(&self) -> usize {
        self.param_chunks().len()
    }

    /// The entity's Parameter cards, as `IGESData_IGESWriter` builds them.
    ///
    /// `OwnParams` clears the card buffer and writes the type number, without a
    /// separator (`IGESData_IGESWriter.cxx:418-434`); every own parameter is
    /// then written by a `Send`, i.e. the separator followed by the value
    /// (`cxx:560-565` for integers, `:586-589` for reals, `:610-620` for
    /// pointers — a *void* parameter writes the separator alone, `:557-560`);
    /// `EndEntity` closes the record with `theendm` and dumps the last card
    /// (`:471-483`). Splitting the parameter text on the separator and feeding
    /// the pieces to [`CardBuffer`] is that same sequence, so the cards break
    /// where OCCT's break instead of every 64 columns.
    pub(super) fn param_chunks(&self) -> Vec<String> {
        let body = self.params.strip_suffix(';').unwrap_or(self.params.as_str());
        let mut buf = CardBuffer::new(MAXCARS_P);
        for (i, field) in body.split(',').enumerate() {
            if i > 0 {
                buf.add_char(b',');
            }
            buf.add_string(field);
        }
        buf.add_char(b';');
        if !buf.line.is_empty() {
            buf.flush();
        }
        buf.cards
    }

    /// The two Directory Entry cards of the entity
    /// (`IGESData_IGESWriter.cxx:276-365` filling `v[0..16]`, `:809-880` laying
    /// the cards out): eight 8-column fields plus four 2-column status fields on
    /// the first card and nine 8-column fields on the second - 72 data columns
    /// each. The four status fields are printed `%2.2d`
    /// (`IGESData_IGESWriter.cxx:836-840`), i.e. zero filled, not blank filled:
    /// a status of 0 is `00`, so card 1 of a default entity ends on `00000000`.
    /// Card 1 field 2 is the sequence number of the entity's **first
    /// parameter card** (`v[1] = thepnum.Value(i)`, `cxx:834`), field 7 the
    /// transformation matrix pointer (`cxx:324-331`); card 2 field 4 is the number
    /// of parameter cards (`v[15]`, `cxx:835`), field 5 the form number
    /// (`v[16] = anent->FormNumber()`, `cxx:364`).
    ///
    /// The status number (`v[8..11]`) is [`Ent::status`], computed on the
    /// finished model by [`compute_status`].
    ///
    /// The remaining fields carry `IGESData_IGESEntity`'s defaults: line weight
    /// `theLWeightNum = 0` and color `DefColor() == DefVoid` ⇒ 0
    /// (`IGESData_IGESEntity.cxx:53-63`, `IGESData_IGESWriter.cxx:347-361`), a
    /// blank label (`theShortLabel` is null unless `SetShortLabel` is called,
    /// `cxx:366-379`) and the subscript number `theSubScriptN = 0`
    /// (`IGESData_IGESEntity.cxx:60`, written right-justified by `cxx:380-391`).
    pub(super) fn directory_lines(&self, p_start: usize) -> (String, String) {
        let ty = field8(&self.ty.to_string());
        let pstart = field8(&p_start.to_string());
        let pcount = field8(&self.param_line_count().to_string());
        let form = field8(&self.form.to_string());
        let trsf = field8(&self.trsf.map_or("0".to_string(), |t| t.to_string()));
        let l1 = format!(
            "{ty}{pstart}{}{}{}{}{trsf}{}{:02}{:02}{:02}{:02}",
            field8("0"),
            field8("0"),
            field8("0"),
            field8("0"),
            field8("0"),
            self.status[0],
            self.status[1],
            self.status[2],
            self.status[3]
        );
        let l2 = format!(
            "{ty}{}{}{pcount}{form}{}{}{}{}",
            field8("0"),
            field8("0"),
            field8(""),
            field8(""),
            field8(""),
            field8("0")
        );
        (l1, l2)
    }
}

/// `IGESData_BasicEditor::ComputeStatus` (`IGESData_BasicEditor.cxx:216-328`)
/// followed by `AutoCorrectModel` -> `AutoCorrect` (`:330-433`), which
/// `IGESControl_Writer::ComputeModel` runs before printing
/// (`IGESControl_Writer.cxx:256-262`). Both work on the numbered model, i.e.
/// after [`IgesWriter::final_entities`], and fill DE card 1's status number.
///
/// ComputeStatus has two phases, documented as such at `cxx:231-255`:
///
/// * **Subordinate status** (`cxx:257-281`): every entity ORs `1` into the
///   status of each entity its *own* parameters reference
///   (`gmod->OwnSharedCase` = [`Ent::refs`]), and `2` instead of `1` when the
///   referring entity is an associativity (type 402 or 404). DE-part
///   references are not `OwnShared`, so the transformation matrix
///   ([`Ent::trsf`], card 1 field 7) contributes nothing - the oracle confirms
///   it: all 347 type-124 entries are `00000000`, while every entity
///   referenced from a parameter record is `00010000` or `00020000`.
/// * **UseFlag** (`cxx:283-309`): not ported. It propagates only for types
///   200-299 (`G.GetFromEntity(ent, true, 1)`), 134, 116 and 132
///   (`G.GetFromEntity(ent, true, 4)`); the port's transfer classes create
///   none of those, and the phase reads an `Interface_Graph` status array
///   whose shared sets include the DE-part pointers this writer does not
///   model. `G.Status(i)` is therefore always 0, so the second phase
///   (`cxx:315-327`) keeps each entity's own UseFlag when it is non-zero and
///   takes 0 otherwise.
///
/// `AutoCorrect` then applies the general module's `DirChecker(CN, ent)`:
/// `Correct` forces a status field whenever it differs from the
/// `...Required(val)` the tool declared (`IGESData_DirChecker.cxx:422-451`),
/// while a field left at its `-100` "do not test" default
/// (`IGESData_DirChecker.cxx:33-36`) is untouched (`thegraphier == -100` also
/// disables the DE graphics cleanup, `:400-420`). Of the types this writer
/// emits (100/102/104/108/110/120/124/126/128/140/142/144/402) exactly two declare
/// a status:
///
/// * 142 - `UseFlagRequired(5)` (`IGESGeom_ToolCurveOnSurface.cxx:200-211`),
///   the source of the oracle's `00010500` on every type 142;
/// * 144 - `UseFlagRequired(0)` (`IGESGeom_ToolTrimmedSurface.cxx:264-276`),
///   which is a no-op for an entity whose UseFlag is already 0.
///
/// `AutoCorrect`'s specific-module call `smod->OwnCorrect(CN, ent)`
/// (`cxx:413-418`) is deliberately *not* ported: in OCCT 8.0.0 the 142
/// correction is registered under case 9 but casts to `IGESGeom_Boundary`
/// (`IGESGeom_SpecificModule.cxx:337-345`), so `anent.IsNull()` and the
/// CurveUV/Curve3D of a 142 keep UseFlag 0 - which is what the oracle shows.
pub(super) fn compute_status(entities: &mut [Ent]) {
    let nb = entities.len();
    // NCollection_Array1<int> subs(0, nb); subs.Init(0) (`cxx:226-228`).
    let mut subs = vec![0i32; nb + 1];
    for e in entities.iter() {
        let bit = if e.ty == 402 || e.ty == 404 { 2 } else { 1 };
        for r in &e.refs {
            // themodel->Number(sh.Value()) is 0 for anything outside the model;
            // index 0 of `subs` absorbs those (`cxx:269-280`).
            if let Some(s) = subs.get_mut(*r) {
                *s |= bit;
            }
        }
    }

    for (i, e) in entities.iter_mut().enumerate() {
        // InitStatus(bl, subs.Value(i), uf, hy), `uf == 0` taking G.Status(i)
        // (`cxx:315-327`). Blank and hierarchy status stay as the entity set
        // them - the port never changes them from 0.
        e.status = [e.status[0], subs[i + 1], e.status[2], e.status[3]];

        // DirChecker::Correct (`IGESData_DirChecker.cxx:422-451`).
        let required = match e.ty {
            142 => Some(5), // IGESGeom_ToolCurveOnSurface.cxx:209
            144 => Some(0), // IGESGeom_ToolTrimmedSurface.cxx:272
            _ => None,
        };
        if let Some(val) = required {
            if val != e.status[2] {
                e.status[2] = val;
            }
        }
    }
}

