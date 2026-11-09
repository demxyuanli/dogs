//! `HLRBRep_Hider`
//! (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_Hider.{hxx,cxx}`,
//! `hxx:31-53`, `cxx:30-852`): removes from the edges of one `HLRBRep_Data`
//! the parts hidden by hiding face number `FI`.
//!
//! The pass, in OCCT order:
//!
//! 1. `InitEdge(FI, MST)` / `MoreEdge()` / `IsBadFace()` (`cxx:88-95`): the
//!    edge iterator of `HLRBRep_Data` does the face-side rejections.
//! 2. Per edge: collect the interferences into two sorted lists, `ILHidden`
//!    (`TopAbs_IN`) and `ILOn` (`TopAbs_ON`), remembering `hasOut` for a
//!    rejected interference above a simple face (`cxx:114-139`).
//! 3. Merge the `SegIndex` pairs `(+n, -n)` of `ILHidden` (`cxx:141-227`).
//! 4. Resolve `ILOn`: complex transitions first, then the `State3D` table of
//!    `cxx:245-381` that moves `ON` interferences to `IN`/`OUT`.
//! 5. If both lists are empty and `hasOut` is false, a single
//!    `Compare` decides whether the whole edge is hidden (`cxx:386-395`).
//! 6. Otherwise the `ILHidden` pass: complex transitions, hiding start level,
//!    the "excess interferences" filter of `cxx:419-462`, the multi-inside
//!    filter of `cxx:468-511`, then the `EdgeBuilder` over `Builds(IN)` and
//!    `Builds(ON)` with the `allInt` clamp and a `Classify` confirmation.
//! 7. Then the same two passes over `ILOn` (`cxx:690-836`), without the
//!    `allInt` clamp and without `Classify`.
//!
//! Deviations:
//! * `myDS` (`hxx:52`) is a `handle<HLRBRep_Data>`; the port's hider is
//!   stateless and takes `&mut HlrBRepData`, so the caller keeps ownership
//!   (`HlrBRepInternalAlgo::my_ds`). The OCCT constructor (`cxx:30-33`) only
//!   stores the handle, and the single `Cache` instance of `HideSelected` is
//!   reused across faces - a stateless port is equivalent.
//! * `try { ... } catch (Standard_Failure)` (`cxx:104-850`) becomes
//!   `Result`: the `Err` arm of [`HlrBRepHider::hide_edge`] is the caught
//!   exception and the loop continues with the next edge.
//! * The `#ifdef OCCT_DEBUG` printing (`cxx:432-439`, `:259-262`, `:279-282`,
//!   `:319-321`, `:341-343`, `:365-367`, `:374-376`, `:842-849`) is not
//!   ported.
//! * In the `INTERNAL` / `stbef == UNKNOWN` arm (`cxx:374-379`) OCCT leaves the
//!   `ILOn` cursor untouched, which would spin; the port ends the `ILOn` pass
//!   there. `EdgeState` (`HLRBRep_Data.cxx:1580-1656`) only ever writes
//!   `IN`/`OUT`/`ON`, so the arm is unreachable.

use occt_topo::abs::{Orientation, State};
use crate::hlr_algo::HlrAlgoInterference;
use crate::hlrbrep_data::HlrBRepData;
use crate::hlrbrep_edge_builder::HlrBRepEdgeBuilder;
use crate::hlrbrep_edge_i_list::{add_interference, process_complex};
use crate::hlrbrep_edge_interference_tool::HlrBRepEdgeInterferenceTool;
use crate::hlrbrep_vertex_list::HlrBRepVertexList;
use crate::hlrtopobrep::out_liner::ShapeToTopAdaptorToolMap;

/// `HLRBRep_Hider` (`hxx:31-53`), stateless in the port.
#[derive(Debug, Clone, Default)]
pub struct HlrBRepHider;

impl HlrBRepHider {
    /// `HLRBRep_Hider(DS)` (`cxx:30-33`): only stores the data handle, which
    /// the port passes per call.
    pub fn new() -> Self {
        Self
    }

    /// `OwnHiding(FI)` (`cxx:37`): an empty body in OCCT.
    pub fn own_hiding(&self, _fi: i32) {}

    /// `Hide(FI, MST)` (`cxx:42-852`).
    ///
    /// `InitEdge` raises the hide level by two, so a face is never run twice on
    /// the same edges.
    pub fn hide(&mut self, ds: &mut HlrBRepData, fi: i32, mst: &mut ShapeToTopAdaptorToolMap) {
        ds.init_edge(fi, mst);
        if !ds.more_edge() {
            // there is nothing to do
            return;
        }
        if ds.is_bad_face() {
            return;
        }
        // `HLRBRep_EdgeInterferenceTool EIT(myDS)` (`cxx:97`).
        let mut eit = HlrBRepEdgeInterferenceTool::new();

        while ds.more_edge() {
            // loop on the Edges
            let e = ds.edge();
            // `try { ... } catch (Standard_Failure)` (`cxx:104-850`)
            let _ = Self::hide_edge(ds, e, &mut eit);
            // `for (...; myDS->NextEdge())` - `NextEdge()` defaults to
            // `skip = true` (`HLRBRep_Data.hxx:141`), so the increment
            // advances the iterator before re-testing.
            ds.next_edge(true);
        }
    }

    /// The body of the per-edge `try` block (`cxx:106-839`). `Err` is the
    /// caught `Standard_Failure`.
    fn hide_edge(
        ds: &mut HlrBRepData,
        e: i32,
        eit: &mut HlrBRepEdgeInterferenceTool,
    ) -> Result<(), &'static str> {
        let mut has_out = false;
        let mut il_hidden: Vec<HlrAlgoInterference> = Vec::new();
        let mut il_on: Vec<HlrAlgoInterference> = Vec::new();
        eit.load_edge(ds);

        // intersections with face-edges (`cxx:112-139`)
        ds.init_interference();
        while ds.more_interference() {
            if ds.rejected_interference() {
                if ds.above_interference() && ds.simple_hiding_face() {
                    has_out = true;
                }
            } else {
                let int = *ds.interference();
                match int.intersection().state() {
                    State::In => add_interference(&mut il_hidden, &int, eit),
                    State::On => add_interference(&mut il_on, &int, eit),
                    State::Out | State::Unknown => {}
                }
            }
            ds.next_interference();
        }

        // merge the SegIndex pairs (+n, -n) of ILHidden (`cxx:141-227`)
        loop {
            let mut modif = false;
            let mut i1 = 0;
            while i1 < il_hidden.len() && !modif {
                let numseg1 = il_hidden[i1].intersection().seg_index();
                if numseg1 != 0 {
                    let mut i2 = 0;
                    while i2 < il_hidden.len() && !modif {
                        let numseg2 = il_hidden[i2].intersection().seg_index();
                        if numseg1 + numseg2 == 0 {
                            let (stbef1, staft1) = il_hidden[i1].boundary().state_3d();
                            let (stbef2, staft2) = il_hidden[i2].boundary().state_3d();
                            if il_hidden[i1].orientation() == il_hidden[i2].orientation()
                                && il_hidden[i1].transition() == il_hidden[i2].transition()
                                && stbef1 == stbef2
                                && staft1 == staft2
                                && stbef1 != State::On
                                && staft1 != State::On
                            {
                                let idx1 = il_hidden[i1].intersection().index();
                                let idx2 = il_hidden[i2].intersection().index();
                                let mut nind: i32 = -1;
                                if idx1 != 0 {
                                    nind = idx1;
                                }
                                if idx2 != 0 {
                                    if nind != -1 {
                                        if idx1 != idx2 {
                                            nind = -1;
                                        }
                                    } else {
                                        nind = idx2;
                                    }
                                }
                                if idx1 == 0 && idx2 == 0 {
                                    nind = 0;
                                }

                                if nind != -1 {
                                    let p1 = il_hidden[i1].intersection().parameter();
                                    let p2 = il_hidden[i2].intersection().parameter();
                                    il_hidden[i1].change_intersection().set_seg_index(nind);
                                    il_hidden[i1]
                                        .change_intersection()
                                        .set_parameter((p1 + p2) * 0.5);
                                    il_hidden[i1].set_boundary_transition(Orientation::External);
                                    il_hidden.remove(i2);
                                    modif = true;
                                }
                            }
                        }
                        if !modif {
                            i2 += 1;
                        }
                    }
                }
                if !modif {
                    i1 += 1;
                }
            }
            if !modif {
                break;
            }
        }

        // process the interferences on ILOn (`cxx:230-384`)
        if !il_on.is_empty() {
            // complex transition on ILOn
            process_complex(&mut il_on, ds, eit)?;

            let mut it = 0;
            while it < il_on.len() {
                let int = il_on[it];
                // read the 3d states
                let (stbef, staft) = int.boundary().state_3d();
                match int.transition() {
                    Orientation::Forward => match staft {
                        State::Out => {
                            il_on.remove(it);
                        }
                        State::In => {
                            add_interference(&mut il_hidden, &int, eit);
                            il_on.remove(it);
                        }
                        State::Unknown | State::On => it += 1,
                    },
                    Orientation::Reversed => match stbef {
                        State::Out => {
                            il_on.remove(it);
                        }
                        State::In => {
                            add_interference(&mut il_hidden, &int, eit);
                            il_on.remove(it);
                        }
                        State::Unknown | State::On => it += 1,
                    },
                    Orientation::External => {
                        il_on.remove(it);
                    }
                    Orientation::Internal => match stbef {
                        State::In => match staft {
                            State::In => {
                                add_interference(&mut il_hidden, &int, eit);
                                il_on.remove(it);
                            }
                            State::On => {
                                // FORWARD in ILOn, REVERSED in ILHidden
                                il_on[it].set_transition(Orientation::Forward);
                                let mut rev = int;
                                rev.set_transition(Orientation::Reversed);
                                add_interference(&mut il_hidden, &rev, eit);
                                it += 1;
                            }
                            State::Out => {
                                // set REVERSED
                                il_on[it].set_transition(Orientation::Reversed);
                                let rev = il_on[it];
                                add_interference(&mut il_hidden, &rev, eit);
                                il_on.remove(it);
                            }
                            State::Unknown => it += 1,
                        },
                        State::On => {
                            match staft {
                                State::In => {
                                    // REVERSED in ILOn, FORWARD in ILHidden
                                    il_on[it].set_transition(Orientation::Reversed);
                                    let mut fwd = int;
                                    fwd.set_transition(Orientation::Forward);
                                    add_interference(&mut il_hidden, &fwd, eit);
                                }
                                State::On => {}
                                State::Out => {
                                    il_on[it].set_transition(Orientation::Reversed);
                                }
                                State::Unknown => {}
                            }
                            it += 1;
                        }
                        State::Out => match staft {
                            State::In => {
                                // set FORWARD
                                il_on[it].set_transition(Orientation::Forward);
                                let fwd = il_on[it];
                                add_interference(&mut il_hidden, &fwd, eit);
                                il_on.remove(it);
                            }
                            State::On => {
                                // FORWARD in ILOn
                                il_on[it].set_transition(Orientation::Forward);
                                it += 1;
                            }
                            State::Out => {
                                il_on.remove(it);
                            }
                            State::Unknown => it += 1,
                        },
                        State::Unknown => {
                            // OCCT leaves the cursor untouched here
                            // (`cxx:374-379`); the arm is unreachable because
                            // `EdgeState` only writes IN/OUT/ON.
                            it = il_on.len();
                        }
                    },
                }
            }
        }

        if il_hidden.is_empty() && il_on.is_empty() && !has_out {
            // classification (`cxx:386-395`)
            let st = ds.compare(e);
            if st == State::In || st == State::On {
                ds.e_data_mut(e).status_mut().hide_all();
            }
            return Ok(());
        }

        let mut found_hidden = false;

        if !il_hidden.is_empty() {
            // complex transition on ILHidden
            process_complex(&mut il_hidden, ds, eit)?;
            let mut level: i32 = 0;
            if !ds.simple_hiding_face() {
                // Level at Start
                level = ds.hiding_start_level(e, &il_hidden);
            }

            if ds.simple_hiding_face() {
                // remove excess interferences (`cxx:419-462`)
                let mut to_remove: Vec<f64> = Vec::new();
                let mut prev_trans = Orientation::External;
                let mut prev_param = 0.0f64;
                for int in il_hidden.iter() {
                    let a_trans = int.transition();
                    if a_trans == prev_trans {
                        if a_trans == Orientation::Forward {
                            to_remove.push(int.intersection().parameter());
                        } else if a_trans == Orientation::Reversed {
                            to_remove.push(prev_param);
                        }
                    }
                    prev_trans = a_trans;
                    prev_param = int.intersection().parameter();
                }
                let mut it = 0;
                while it < il_hidden.len() {
                    let a_param = il_hidden[it].intersection().parameter();
                    let mut found = false;
                    for k in 0..to_remove.len() {
                        if a_param == to_remove[k] {
                            found = true;
                            il_hidden.remove(it);
                            to_remove.remove(k);
                            break;
                        }
                    }
                    if !found {
                        it += 1;
                    }
                }
            }

            // suppress multi-inside intersections (`cxx:468-511`)
            let mut it = 0;
            while it < il_hidden.len() {
                match il_hidden[it].transition() {
                    Orientation::Forward => {
                        let decal = il_hidden[it].intersection().level();
                        if level > 0 {
                            il_hidden.remove(it);
                        } else {
                            it += 1;
                        }
                        level += decal;
                    }
                    Orientation::Reversed => {
                        level -= il_hidden[it].intersection().level();
                        if level > 0 {
                            il_hidden.remove(it);
                        } else {
                            it += 1;
                        }
                    }
                    _ => it += 1,
                }
            }

            if il_hidden.is_empty() {
                // edge hidden
                ds.e_data_mut(e).status_mut().hide_all();
            } else {
                found_hidden = true;
            }
        }

        if !il_hidden.is_empty() {
            // `cxx:524-689`: `allInt` and `IsSuspicion` are both live here
            let mut pmin = f64::MAX;
            let mut pmax = -pmin;
            let mut all_int = false;

            if il_hidden.len() > 1 {
                all_int = true;
                for int in il_hidden.iter() {
                    let p = int.intersection().parameter();
                    all_int = all_int && int.transition() == Orientation::Internal;
                    if p < pmin {
                        pmin = p;
                    }
                    if p > pmax {
                        pmax = p;
                    }
                }
            }

            let clamp = if all_int { Some((pmin, pmax)) } else { None };
            // build hidden parts
            Self::build_parts(ds, e, eit, &il_hidden, State::In, false, false, clamp, true)?;
            // build parts under the boundary
            Self::build_parts(ds, e, eit, &il_hidden, State::On, false, true, clamp, true)?;
        }

        if !il_on.is_empty() {
            let mut level: i32 = 0;
            if !ds.simple_hiding_face() {
                // Level at Start
                level = ds.hiding_start_level(e, &il_on);
            }
            if level > 0 {
                // suppress multi-inside intersections (`cxx:701-741`)
                let mut it = 0;
                while it < il_on.len() {
                    match il_on[it].transition() {
                        Orientation::Forward => {
                            let decal = il_on[it].intersection().level();
                            if level > 0 {
                                il_on.remove(it);
                            } else {
                                it += 1;
                            }
                            level += decal;
                        }
                        Orientation::Reversed => {
                            level -= il_on[it].intersection().level();
                            if level > 0 {
                                il_on.remove(it);
                            } else {
                                it += 1;
                            }
                        }
                        _ => it += 1,
                    }
                }
                if il_on.is_empty() && !found_hidden {
                    // edge hidden
                    ds.e_data_mut(e).status_mut().hide_all();
                }
            }
        }

        if !il_on.is_empty() {
            // `cxx:747-836`: no `allInt` clamp and no `Classify` here
            // build parts on the Face
            Self::build_parts(ds, e, eit, &il_on, State::In, true, false, None, false)?;
            // build hidden parts under the boundary
            Self::build_parts(ds, e, eit, &il_on, State::On, true, true, None, false)?;
        }

        Ok(())
    }

    /// One `while (EB.MoreEdges())` loop of `cxx:563-633` (`ILHidden`) or
    /// `cxx:751-793` / `cxx:797-834` (`ILOn`).
    ///
    /// `clamp` is the `if (allInt) { p1 = Max(p1, pmin); p2 = Min(p2, pmax); }`
    /// arm of `cxx:599-610`, present only for `ILHidden`; `classify` is the
    /// `if (IsSuspicion) aTestState = Classify(...)` gate of `cxx:611-620`,
    /// also only for `ILHidden` (the `ILOn` loops hide unconditionally).
    #[allow(clippy::too_many_arguments)]
    fn build_parts(
        ds: &mut HlrBRepData,
        e: i32,
        eit: &HlrBRepEdgeInterferenceTool,
        il: &[HlrAlgoInterference],
        to_build: State,
        on_face: bool,
        on_boundary: bool,
        clamp: Option<(f64, f64)>,
        classify: bool,
    ) -> Result<(), &'static str> {
        let mut vlist = HlrBRepVertexList::new(eit, il);
        let mut eb = HlrBRepEdgeBuilder::new(&mut vlist)?;
        eb.builds(to_build);
        while eb.more_edges() {
            let mut p1 = 0.0f64;
            let mut p2 = 0.0f64;
            let mut tol1 = 0.0f32;
            let mut tol2 = 0.0f32;
            let mut mask = 0i32;
            while eb.more_vertices() {
                match eb.orientation() {
                    Orientation::Forward => {
                        let cur = eb.current().ok_or("HLRBRep_Hider: no current vertex")?;
                        p1 = cur.parameter();
                        tol1 = cur.tolerance();
                        mask |= 1;
                    }
                    Orientation::Reversed => {
                        let cur = eb.current().ok_or("HLRBRep_Hider: no current vertex")?;
                        p2 = cur.parameter();
                        tol2 = cur.tolerance();
                        mask |= 2;
                    }
                    Orientation::Internal | Orientation::External => {}
                }
                eb.next_vertex().ok_or("HLRBRep_Hider: no current edge")?;
            }

            if mask != 3 || p2 - p1 <= 1.0e-7 {
                eb.next_edge();
                continue;
            }

            if let Some((pmin, pmax)) = clamp {
                if p1 < pmin {
                    p1 = pmin;
                }
                if p2 > pmax {
                    p2 = pmax;
                }
            }

            if classify {
                let mut tmplevel = 0;
                let a_test_state = ds.classify(e, true, &mut tmplevel, (p1 + p2) / 2.);
                if a_test_state == State::Out {
                    eb.next_edge();
                    continue;
                }
            }

            ds.e_data_mut(e)
                .status_mut()
                .hide(p1, tol1, p2, tol2, on_face, on_boundary);
            eb.next_edge();
        }
        eb.destroy();
        Ok(())
    }
}
