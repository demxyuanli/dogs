//! `BRepFill_Sweep::EvalExtrapol` (`BRepFill_Sweep.cxx:3910-3978`) and
//! `BRepFill_Sweep::PerformCorner` (`cxx:3587-3906`).
//!
//! Both are ported end to end: the bisector frame construction
//! (`cxx:3611-3697`), `BRepFill_TrimShellCorner` with its `Modified`
//! write-back (`cxx:3715-3766`), the `IsUClosed` seam glue (`cxx:3771-3779`)
//! and the `Round` filling block (`cxx:3781-3887`).
//!
//! Translation notes:
//! - `gp_Mat::Column(3)` is [`GpMat::column`] (1-based) wrapped by
//!   [`GpVec::from_xyz`].
//! - `Bnd_Box` is a plain min/max over the section law's poles (`Box`,
//!   `cxx:193-206`).
//! - `M_PI` is [`std::f64::consts::PI`].
//! - `TopExp::CommonVertex` is [`crate::topo_tools::common_vertex`].

use occt_core::gp::{GpAx2, GpDir, GpVec};

use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::is_same;

use super::filling;
use super::helpers::{self, ShapeTable};
use super::{BrepFillSweep, TransitionStyle};

/// `anAngularTol` of `PerformCorner` (`cxx:3594`).
const CORNER_ANGULAR_TOL: f64 = 0.025;

impl BrepFillSweep {
    /// `EvalExtrapol(Index, Transition)` (`cxx:3910-3978`): the length the
    /// first / last sweep surface must be extended for a `Right` transition.
    pub fn eval_extrapol(&self, index: usize, transition: TransitionStyle) -> Result<f64, String> {
        let mut extrap = 0.0;
        if transition != TransitionStyle::Right {
            return Ok(extrap);
        }
        let nb_law = self.my_loc.base().nb_law();
        let i1;
        let i2;
        if index == 1 || index == nb_law + 1 {
            if !self.my_loc.base().is_closed() || !self.my_sec.base().is_v_closed() {
                return Ok(extrap);
            }
            i1 = nb_law;
            i2 = 1;
        } else {
            i1 = index - 1;
            i2 = index;
        }

        let law1 = self
            .my_loc
            .base()
            .law(i1)
            .ok_or_else(|| format!("EvalExtrapol: null location law {i1}"))?;
        let law2 = self
            .my_loc
            .base()
            .law(i2)
            .ok_or_else(|| format!("EvalExtrapol: null location law {i2}"))?;
        let (f, l) = law1.get_domain();
        let _ = f;
        let (m1, _v1) = law1.d0(l)?;
        let t1 = GpVec::from_xyz(&m1.column(3));
        let (f, _l) = law2.get_domain();
        let (m2, _v2) = law2.d0(f)?;
        let t2 = GpVec::from_xyz(&m2.column(3));

        let alpha = angle_between(&t1, &t2);
        if alpha > self.my_ang_max || alpha < self.my_ang_min {
            return Ok(extrap);
        }

        let sec = self
            .my_sec
            .concatened_law()
            .ok_or_else(|| "EvalExtrapol: null ConcatenedLaw (BRepFill_Sweep.cxx:3952)".to_string())?;

        // `cxx:3954-3961`: the parameter of the section at the corner.
        let (_lf, length) = self.my_loc.base().curvilinear_bounds(nb_law)?;
        let law = self
            .my_sec
            .base().law(1)
            .ok_or_else(|| "EvalExtrapol: null section law 1".to_string())?;
        let (sec_first, sec_last) = law.get_domain();
        let sec_len = sec_last - sec_first;
        let (_lf, ll) = self.my_loc.base().curvilinear_bounds(i1)?;
        let u = sec_first + (ll / length) * sec_len;

        // `Box(Sec, U, box)` (`cxx:3962-3965`).
        let (nb_poles, _, _) = sec.section_shape();
        let (poles, _w) = sec.d0(u)?;
        let mut xmin = f64::INFINITY;
        let mut ymin = f64::INFINITY;
        let mut zmin = f64::INFINITY;
        let mut xmax = f64::NEG_INFINITY;
        let mut ymax = f64::NEG_INFINITY;
        let mut zmax = f64::NEG_INFINITY;
        for p in poles.iter().take(nb_poles.max(0) as usize) {
            xmin = xmin.min(p.x());
            ymin = ymin.min(p.y());
            zmin = zmin.min(p.z());
            xmax = xmax.max(p.x());
            ymax = ymax.max(p.y());
            zmax = zmax.max(p.z());
        }

        // `cxx:3967-3977`.
        let mut r = xmin.abs().max(xmax.abs()).max(ymin.abs()).max(ymax.abs());
        let coef = 2.0;
        r *= coef;
        extrap = zmin.abs().max(zmax.abs()) + 100.0 * self.my_tol_3d;
        extrap += r * (alpha / 2.0).tan();
        Ok(extrap)
    }

    /// `PerformCorner(Index, Transition, Bounds)` (`cxx:3587-3906`).
    pub fn perform_corner(
        &mut self,
        index: usize,
        transition: TransitionStyle,
        bounds: &ShapeTable,
    ) -> Result<bool, String> {
        // `cxx:3591-3594`.
        if transition == TransitionStyle::Modified {
            return Ok(true);
        }

        // `cxx:3611-3618`.
        let (i1, i2) = if index > 1 {
            (index - 1, index)
        } else {
            (self.my_loc.base().nb_law(), 1)
        };

        // `cxx:3620-3628`.
        let law1 = self
            .my_loc
            .base()
            .law(i1)
            .ok_or_else(|| format!("PerformCorner: null location law {i1}"))?;
        let law2 = self
            .my_loc
            .base()
            .law(i2)
            .ok_or_else(|| format!("PerformCorner: null location law {i2}"))?;
        let (_f, l) = law1.get_domain();
        let curve1 = law1
            .get_curve()
            .ok_or_else(|| "PerformCorner: null location curve".to_string())?;
        let (p1, mut t1) = curve1.d1(l);
        t1.normalize();
        let (f, _l) = law2.get_domain();
        let curve2 = law2
            .get_curve()
            .ok_or_else(|| "PerformCorner: null location curve".to_string())?;
        let (p2, mut t2) = curve2.d1(f);
        t2.normalize();

        // `cxx:3630-3651`: the two laws are tangent, so there is no corner.
        let mut is_tangent = false;
        let mut sortant_t = GpVec::zero();
        if angle_between(&t1, &t2) < self.my_ang_min {
            is_tangent = true;
            let (m1, _v) = law1.d0(l)?;
            let t1_axis = GpVec::from_xyz(&m1.column(3));
            let (m2, _v) = law2.d0(l)?;
            let t2_axis = GpVec::from_xyz(&m2.column(3));
            if angle_between(&t1_axis, &t2_axis) < self.my_ang_min {
                return Ok(true);
            }
            // `cxx:3649`.
            sortant_t = t2_axis.subtracted(&t1_axis);
        }

        // `cxx:3654-3658`.
        let angle = angle_between(&t1, &t2);
        if angle >= std::f64::consts::PI - CORNER_ANGULAR_TOL {
            return Ok(false);
        }
        // `cxx:3656-3659`: a `Right` corner wider than `myAngMax` becomes a
        // `Round` one.
        let mut the_transition = transition;
        if the_transition == TransitionStyle::Right && angle > self.my_ang_max {
            the_transition = TransitionStyle::Round;
        }

        // `cxx:3660-3678`.
        let mut tang = t1.added(&t2); // `Tang = T1 + T2`
        let normal_of_bis_plane = GpDir::new(tang.x(), tang.y(), tang.z())
            .map_err(|e| format!("PerformCorner: gp_Dir(Tang) is degenerate: {e} (cxx:3663)"))?;
        let intersect_cross_dir = t1.crossed(&t2);
        let sortant: GpVec;
        if is_tangent {
            // `Sortant -= Tang.Dot(Tang) * Tang` (`cxx:3667`).
            let d = tang.dot(&tang);
            sortant = sortant_t.subtracted(&tang.multiplied_scalar(d));
        } else {
            // `cxx:3670-3673`: `Sortant = T2 - T1`, negated, and `Tang` is
            // made orthogonal to `T2`.
            sortant = t2.subtracted(&t1).reversed();
            tang = tang.subtracted(&t2.multiplied_scalar(tang.dot(&t2)));
        }

        // `cxx:3675` — `P1.BaryCenter(0.5, P2, 0.5)`.
        let bary = p1.xyz().multiplied(0.5);
        let p2xyz = p2.xyz().multiplied(0.5);
        let bary_pnt = occt_core::gp::GpPnt::from_xyz(&bary.added(&p2xyz));

        // `cxx:3676-3680`.
        let n = GpDir::new(sortant.x(), sortant.y(), sortant.z())
            .map_err(|e| format!("PerformCorner: gp_Dir(Sortant) is degenerate: {e} (cxx:3676)"))?;
        let dx = GpDir::new(tang.x(), tang.y(), tang.z())
            .map_err(|e| format!("PerformCorner: gp_Dir(Tang) is degenerate: {e} (cxx:3677)"))?;
        let axe = GpAx2::new(bary_pnt, n, dx)
            .map_err(|e| format!("PerformCorner: gp_Ax2(P1, N, Dx) has parallel axes: {e} (cxx:3679)"))?;
        let axe_of_bis_plane = GpAx2::new(bary_pnt, normal_of_bis_plane, normal_of_bis_plane)
            .map_err(|e| format!("PerformCorner: gp_Ax2(P1, NormalOfBisPlane): {e} (cxx:3680)"))?;

        // `cxx:3683-3697` — `aFaces` / `aUEdges` are the two columns the trim
        // reads; the `UEdges` copy at `cxx:3683-3686` is never read afterwards.
        let nb_law = self.my_sec.base().nb_law();
        let mut a_faces = ShapeTable::new(self.my_faces.nb_rows(), 2);
        helpers::translate(&self.my_faces, i1, &mut a_faces, 1);
        helpers::translate(&self.my_faces, i2, &mut a_faces, 2);
        let mut a_u_edges = ShapeTable::new(self.my_u_edges.nb_rows(), 2);
        helpers::translate(&self.my_u_edges, i1, &mut a_u_edges, 1);
        helpers::translate(&self.my_u_edges, i2, &mut a_u_edges, 2);

        // `cxx:3699-3712` builds `aPlaneF` from `gp_Pln(P1, gp_Dir(T2 + T1))`
        // and raises its tolerance to `10 * Precision::Confusion()`. In 7.9.3
        // no later statement reads `aPlaneF`, so it has no port counterpart.

        // `cxx:3715-3716` — `BRepFill_TrimShellCorner aTrim(aFaces, Transition,
        // AxeOfBisPlane, anIntersectPointCrossDirection)` plus
        // `AddBounds(Bounds)` / `AddUEdges(aUEdges)` / `AddVEdges(myVEdges, Index)`.
        let (trim_done, trim_has_section) = self.trim_shell_corner(
            &a_faces,
            the_transition,
            &axe_of_bis_plane,
            &intersect_cross_dir,
            &a_u_edges,
            bounds,
            index,
            i1,
            i2,
        )?;

        // `cxx:3767-3770`: a failed trim leaves the data untouched for a
        // `Right` transition and for a trim that found a section.
        if !trim_done && (the_transition == TransitionStyle::Right || trim_has_section) {
            return Ok(true);
        }

        // `cxx:3771-3779` — `mySec->IsUClosed()` seam glue, for any transition.
        if self.my_sec.base().is_u_closed() {
            for col in [i1, i2] {
                let v = self.my_u_edges.value(nb_law + 1, col).cloned();
                self.my_u_edges.set_value(1, col, v);
            }
        }

        // `cxx:3781-3887`.
        if the_transition == TransitionStyle::Round {
            self.corner_round_filling(i1, i2, &axe, &t1, bounds)?;
        }

        Ok(true)
    }

    /// `BRepFill_TrimShellCorner` construction and `Perform` call
    /// (`BRepFill_Sweep.cxx:3715-3719`), the `Modified` write-back
    /// (`cxx:3721-3766`), then answering the reference's `IsDone()` and
    /// `HasSection()`.
    ///
    /// The class itself lives in [`super::trim_shell_corner`].
    fn trim_shell_corner(
        &mut self,
        a_faces: &ShapeTable,
        transition: TransitionStyle,
        axe_of_bis_plane: &GpAx2,
        cross_dir: &GpVec,
        a_u_edges: &ShapeTable,
        bounds: &ShapeTable,
        index: usize,
        i1: usize,
        i2: usize,
    ) -> Result<(bool, bool), String> {
        let mut trim =
            super::trim_shell_corner::TrimShellCorner::new(
                a_faces,
                transition,
                *axe_of_bis_plane,
                cross_dir.clone(),
            );
        trim.add_bounds(bounds);
        trim.add_u_edges(a_u_edges);
        trim.add_v_edges(&self.my_v_edges, index);
        trim.perform()?;
        let done = trim.is_done();
        let has_section = trim.has_section();

        // `cxx:3721-3739`: a successful trim replaces every `v`-iso edge of the
        // corner column, dropping the ones it did not keep.
        if done {
            let nb_law = self.my_sec.base().nb_law();
            for ii in 1..=nb_law {
                let current = self.my_v_edges.value(ii, index).cloned();
                let modified = match &current {
                    Some(s) => trim.modified(s),
                    None => Vec::new(),
                };
                self.my_v_edges
                    .set_value(ii, index, modified.first().cloned());
            }

            // `cxx:3740-3766`: the two face columns and their `u`-iso edges.
            for col in [i1, i2] {
                for ii in 1..=nb_law {
                    if let Some(f) = self.my_faces.value(ii, col).cloned() {
                        let modified = trim.modified(&f);
                        if let Some(s) = modified.first() {
                            self.my_faces.set_value(ii, col, Some(s.clone()));
                        }
                    }
                }
                for row in 1..=self.my_u_edges.nb_rows() {
                    if let Some(e) = self.my_u_edges.value(row, col).cloned() {
                        let modified = trim.modified(&e);
                        if let Some(s) = modified.first() {
                            self.my_u_edges.set_value(row, col, Some(s.clone()));
                        }
                    }
                }
            }
        }

        Ok((done, has_section))
    }

    /// The `Round` filling block of `PerformCorner` (`cxx:3781-3887`),
    /// reached after `BRepFill_TrimShellCorner` has rewritten
    /// `myVEdges` / `myFaces` / `myUEdges`.
    ///
    /// `t1` is the reference's `T1` (the normalized tangent at the end of
    /// location law `I1`) and `axe` its `Axe` (`cxx:3679`).
    fn corner_round_filling(
        &mut self,
        i1: usize,
        i2: usize,
        axe: &GpAx2,
        t1: &GpVec,
        bounds: &ShapeTable,
    ) -> Result<(), String> {

        // `cxx:3781-3787`.
        let mut bord1: Option<Edge> = None;
        let mut bord2: Option<Edge> = None;
        let mut bord_first: Option<Edge> = None;
        let mut has_filling = false;

        let nb_law = self.my_sec.base().nb_law();
        let tol3d = self.my_tol_3d;
        // `Bounds` (`cxx:3795-3796`): the per-section limiting wires.

        for ii in 1..=nb_law {
            // `cxx:3794-3795` — `KeepEdge(myFaces->Value(ii, I), Bounds->Value(ii, J), list)`.
            let mut list1 = Vec::new();
            let mut list2 = Vec::new();
            if let (Some(face), Some(bound)) = (self.my_faces.value(ii, i1), bounds.value(ii, 1)) {
                list1 = helpers::keep_edge(face, &Edge(bound.clone()));
            }
            if let (Some(face), Some(bound)) = (self.my_faces.value(ii, i2), bounds.value(ii, 2)) {
                list2 = helpers::keep_edge(face, &Edge(bound.clone()));
            }

            // `cxx:3797` — only symmetric free borders are filled.
            if list1.len() != list2.len() {
                continue;
            }

            for (v1, v2) in list1.iter().zip(list2.iter()) {
                let e = Edge(v1.clone());

                // `cxx:3803-3845`: the transversal choice of constraints.
                if has_filling {
                    let mut e1 = None;
                    let mut e2 = None;
                    let vertices = crate::topo_tools::edge_vertices(&e);
                    let vf = vertices.0.clone();
                    for bord in [bord1.as_ref(), bord2.as_ref(), bord_first.as_ref()] {
                        let Some(bord) = bord else { continue };
                        let Some(vc) = crate::topo_tools::common_vertex(&e, bord) else { continue };
                        let is_first = match vf.as_ref() {
                            Some(vf) => is_same(&vc, &vf.0),
                            None => false,
                        };
                        if is_first {
                            e1 = Some(bord.clone());
                        } else {
                            e2 = Some(bord.clone());
                        }
                    }
                    bord1 = e1;
                    bord2 = e2;
                }

                // `cxx:3849-3860`.
                let face1 = self
                    .my_faces
                    .value(ii, i1)
                    .cloned()
                    .ok_or_else(|| format!("PerformCorner: null face ({ii}, {i1}) before Filling (cxx:3850)"))?;
                let face2 = self
                    .my_faces
                    .value(ii, i2)
                    .cloned()
                    .ok_or_else(|| format!("PerformCorner: null face ({ii}, {i2}) before Filling (cxx:3852)"))?;
                let fr = filling::filling(
                    &e,
                    &Face(face1),
                    &Edge(v2.clone()),
                    &Face(face2),
                    &mut self.my_v_edges_modified,
                    tol3d,
                    axe,
                    t1,
                    bord1.clone(),
                    bord2.clone(),
                )?;
                // `Filling`'s `Aux1` / `Aux2` are in/out: the reference keeps
                // the updated values in `Bord1` / `Bord2`.
                bord1 = fr.aux1;
                bord2 = fr.aux2;

                // `cxx:3861-3877`.
                if fr.done {
                    if let Some(ff) = fr.result {
                        self.my_aux_shape.push(ff.0.clone());

                        // `myVEdges->ChangeValue(ii, I2)` becomes a compound of
                        // the previous value and the filling face.
                        let prev = self.my_v_edges.value(ii, i2).cloned();
                        let mut members: Vec<TopoShape> = Vec::new();
                        if let Some(p) = prev {
                            members.push(p);
                        }
                        members.push(ff.0.clone());
                        let comp = TopoBuilder::new().make_compound_of(&members);
                        self.my_v_edges.set_value(ii, i2, Some(comp.0.clone()));

                        // `BB.Add(myTapes->ChangeValue(ii), FF)`.
                        let tape = self.my_tapes.get_mut(ii - 1).and_then(|t| t.as_mut());
                        match tape {
                            Some(t) => TopoBuilder::new().add(t, &ff.0),
                            None => {
                                return Err(format!(
                                    "PerformCorner: change value of myTapes({ii}) is null (BRepFill_Sweep.cxx:3875)"
                                ))
                            }
                        }
                        has_filling = true;
                    }
                }

                // `cxx:3879-3881`.
                if ii == 1 {
                    bord_first = bord1.clone();
                }
            }
        }

        Ok(())
    }
}

/// `gp_Vec::Angle` (`gp_Vec.cxx`).
fn angle_between(a: &GpVec, b: &GpVec) -> f64 {
    a.angle(b)
}
