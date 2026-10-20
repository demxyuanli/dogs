use occt_core::gp::GpTrsf2d;
use occt_core::precision::CONFUSION;

use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape};
use crate::shhealing::{first_vertex, fix_same_parameter, last_vertex, ShapeFixEdge};
use crate::tgeometry::GeometryRegistry;

use super::super::iges_boundary::connect_next_wire;
use super::super::iges_to_brep_curve_and_surface::Session;

use super::boundary::{connect_wire, rebuild_wire_max_tolerance};
use super::curve_on_face::circular_bad_case;
use super::pcurve_edge::reverse_curves_2d;
use super::pcurve_transfer::{curve_sequence, edges_from_2d, transfer_pcurve};

/// `IGESToBRep_IGESBoundary` for the `usescurve == true` entry point
/// (`IGESToBRep_IGESBoundary.cxx:94-113`), whose only caller is
/// `IGESToBRep_BRepEntity::TransferLoop` (`IGESToBRep_BRepEntity.cxx:415-421`).
/// `mysewd`, `mysewd3d`, and `mysewd2d` become `out`, `wire3`, and `wire2`;
/// `myCS` is the `Session`; `ConnectNextWire`'s `distmin` / `revsewd` /
/// `revnextsewd` are never read back by OCCT and are not kept.
pub(in crate::iges_read) struct WireDataBoundary {
    pub precision: f64,
    pub maxtol: f64,
    /// `Init`'s `filepreference` (`cxx:52-63`). `TransferLoop` passes 3.
    pub file_preference: i32,
    pub out: Vec<Edge>,
    pub wire3: Vec<Edge>,
    pub wire2: Vec<Edge>,
    pub ok: bool,
    pub ok3: bool,
    pub ok2: bool,
}

impl WireDataBoundary {
    pub(in crate::iges_read) fn new(precision: f64, maxtol: f64, file_preference: i32) -> Self {
        Self {
            precision,
            maxtol,
            file_preference,
            out: Vec::new(),
            wire3: Vec::new(),
            wire2: Vec::new(),
            ok: true,
            ok3: true,
            ok2: true,
        }
    }

    /// One `Transfer` call (`cxx:146-536`). `usescurve` is true, so the 3D side
    /// is the single prebuilt edge `edge3d` and `toreverse3d` is false
    /// (`IGESToBRep_IGESBoundary.cxx:105-111`).
    ///
    /// Returns the `lsewd` of this call.
    pub(in crate::iges_read) fn transfer(
        &mut self,
        session: &mut Session<'_>,
        face: &TopoShape,
        trans: &GpTrsf2d,
        u_fact: f64,
        edge3d: &[Edge],
        param_curves: &[Option<usize>],
        toreverse2d: bool,
        number: usize,
    ) -> Result<Vec<Edge>, String> {
        let len3d = edge3d.len();
        // `cxx:190-197`: every entry of `curves2d` is flattened into `seq2d`.
        let mut seq2d: Vec<usize> = Vec::new();
        for child in param_curves.iter().flatten() {
            seq2d.extend(curve_sequence(session.model, *child)?);
        }
        let len2d = seq2d.len();

        // `cxx:199-207` with `surfcurv == 0` (`read.surfacecurve.mode`).
        let mut translate3d = true;
        let mut translate2d = true;
        if len3d == 0 {
            translate3d = false;
        } else if len2d == 0 {
            translate2d = false;
        }
        // `cxx:210-233` with `surfcurv == 0`: the file preference decides.
        let mut preferred3d = true;
        let mut preferred2d = true;
        if translate3d && translate2d {
            if self.file_preference == 2 {
                preferred3d = false;
            } else if self.file_preference == 3 {
                preferred2d = false;
            } else {
                preferred3d = false;
            }
        }
        // `cxx:235-239`.
        if translate3d && translate2d && len3d != len2d {
            translate3d = preferred3d;
            translate2d = preferred2d;
        }

        let mut lsewd3d = edge3d.to_vec();
        let mut lsewd2d = Vec::new();

        if translate3d && !translate2d {
            // `cxx:240-265`, `usescurve` arm: `Gsewd->Add(scurve3d->Wire())`.
        } else if !translate3d && translate2d {
            // `cxx:266-281`.
            for child in param_curves.iter().flatten() {
                if let Ok(edges) = edges_from_2d(session, *child, face, trans, u_fact) {
                    lsewd2d.extend(edges);
                }
            }
            if toreverse2d {
                lsewd2d = reverse_curves_2d(&lsewd2d, face);
            }
        } else {
            // `cxx:282-437`: both representations with `len3d == 1`.
            // `TC.SetBadCase(false)` / `BadCase()` (`cxx:291-311`): the 508
            // caller hands over a prebuilt edge, so nothing reports the
            // closed-arc bad case for the 3D side.
            let mut lsewd3d_c = std::mem::take(&mut lsewd3d);
            self.ok3 &= connect_next_wire(
                &mut self.wire3,
                &mut lsewd3d_c,
                self.precision,
                self.maxtol,
            );
            lsewd3d = lsewd3d_c;
            let mut bad2d = false;
            let mut local2d = false;
            if let Some(&child2) = seq2d.first() {
                bad2d = circular_bad_case(session, child2);
                match edges_from_2d(session, child2, face, trans, u_fact) {
                    Ok(edges) if !edges.is_empty() => {
                        lsewd2d = edges;
                        if toreverse2d {
                            lsewd2d = reverse_curves_2d(&lsewd2d, face);
                        }
                        self.ok2 &= connect_next_wire(
                            &mut self.wire2,
                            &mut lsewd2d,
                            self.precision,
                            self.maxtol,
                        );
                        local2d = true;
                    }
                    _ => {}
                }
            }
            // `cxx:342-358`. `usescurve` clears `isBSpline`.
            let mut local3d = true;
            if local3d && local2d && (lsewd3d.len() != 1 || lsewd2d.len() != 1 || bad2d) {
                local3d = preferred3d;
                local2d = preferred2d;
            }
            if local3d && local2d {
                // `cxx:369-434`: copy the pcurve, `FixReversed2d`, then
                // `FixSameParameter`. Above `maxtol`, `SetTolerance` forces
                // `Precision::Confusion` and the pcurves are removed;
                // `Preferred3d` false keeps the 2D edge.
                let reg = GeometryRegistry::global();
                let face_fix = Face(face.clone());
                let mut use_2d = false;
                for (k, edge2) in lsewd2d.iter().enumerate() {
                    if k >= lsewd3d.len() {
                        break;
                    }
                    if !transfer_pcurve(edge2, &lsewd3d[k], face) {
                        continue;
                    }
                    let _ = ShapeFixEdge::new().fix_reversed_2d(&lsewd3d[k], face);
                    fix_same_parameter(&lsewd3d[k], &face_fix);
                    let maxdev = reg.edge_tolerance(&lsewd3d[k].0);
                    if maxdev > self.maxtol {
                        reg.set_edge_tolerance(&lsewd3d[k].0, CONFUSION);
                        if let Some(v) = first_vertex(&lsewd3d[k]) {
                            v.set_tolerance(CONFUSION);
                        }
                        if let Some(v) = last_vertex(&lsewd3d[k]) {
                            v.set_tolerance(CONFUSION);
                        }
                        if preferred3d {
                            for edge in lsewd3d.iter().take(k + 1) {
                                reg.remove_pcurves_on_surface(&edge.0, face);
                            }
                        } else {
                            use_2d = true;
                        }
                        break;
                    }
                }
                if use_2d {
                    lsewd3d = lsewd2d.clone();
                }
            } else if local2d {
                lsewd3d = lsewd2d.clone();
            }
        }

        let mut lsewd = lsewd3d.clone();
        // `cxx:437-439`.
        if !connect_next_wire(&mut self.out, &mut lsewd, self.precision, self.maxtol) {
            // `cxx:446-508`: re-transfer every segment, run
            // `ShapeFix_Wire::Perform` on each wire, and keep the one whose
            // maximum tolerance is smaller. `Gsewd2d` is then filled with the
            // 3D wire (`cxx:496`, a 3D wire in both slots).
            let mut rebuilt3d = edge3d.to_vec();
            let mut rebuilt2d = Vec::new();
            for &child in &seq2d {
                if let Ok(edges) = edges_from_2d(session, child, face, trans, u_fact) {
                    rebuilt2d.extend(edges);
                }
            }
            let builder = TopoBuilder::new();
            let face_fix = Face(face.clone());
            let mut wire3 = builder.make_wire(&rebuilt3d);
            let mut wire2 = builder.make_wire(&rebuilt2d);
            crate::shhealing::perform_wire_constructor_defaults(
                &mut wire3,
                &face_fix,
                self.precision,
            );
            crate::shhealing::perform_wire_constructor_defaults(
                &mut wire2,
                &face_fix,
                self.precision,
            );
            rebuilt3d = crate::topo_tools::edges_of_wire(&wire3);
            rebuilt2d = crate::topo_tools::edges_of_wire(&wire2);
            let heal3 = rebuild_wire_max_tolerance(&rebuilt3d);
            let heal2 = rebuild_wire_max_tolerance(&rebuilt2d);
            let keep = if heal3 < heal2 {
                rebuilt3d.clone()
            } else {
                rebuilt2d
            };
            self.wire3 = rebuilt3d;
            self.wire2 = self.wire3.clone();
            self.out.extend(keep.clone());
            self.ok = true;
            self.ok3 = true;
            self.ok2 = true;
            lsewd = keep;
        }

        if number > 1 {
            // `cxx:512-529`.
            self.ok &= connect_wire(
                &mut self.out,
                &mut lsewd,
                len3d > 1 || len2d > 1,
                self.precision,
                self.maxtol,
            );
            let mut w3 = std::mem::take(&mut self.wire3);
            self.ok3 &= connect_wire(
                &mut w3,
                &mut lsewd3d,
                len3d > 1,
                self.precision,
                self.maxtol,
            );
            self.wire3 = w3;
            let mut w2 = std::mem::take(&mut self.wire2);
            self.ok2 &= connect_wire(
                &mut w2,
                &mut lsewd2d,
                len2d > 1,
                self.precision,
                self.maxtol,
            );
            self.wire2 = w2;
        } else {
            // `cxx:530-535`.
            self.out = lsewd.clone();
            self.wire3 = lsewd3d.clone();
            self.wire2 = lsewd2d;
        }
        Ok(lsewd)
    }
}
