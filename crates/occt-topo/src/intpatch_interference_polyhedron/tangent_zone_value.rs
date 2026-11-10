//! `IntPatch_InterferencePolyhedron::TangentZoneValue`.
//!
//! Source: `IntPatch_InterferencePolyhedron.cxx:842-1117`.
//!
//! Called from `Intersect` when the six plane distances of the two triangles
//! are all within `Tolerance` and the planes really cross. It fills `TheTZ`
//! with the section points that bound the tangency zone:
//!
//! * the position of each point of `Tri1` in `Tri2`
//!   (`cxx:866-937`, vertex / edge / face)
//! * when `Tri1` is not fully inside `Tri2`, the position of each point of
//!   `Tri2` in `Tri1` (`cxx:941-1004`)
//! * otherwise the edge/edge crossings of the two triangles
//!   (`cxx:1019-1111`)
//!
//! The `TheTZ.Insert(...)` retry loop at `cxx:1081-1111` is transliterated even
//! though `Intf_TangentZone::Insert` always returns `false`
//! (`Intf_TangentZone.cxx:79-118`); it therefore falls through to the `Dump`
//! tail, which writes to `std::cout` and is not translated.

use occt_core::gp::GpVec;
use occt_core::intf::{contain, IntfPIType, IntfSectionPoint, IntfTangentZone};

use super::{IntPatchInterferencePolyhedron, POURCENT3};
use crate::intpatch_polyhedron::IntPatchPolyhedron;

impl IntPatchInterferencePolyhedron {
    /// `TangentZoneValue(TheTZ, FirstPol, Tri1, SeconPol, Tri2)`
    /// (`cxx:842-1117`).
    pub(crate) fn tangent_zone_value(
        &self,
        the_tz: &mut IntfTangentZone,
        first_pol: &IntPatchPolyhedron,
        tri1: i32,
        secon_pol: &IntPatchPolyhedron,
        tri2: i32,
    ) -> bool {
        let tolerance = self.base.tolerance;

        let mut finished = false;

        let mut t_op = [IntfPIType::External; 3];
        let mut t_tp = [IntfPIType::External; 3];

        let mut nbp_int = 0i32;
        let mut tpi: Vec<IntfSectionPoint> = Vec::new();

        // Compute the positions of the points of <Tri1> in the triangle <Tri2>
        // (`cxx:866-937`).
        for nob in 0..=2usize {
            let _nob2 = POURCENT3[nob + 1];
            for nou in 0..=2usize {
                let nou2 = POURCENT3[nou + 1];
                if self.dp_op_t[nob][nou] <= tolerance {
                    tpi.push(IntfSectionPoint::with_3d(
                        &first_pol.point(self.oi[nob]),
                        IntfPIType::Vertex,
                        self.oi[nob],
                        0,
                        0.0,
                        IntfPIType::Vertex,
                        self.ti[nou],
                        0,
                        0.0,
                        1.0,
                    ));
                    t_op[nob] = IntfPIType::Vertex;
                    t_tp[nou] = IntfPIType::Vertex;
                    nbp_int += 1;
                    break;
                } else if self.dp_oe_t[nob][nou].abs() <= tolerance {
                    if self.dp_op_t[nob][nou] + self.dp_op_t[nob][nou2]
                        < self.vtt[nou].modulus()
                    {
                        let mut par = self.dp_op_t[nob][nou]
                            / (self.dp_op_t[nob][nou] + self.dp_op_t[nob][nou2]);
                        if self.ti[nou] > self.ti[nou2] {
                            par = 1.0 - par;
                        }
                        tpi.push(IntfSectionPoint::with_3d(
                            &first_pol.point(self.oi[nob]),
                            IntfPIType::Vertex,
                            self.oi[nob],
                            0,
                            0.0,
                            IntfPIType::Edge,
                            self.ti[nou].min(self.ti[nou2]),
                            self.ti[nou].max(self.ti[nou2]),
                            par,
                            1.0,
                        ));
                        t_op[nob] = IntfPIType::Edge;
                        nbp_int += 1;
                        break;
                    }
                }
            }
            if t_op[nob] == IntfPIType::External {
                if contain(
                    &secon_pol.point(self.ti[0]),
                    &secon_pol.point(self.ti[1]),
                    &secon_pol.point(self.ti[2]),
                    &first_pol.point(self.oi[nob]),
                ) {
                    tpi.push(IntfSectionPoint::with_3d(
                        &first_pol.point(self.oi[nob]),
                        IntfPIType::Vertex,
                        self.oi[nob],
                        0,
                        0.0,
                        IntfPIType::Face,
                        tri2,
                        0,
                        0.0,
                        1.0,
                    ));
                    t_op[nob] = IntfPIType::Face;
                    nbp_int += 1;
                }
            }
        }

        // If the three points of <Tri1> are in <Tri2> the triangle Tri1 is
        // itself the tangent zone, else compute the positions of the points of
        // <Tri2> in <Tri1> (`cxx:939-1006`).
        if nbp_int < 3 {
            for nou in 0..=2usize {
                let _nou2 = POURCENT3[nou + 1];
                if t_tp[nou] == IntfPIType::External {
                    for nob in 0..=2usize {
                        let nob2 = POURCENT3[nob + 1];
                        if self.de_op_t[nob][nou].abs() <= tolerance {
                            if self.dp_op_t[nob][nou] + self.dp_op_t[nob2][nou]
                                < self.voo[nob].modulus()
                            {
                                let mut par = self.dp_op_t[nob][nou]
                                    / (self.dp_op_t[nob][nou] + self.dp_op_t[nob2][nou]);
                                if self.oi[nob] > self.oi[nob2] {
                                    par = 1.0 - par;
                                }
                                tpi.push(IntfSectionPoint::with_3d(
                                    &secon_pol.point(self.ti[nou]),
                                    IntfPIType::Edge,
                                    self.oi[nob].min(self.oi[nob2]),
                                    self.oi[nob].max(self.oi[nob2]),
                                    par,
                                    IntfPIType::Vertex,
                                    self.ti[nou],
                                    0,
                                    0.0,
                                    1.0,
                                ));
                                t_tp[nou] = IntfPIType::Edge;
                                nbp_int += 1;
                                break;
                            }
                        }
                    }
                    if t_tp[nou] == IntfPIType::External {
                        if contain(
                            &first_pol.point(self.oi[0]),
                            &first_pol.point(self.oi[1]),
                            &first_pol.point(self.oi[2]),
                            &secon_pol.point(self.ti[nou]),
                        ) {
                            tpi.push(IntfSectionPoint::with_3d(
                                &secon_pol.point(self.ti[nou]),
                                IntfPIType::Face,
                                tri1,
                                0,
                                0.0,
                                IntfPIType::Vertex,
                                self.ti[nou],
                                0,
                                0.0,
                                1.0,
                            ));
                            t_tp[nou] = IntfPIType::Face;
                            nbp_int += 1;
                        }
                    }
                }
            }
            if t_tp[0] != IntfPIType::External
                && t_tp[1] != IntfPIType::External
                && t_tp[2] != IntfPIType::External
            {
                finished = true;
            }
        } else {
            finished = true;
        }

        // Insertion of the points of intersection in the zone of tangency
        // (`cxx:1009-1016`).
        for nob in 1..=nbp_int {
            the_tz.append_point(&tpi[(nob - 1) as usize]);
        }

        if !finished {
            // If one of the triangles is not in the zone of tangency, find the
            // points of intersection edge/edge (`cxx:1019-1111`).
            let mut par_o = [0.0f64; 10];
            let mut par_t = [0.0f64; 10];
            let mut nb_no_inserted = 0i32;
            let mut pi_to_insert = [0i32; 17]; // for GCC 4.9

            for nob in 0..3usize {
                // processing of the object segment P[nob], P[nob+1]
                let nob2 = POURCENT3[nob + 1];

                for nou in 0..3usize {
                    // processing of the segment of the tool P[nou], P[nou+1]
                    let nou2 = POURCENT3[nou + 1];

                    if self.dp_oe_t[nob][nou] * self.dp_oe_t[nob2][nou] < 0.0
                        && self.de_op_t[nob][nou] * self.de_op_t[nob][nou2] < 0.0
                    {
                        if nbp_int >= 5 {
                            break;
                        } else {
                            let idx = nbp_int as usize;
                            par_o[idx] = self.dp_oe_t[nob][nou]
                                / (self.dp_oe_t[nob][nou] - self.dp_oe_t[nob2][nou]);
                            par_t[idx] = self.de_op_t[nob][nou]
                                / (self.de_op_t[nob][nou] - self.de_op_t[nob][nou2]);
                            let shift =
                                GpVec::from_xyz(&self.vtt[nou].multiplied(par_t[idx]));
                            let lepi = secon_pol
                                .point(self.ti[nou])
                                .translated_vec(&shift);
                            if self.oi[nob] > self.oi[nob2] {
                                par_o[idx] = 1.0 - par_o[idx];
                            }
                            if self.ti[nou] > self.ti[nou2] {
                                par_t[idx] = 1.0 - par_t[idx];
                            }
                            tpi.push(IntfSectionPoint::with_3d(
                                &lepi,
                                IntfPIType::Edge,
                                self.oi[nob].min(self.oi[nob2]),
                                self.oi[nob].max(self.oi[nob2]),
                                par_o[idx],
                                IntfPIType::Edge,
                                self.ti[nou].min(self.ti[nou2]),
                                self.ti[nou].max(self.ti[nou2]),
                                par_t[idx],
                                self.incidence,
                            ));
                            nbp_int += 1;
                            if !the_tz.insert(&tpi[(nbp_int - 1) as usize]) {
                                pi_to_insert[nb_no_inserted as usize] = nbp_int;
                                nb_no_inserted += 1;
                            }
                        }
                    }
                }
                if nbp_int >= 5 {
                    break; // Number of pi passed in TZ !
                }
            }

            // `cxx:1081-1111`: retry the refused points, newest first. The
            // ported `Intf_TangentZone::Insert` always answers `false`, so the
            // walk drains `nob` to -1 and the tail dump (no port) is skipped.
            let mut nob = nb_no_inserted - 1;
            while nob >= 0 {
                while !the_tz.insert(&tpi[(pi_to_insert[nob as usize] - 1) as usize]) {
                    nob -= 1;
                    if nob < 0 {
                        break;
                    }
                }
                if nob >= 0 {
                    nb_no_inserted -= 1;
                    while nob < nb_no_inserted {
                        pi_to_insert[nob as usize] =
                            pi_to_insert[(nob + 1) as usize];
                        nob += 1;
                    }
                    nob = nb_no_inserted - 1;
                }
            }
        }

        if nbp_int < 3 {
            nbp_int = 0;
        }
        nbp_int > 0
    }
}
