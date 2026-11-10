//! `IntPatch_InterferencePolyhedron::Intersect`.
//!
//! Source: `IntPatch_InterferencePolyhedron.cxx:151-838`.
//!
//! The body is one long cascade of singularities, in the exact order OCCT
//! tests them, each with its own bookkeeping of which edges of the two
//! triangles have been consumed:
//!
//! 1. `VERTEX VERTEX`      `cxx:288-311`
//! 2. `VERTEX EDGE`        `cxx:313-354`
//! 3. `EDGE VERTEX`        `cxx:357-394`
//! 4. `FACE VERTEX`        `cxx:397-421`
//! 5. `VERTEX FACE`        `cxx:423-446`
//! 6. `EDGE EDGE`          `cxx:448-551`
//! 7. `EDGE FACE`          `cxx:553-584`
//! 8. `FACE EDGE`          `cxx:586-615`
//! 9. sort by the section direction `ONor ^ TNor`   `cxx:619-655`
//! 10. the section-line insertion rules             `cxx:658-837`
//!
//! All `Intf_SectionPoint` addresses are the OCCT 1-based point indices
//! (`OI` / `TI`), so `ParamOnFirst` / `ParamOnSecond` keep their meaning.
//! `piOT` is kept as a plain `Vec` with explicit `index - 1` conversions where
//! the OCCT code indexes `piOT(i)` directly.

use occt_core::gp::{GpPnt, GpXyz};
use occt_core::intf::{plane_equation, IntfPIType, IntfSectionPoint, IntfTangentZone};

use super::{IntPatchInterferencePolyhedron, POURCENT3};
use crate::intpatch_polyhedron::IntPatchPolyhedron;

impl IntPatchInterferencePolyhedron {
    /// `Intersect(Tri1, FirstPol, Tri2, SeconPol)` (`cxx:151-838`).
    pub(crate) fn intersect(
        &mut self,
        tri1: i32,
        first_pol: &IntPatchPolyhedron,
        tri2: i32,
        secon_pol: &IntPatchPolyhedron,
    ) {
        let (o0, o1, o2) = first_pol.triangle(tri1 as usize);
        self.oi = [o0, o1, o2];
        let (t0, t1, t2) = secon_pol.triangle(tri2 as usize);
        self.ti = [t0, t1, t2];

        // If there is an intersection of a polyhedron with itself, the
        // intersections are excluded from a triangle with connected triangles
        // (`cxx:161-172`).
        if self.base.self_intf {
            for i in 0..3 {
                for j in 0..3 {
                    if self.oi[i] == self.ti[j] {
                        return;
                    }
                }
            }
        }

        // `double floatGap = 1e-13; //-- Epsilon(1000.);` (`cxx:186`).
        let float_gap = 1e-13;

        // Equation of the triangle plane of the objet (`cxx:188-195`).
        let (onor, odp) = plane_equation(
            &first_pol.point(self.oi[0]),
            &first_pol.point(self.oi[1]),
            &first_pol.point(self.oi[2]),
        );

        // Equation of the triangle plane of the tool (`cxx:197-204`).
        let (tnor, tdp) = plane_equation(
            &secon_pol.point(self.ti[0]),
            &secon_pol.point(self.ti[1]),
            &secon_pol.point(self.ti[2]),
        );

        // Scalar product of two normalized vectors -> cosinus of the angle.
        self.incidence = tnor.dot(&onor).abs();

        // Distance of the plane of the triangle from the object by three points
        // of SeconPol (`cxx:209-212`).
        let mut df_op_t = [0.0f64; 3];
        for i in 0..3 {
            df_op_t[i] = onor.dot(secon_pol.point(self.ti[i]).xyz()) - odp;
        }

        // Distance of the plane of the triangle from the tool by three points
        // of FirstPol (`cxx:215-218`).
        let mut dp_of_t = [0.0f64; 3];
        for i in 0..3 {
            dp_of_t[i] = tnor.dot(first_pol.point(self.oi[i]).xyz()) - tdp;
        }

        // Values defining the couple of triangle dpOpT, dpOeT, deOpT
        // (`cxx:221`).
        self.couple_characteristics(first_pol, secon_pol);

        let tolerance = self.base.tolerance;

        // If three points of the triangle of <SeconPol> are in the plane of the
        // triangle of <Obje> within <Tolerance> the eventual tangency zone is
        // found (`cxx:227-236`).
        let mut is_tangent = true;
        for i in 0..3 {
            if !(df_op_t[i].abs() <= tolerance) || !(dp_of_t[i].abs() <= tolerance) {
                is_tangent = false;
            }
        }
        if is_tangent {
            let df_sum = df_op_t[0] + df_op_t[1] + df_op_t[2];
            if df_sum.abs() == df_op_t[0].abs() + df_op_t[1].abs() + df_op_t[2].abs() {
                is_tangent = false;
            }
            let dp_sum = dp_of_t[0] + dp_of_t[1] + dp_of_t[2];
            if dp_sum.abs() == dp_of_t[0].abs() + dp_of_t[1].abs() + dp_of_t[2].abs() {
                is_tangent = false;
            }
        }

        if is_tangent {
            let mut the_tz = IntfTangentZone::new();
            if self.tangent_zone_value(&mut the_tz, first_pol, tri1, secon_pol, tri2) {
                self.store_tangent_zone(&the_tz);
            }
            return;
        }

        // Otherwise line of section is calculated (`cxx:246-...`).
        let mut nbp_ot = 0i32;
        let mut nbp_o = 0i32;
        let mut nbp_t = 0i32;
        let mut pi_ot: Vec<IntfSectionPoint> = Vec::new();

        let mut par_o = [-1.0f64; 3];
        let mut par_t = [-1.0f64; 3];
        let mut ed_ot = [1i32; 3];
        let mut ed_tt = [1i32; 3];

        let incidence = self.incidence;

        // Singularite VERTEX VERTEX (`cxx:288-311`).
        for i_obj in 0..3 {
            for i_too in 0..3 {
                if self.dp_op_t[i_obj][i_too] <= float_gap {
                    pi_ot.push(IntfSectionPoint::with_3d(
                        &first_pol.point(self.oi[i_obj]),
                        IntfPIType::Vertex,
                        self.oi[i_obj],
                        0,
                        0.0,
                        IntfPIType::Vertex,
                        self.ti[i_too],
                        0,
                        0.0,
                        incidence,
                    ));
                    par_o[i_obj] = 0.0;
                    par_t[i_too] = 0.0;
                    ed_ot[POURCENT3[i_obj + 2]] = 0;
                    ed_ot[i_obj] = 0;
                    ed_tt[POURCENT3[i_too + 2]] = 0;
                    ed_tt[i_too] = 0;
                    nbp_ot += 1;
                    nbp_o += 1;
                    nbp_t += 1;
                }
            }
        }

        // Singularite VERTEX EDGE (`cxx:313-354`).
        for i_obj in 0..3 {
            if par_o[i_obj] == -1.0 {
                for i_too in 0..3 {
                    let inext = POURCENT3[i_too + 1];
                    if ed_tt[i_too] == 1 {
                        if self.dp_oe_t[i_obj][i_too] <= float_gap
                            && self.dp_oe_t[i_obj][i_too] >= -float_gap
                        {
                            if self.dp_op_t[i_obj][i_too] + self.dp_op_t[i_obj][inext]
                                < self.vtt[i_too].modulus()
                            {
                                par_t[i_too] = self.dp_op_t[i_obj][i_too]
                                    / (self.dp_op_t[i_obj][i_too] + self.dp_op_t[i_obj][inext]);
                                if self.ti[i_too] > self.ti[inext] {
                                    par_t[i_too] = 1.0 - par_t[i_too];
                                }
                                pi_ot.push(IntfSectionPoint::with_3d(
                                    &first_pol.point(self.oi[i_obj]),
                                    IntfPIType::Vertex,
                                    self.oi[i_obj],
                                    0,
                                    0.0,
                                    IntfPIType::Edge,
                                    self.ti[i_too].min(self.ti[inext]),
                                    self.ti[i_too].max(self.ti[inext]),
                                    par_t[i_too],
                                    incidence,
                                ));
                                par_o[i_obj] = 0.0;
                                ed_ot[POURCENT3[i_obj + 2]] = 0;
                                ed_ot[i_obj] = 0;
                                ed_tt[i_too] = 0;
                                nbp_ot += 1;
                                nbp_o += 1;
                                nbp_t += 1;
                            }
                        }
                    }
                }
            }
        }

        // Singularite EDGE VERTEX (`cxx:357-394`).
        for i_too in 0..3 {
            if par_t[i_too] == -1.0 {
                for i_obj in 0..3 {
                    let inext = POURCENT3[i_obj + 1];
                    if ed_ot[i_obj] == 1 {
                        if self.de_op_t[i_obj][i_too].abs() <= float_gap {
                            if self.dp_op_t[i_obj][i_too] + self.dp_op_t[inext][i_too]
                                < self.voo[i_obj].modulus()
                            {
                                par_o[i_obj] = self.dp_op_t[i_obj][i_too]
                                    / (self.dp_op_t[i_obj][i_too] + self.dp_op_t[inext][i_too]);
                                if self.oi[i_obj] > self.oi[inext] {
                                    par_o[i_obj] = 1.0 - par_o[i_obj];
                                }
                                pi_ot.push(IntfSectionPoint::with_3d(
                                    &secon_pol.point(self.ti[i_too]),
                                    IntfPIType::Edge,
                                    self.oi[i_obj].min(self.oi[inext]),
                                    self.oi[i_obj].max(self.oi[inext]),
                                    par_o[i_obj],
                                    IntfPIType::Vertex,
                                    self.ti[i_too],
                                    0,
                                    0.0,
                                    incidence,
                                ));
                                par_t[i_too] = 0.0;
                                ed_ot[i_obj] = 0;
                                ed_tt[POURCENT3[i_too + 2]] = 0;
                                ed_tt[i_too] = 0;
                                nbp_ot += 1;
                                nbp_o += 1;
                                nbp_t += 1;
                            }
                        }
                    }
                }
            }
        }

        // Singularite FACE VERTEX (`cxx:397-421`).
        for i_too in 0..3 {
            if par_t[i_too] != 0.0 {
                if df_op_t[i_too].abs() <= float_gap {
                    pi_ot.push(IntfSectionPoint::with_3d(
                        &secon_pol.point(self.ti[i_too]),
                        IntfPIType::Face,
                        tri1,
                        0,
                        0.0,
                        IntfPIType::Vertex,
                        self.ti[i_too],
                        0,
                        0.0,
                        incidence,
                    ));
                    par_t[i_too] = 0.0;
                    ed_tt[POURCENT3[i_too + 2]] = 0;
                    ed_tt[i_too] = 0;
                    nbp_ot += 1;
                    nbp_t += 1;
                }
            }
        }

        // Singularite VERTEX FACE (`cxx:423-446`).
        for i_obj in 0..3 {
            if par_o[i_obj] != 0.0 {
                if dp_of_t[i_obj].abs() <= float_gap {
                    pi_ot.push(IntfSectionPoint::with_3d(
                        &first_pol.point(self.oi[i_obj]),
                        IntfPIType::Vertex,
                        self.oi[i_obj],
                        0,
                        0.0,
                        IntfPIType::Face,
                        tri2,
                        0,
                        0.0,
                        incidence,
                    ));
                    par_o[i_obj] = 0.0;
                    ed_ot[POURCENT3[i_obj + 2]] = 0;
                    ed_ot[i_obj] = 0;
                    nbp_ot += 1;
                    nbp_o += 1;
                }
            }
        }

        // Singularite EDGE EDGE (`cxx:448-551`).
        for i_obj in 0..3 {
            let inext = POURCENT3[i_obj + 1];
            if ed_ot[i_obj] == 1 && (dp_of_t[i_obj] * dp_of_t[inext]) < 0.0 {
                let mut lg = dp_of_t[i_obj] / (dp_of_t[i_obj] - dp_of_t[inext]);
                if lg > 0.0 && lg < 1.0 {
                    for i_too in 0..3 {
                        let jnext = POURCENT3[i_too + 1];
                        if ed_tt[i_too] == 1 && (df_op_t[i_too] * df_op_t[jnext]) < 0.0 {
                            lg = df_op_t[i_too] / (df_op_t[i_too] - df_op_t[jnext]);
                            if lg > 0.0 && lg < 1.0 {
                                let mut pb = false;
                                let mut pi_o = GpXyz::new(0.0, 0.0, 0.0);
                                let mut pi_t = GpXyz::new(0.0, 0.0, 0.0);

                                if self.oi[i_obj] > self.oi[inext] {
                                    let div =
                                        self.dp_oe_t[inext][i_too] - self.dp_oe_t[i_obj][i_too];
                                    if div > float_gap || div < -float_gap {
                                        par_o[i_obj] = self.dp_oe_t[inext][i_too] / div;
                                        pi_o = first_pol
                                            .point(self.oi[inext])
                                            .xyz()
                                            .added(&self.voo[i_obj].reversed().multiplied(par_o[i_obj]));
                                    } else {
                                        pb = true;
                                    }
                                } else {
                                    let div =
                                        self.dp_oe_t[i_obj][i_too] - self.dp_oe_t[inext][i_too];
                                    if div > float_gap || div < -float_gap {
                                        par_o[i_obj] = self.dp_oe_t[i_obj][i_too]
                                            / (self.dp_oe_t[i_obj][i_too]
                                                - self.dp_oe_t[inext][i_too]);
                                        pi_o = first_pol
                                            .point(self.oi[i_obj])
                                            .xyz()
                                            .added(&self.voo[i_obj].multiplied(par_o[i_obj]));
                                    } else {
                                        pb = true;
                                    }
                                }

                                if self.ti[i_too] > self.ti[jnext] {
                                    let div =
                                        self.de_op_t[i_obj][jnext] - self.de_op_t[i_obj][i_too];
                                    if div > float_gap || div < -float_gap {
                                        par_t[i_too] = self.de_op_t[i_obj][jnext] / div;
                                        pi_t = secon_pol
                                            .point(self.ti[jnext])
                                            .xyz()
                                            .added(&self.vtt[i_too].reversed().multiplied(par_t[i_too]));
                                    } else {
                                        pb = true;
                                    }
                                } else {
                                    let div =
                                        self.de_op_t[i_obj][i_too] - self.de_op_t[i_obj][jnext];
                                    if div > float_gap || div < -float_gap {
                                        par_t[i_too] = self.de_op_t[i_obj][i_too] / div;
                                        pi_t = secon_pol
                                            .point(self.ti[i_too])
                                            .xyz()
                                            .added(&self.vtt[i_too].multiplied(par_t[i_too]));
                                    } else {
                                        pb = true;
                                    }
                                }

                                if !pb {
                                    // `piT -= piO.XYZ(); lg = piT.Modulus();`
                                    let diff = pi_t.subtracted(&pi_o);
                                    lg = diff.modulus();
                                    if lg <= float_gap {
                                        pi_ot.push(IntfSectionPoint::with_3d(
                                            &GpPnt::from_xyz(&pi_o),
                                            IntfPIType::Edge,
                                            self.oi[i_obj].min(self.oi[inext]),
                                            self.oi[i_obj].max(self.oi[inext]),
                                            par_o[i_obj],
                                            IntfPIType::Edge,
                                            self.ti[i_too].min(self.ti[jnext]),
                                            self.ti[i_too].max(self.ti[jnext]),
                                            par_t[i_too],
                                            incidence,
                                        ));
                                        ed_ot[i_obj] = 0;
                                        ed_tt[i_too] = 0;
                                        nbp_ot += 1;
                                        nbp_o += 1;
                                        nbp_t += 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Intersection EDGE FACE (`cxx:553-584`).
        for i_obj in 0..3 {
            let inext = POURCENT3[i_obj + 1];
            if ed_ot[i_obj] == 1 && (dp_of_t[i_obj] * dp_of_t[inext]) < 0.0 {
                let lg = dp_of_t[i_obj] / (dp_of_t[i_obj] - dp_of_t[inext]);
                if lg > 0.0 && lg < 1.0 {
                    par_o[i_obj] = lg;
                    let pi_o = first_pol
                        .point(self.oi[i_obj])
                        .xyz()
                        .added(&self.voo[i_obj].multiplied(par_o[i_obj]));
                    if self.oi[i_obj] > self.oi[inext] {
                        par_o[i_obj] = 1.0 - par_o[i_obj];
                    }
                    pi_ot.push(IntfSectionPoint::with_3d(
                        &GpPnt::from_xyz(&pi_o),
                        IntfPIType::Edge,
                        self.oi[i_obj].min(self.oi[inext]),
                        self.oi[i_obj].max(self.oi[inext]),
                        par_o[i_obj],
                        IntfPIType::Face,
                        tri2,
                        0,
                        0.0,
                        incidence,
                    ));
                    nbp_ot += 1;
                    nbp_o += 1;
                }
            }
        }

        // Intersection FACE EDGE (`cxx:586-615`).
        for i_too in 0..3 {
            let jnext = POURCENT3[i_too + 1];
            if ed_tt[i_too] == 1 && (df_op_t[i_too] * df_op_t[jnext]) < 0.0 {
                let lg = df_op_t[i_too] / (df_op_t[i_too] - df_op_t[jnext]);
                if lg > 0.0 && lg < 1.0 {
                    par_t[i_too] = lg;
                    let pi_o = secon_pol
                        .point(self.ti[i_too])
                        .xyz()
                        .added(&self.vtt[i_too].multiplied(par_t[i_too]));
                    if self.ti[i_too] > self.ti[jnext] {
                        par_t[i_too] = 1.0 - par_t[i_too];
                    }
                    pi_ot.push(IntfSectionPoint::with_3d(
                        &GpPnt::from_xyz(&pi_o),
                        IntfPIType::Face,
                        tri1,
                        0,
                        0.0,
                        IntfPIType::Edge,
                        self.ti[i_too].min(self.ti[jnext]),
                        self.ti[i_too].max(self.ti[jnext]),
                        par_t[i_too],
                        incidence,
                    ));
                    nbp_ot += 1;
                    nbp_t += 1;
                }
            }
        }

        let mut id: Vec<usize> = vec![0; (nbp_ot.max(4)) as usize];
        let mut ideb: i32 = -1;
        let mut ifin: i32 = -2;

        // Sort the <nbpiOT> sections points along the intersection between the
        // two triangles (`cxx:625-655`).
        if nbp_ot > 1 {
            let dir = onor.crossed(&tnor);
            let mut d = vec![0.0f64; nbp_ot as usize];
            for i_pi in 0..nbp_ot as usize {
                d[i_pi] = dir.dot(pi_ot[i_pi].pnt().xyz());
            }

            id[0] = 0;
            for i_pi in 1..nbp_ot as usize {
                id[i_pi] = i_pi;
                let mut i_ps = i_pi as i32 - 1;
                while i_ps >= 0 {
                    let left = id[i_ps as usize];
                    let right = id[(i_ps + 1) as usize];
                    if d[left] > d[right] {
                        id.swap(i_ps as usize, (i_ps + 1) as usize);
                        i_ps -= 1;
                    } else {
                        break;
                    }
                }
            }
        }

        // Possibility of line of section (`cxx:658-837`).
        if nbp_o == 2 && nbp_t == 2 {
            // In the case when an edge is in the plane of the other triangle it
            // is necessary to check if it has not been already processed on a
            // connected triangle (`cxx:663-698`).
            let mut pivo: i32 = -1;
            let mut pedg: i32 = -1;
            if par_o[0] == 0.0 {
                pivo = 0;
                if par_o[1] == 0.0 {
                    pedg = 1;
                } else if par_o[2] == 0.0 {
                    pedg = 2;
                }
            } else if par_o[1] == 0.0 {
                pivo = 1;
                if par_o[2] == 0.0 {
                    pedg = 2;
                }
            }
            if pivo >= 0 && pedg >= 0 {
                // `pivo` and `pedg` are the `theTriCon` / `theOtherP` out
                // parameters of `TriConnex` (`IntPatch_PolyhedronTool.lxx:54-61`).
                let (_, tri_con, _) = first_pol.tri_connex(
                    tri1,
                    self.oi[pivo as usize],
                    self.oi[pedg as usize],
                );
                pivo = tri_con;
                // cxx:718 -- TriConnex also writes `theOtherP`, but OCCT never
                // reads it: `pedg` is reset at cxx:729 before any use.
                if pivo > tri1 {
                    nbp_ot = 0;
                    ideb = -1; // On a deja trouve celle ci
                    ifin = -2;
                }
            }

            // For the tool (`cxx:701-734`).
            pivo = -1;
            pedg = -1;
            if par_t[0] == 0.0 {
                pivo = 0;
                if par_t[1] == 0.0 {
                    pedg = 1;
                } else if par_t[2] == 0.0 {
                    pedg = 2;
                }
            } else if par_t[1] == 0.0 {
                pivo = 1;
                if par_t[2] == 0.0 {
                    pedg = 2;
                }
            }
            if pivo >= 0 && pedg >= 0 {
                let (_, tri_con, _) = secon_pol.tri_connex(
                    tri2,
                    self.ti[pivo as usize],
                    self.ti[pedg as usize],
                );
                pivo = tri_con;
                // cxx:746 -- same as cxx:718: the `theOtherP` out parameter is
                // discarded by OCCT as well.
                if pivo > tri2 {
                    nbp_ot = 0;
                    ideb = -1; // It has been already found
                    ifin = -2;
                }
            }

            if nbp_ot > 0 {
                // If there is a covering up: insert the section line in the
                // existent list or create a new section line (`cxx:737-793`).
                let first_type = pi_ot[id[0]].type_on_first();
                let second_type = pi_ot[id[0]].type_on_second();
                if first_type == IntfPIType::Face {
                    if pi_ot[id[1]].type_on_first() == IntfPIType::Face {
                        ideb = -(id[0] as i32) - 1; // No line of section possible
                        ifin = -(id[1] as i32) - 1;
                    } else if pi_ot[id[1]].type_on_second() != IntfPIType::Face {
                        ideb = id[1] as i32; // No line of section possible
                        ifin = id[1] as i32; // only a pointersec
                    } else if nbp_ot >= 3 {
                        ideb = id[1] as i32; // Retrieve 2 segments of section
                        ifin = id[2] as i32;
                    } else {
                        ideb = -999; // No line of section possible
                        ifin = -999;
                    }
                } else if second_type == IntfPIType::Face {
                    if pi_ot[id[1]].type_on_second() == IntfPIType::Face {
                        ideb = -(id[0] as i32) - 1; // No line of section possible
                        ifin = -(id[1] as i32) - 1;
                    } else if pi_ot[id[1]].type_on_first() != IntfPIType::Face {
                        ideb = id[1] as i32; // No line of section possible
                        ifin = id[1] as i32; // only a pointersec
                    } else if nbp_ot >= 3 {
                        ideb = id[1] as i32; // Recouvrement des 2 segments
                        ifin = id[2] as i32;
                    } else {
                        ideb = -999; // No line of section possible
                        ifin = -999;
                    }
                } else {
                    // Singularity on the first point: there is only two or
                    // three pointersec, so the first is a solution and the
                    // second too.
                    ideb = id[0] as i32;
                    ifin = id[1] as i32;
                }
            }

            // Insertion of the segment found in the existing section lines
            // (`cxx:795-837`).
            if ideb < 0 {
                if ifin < 0 {
                    if ideb != -999 {
                        let p = pi_ot[(-ideb) as usize - 1];
                        let q = pi_ot[(-ifin) as usize - 1];
                        let d = p.pnt().distance(q.pnt());
                        if d < self.base.tolerance {
                            self.base.insert_segment(&p, &q);
                        }
                    }
                }
            } else if ideb >= 0 {
                if ideb != ifin {
                    // `piOT(ideb + 1)` / `piOT(ifin + 1)`: 1-based OCCT index.
                    let p = pi_ot[ideb as usize];
                    let q = pi_ot[ifin as usize];
                    self.base.insert_segment(&p, &q);
                }
            }
        }
    }
}
