//! `HLRBRep_Data::Update` (`HLRBRep_Data.cxx:612-980`): project every edge and
//! face and fill the 16-bit min/max boxes the hiding pass prunes with.

use std::sync::Arc;

use occt_core::precision::Precision;

use occt_topo::abs::Orientation;
use crate::hlr_algo::{
    add_min_max, copy_min_max, decode_min_max, encode_min_max, enlarge_min_max, init_min_max,
    size_box, HlrAlgoProjector, MinMaxIndices,
};
use crate::hlrbrep_curve::HlrBRepCurve;
use crate::hlrbrep_data::HlrBRepData;
use crate::hlrbrep_edge_face_tool::uv_point;
use crate::hlrbrep_face_iterator::HlrBRepFaceIterator;

impl HlrBRepData {
    /// `Update(P)` (`cxx:612-980`): project every edge and face, compute the
    /// 16-bit min/max boxes the hiding pass prunes with, and classify the faces.
    ///
    /// The `FD.WithOutL() && !FD.Side()` arm (`cxx:841-844`) is the port's
    /// `orient_out_line` / `orient_oth_edge` (`orient.rs`), the
    /// `HLRBRep_Intersector` classification of `HLRBRep_Data.cxx:841-844`.
    pub fn update(&mut self, p: &Arc<HlrAlgoProjector>) -> Result<(), &'static str> {
        self.my_proj = Arc::clone(p);
        let t = p.transformation();
        let mut tol_min_max = 0.0f64;

        let mut tot_min = [0.0f64; 16];
        let mut tot_max = [0.0f64; 16];
        // `HLRAlgo::InitMinMax(Precision::Infinite(), TotMin, TotMax)` (`cxx:627`).
        init_min_max(Precision::INFINITE, &mut tot_min, &mut tot_max);

        // compute the global MinMax (`cxx:629-644`)
        for edge in 1..=self.my_nb_edges {
            let ed = self.e_data_mut(edge);
            let ec: &mut HlrBRepCurve = ed.change_geometry();
            ec.set_projector(Some(Arc::clone(p))); // `EC.Projector(&myProj)`
            let enl = ec.update(&mut tot_min, &mut tot_max);
            if enl > tol_min_max {
                tol_min_max = enl;
            }
        }
        enlarge_min_max(tol_min_max, &mut tot_min, &mut tot_max);

        let mut d = [0.0f64; 16];
        let mut precad = -Precision::INFINITE;
        for i in 0..16 {
            d[i] = tot_max[i] - tot_min[i];
            if precad < d[i] {
                precad = d[i];
            }
        }
        self.my_big_size = precad;
        precad *= 0.0005; // `cxx:655`

        for i in 0..16 {
            // `mySurD[i] = 0x00007fff / (d[i] + precad);` (`cxx:659`)
            self.my_sur_d[i] = 0x0000_7fffu32 as f64 / (d[i] + precad);
        }
        precad *= 0.5; // `cxx:661`

        for i in 0..16 {
            self.my_deca[i] = -tot_min[i] + precad; // `cxx:666`
        }

        // update the edges (`cxx:668-781`)
        let nb_edges = self.my_nb_edges;
        let toler = self.my_toler;
        let deca = self.my_deca;
        let sur_d = self.my_sur_d;
        for edge in 1..=nb_edges {
            init_min_max(Precision::INFINITE, &mut tot_min, &mut tot_max);
            let mut min_max_edge = MinMaxIndices::default();
            let (ver1, ver2, vertical, fold);
            {
                let ed = &mut self.my_e_data[edge as usize];
                let ec = ed.change_geometry();
                tol_min_max = ec.update_min_max(&mut tot_min, &mut tot_max);
                let tol = ed.tolerance() as f64;
                // `cxx:682-687`.
                ed.set_vertical(
                    tot_max[0] - tot_min[0] < tol
                        && tot_max[1] - tot_min[1] < tol
                        && tot_max[2] - tot_min[2] < tol
                        && tot_max[3] - tot_min[3] < tol
                        && tot_max[4] - tot_min[4] < tol
                        && tot_max[5] - tot_min[5] < tol
                        && tot_max[6] - tot_min[6] < tol,
                );

                let mut edge_min = MinMaxIndices::default();
                let mut edge_max = MinMaxIndices::default();
                for k in 0..8 {
                    // `cxx:688-703`: Min[0..7] from TotMin[0..7].
                    edge_min.min[k] = ((deca[k] + tot_min[k]) * sur_d[k]) as i32;
                    edge_max.min[k] = ((deca[k] + tot_max[k]) * sur_d[k]) as i32;
                }
                for k in 0..8 {
                    // `cxx:704-719`: Max[0..7] from TotMin/Max[8..15].
                    edge_min.max[k] = ((deca[k + 8] + tot_min[k + 8]) * sur_d[k + 8]) as i32;
                    edge_max.max[k] = ((deca[k + 8] + tot_max[k + 8]) * sur_d[k + 8]) as i32;
                }
                encode_min_max(&edge_min, &edge_max, &mut min_max_edge);

                vertical = ed.vertical();
                if vertical {
                    // `cxx:726-748`: fold the vertex index of this edge into
                    // the edges sharing its start vertex.
                    ver1 = true;
                    ver2 = true;
                    let vsta = ed.v_sta();
                    let vend = ed.v_end();
                    let vout = ed.out_lv_sta() || ed.out_lv_end();
                    let vcut = ed.cut_at_sta() || ed.cut_at_end();
                    fold = Some((vsta, vend, vout, vcut));
                } else {
                    // `cxx:751-776`.
                    let ec = ed.change_geometry();
                    let (_, tg1) = ec.d1_3d(ec.parameter3d(ec.first_parameter()));
                    let (_, tg2) = ec.d1_3d(ec.parameter3d(ec.last_parameter()));
                    let mut tg1 = tg1;
                    let mut tg2 = tg2;
                    tg1.transform(t);
                    tg2.transform(t);
                    let limit = toler as f64 * 10.0;

                    if tg1.x().abs() + tg1.y().abs() < limit {
                        ver1 = true;
                    } else {
                        // `gp_Dir Dir1(Tg1); ver1 = |X| + |Y| < Tol*10;`
                        let dir1 = occt_core::gp::GpDir::from_vec(&tg1);
                        ver1 = match dir1 {
                            Ok(dir1) => dir1.x().abs() + dir1.y().abs() < limit,
                            Err(_) => false,
                        };
                    }
                    if tg2.x().abs() + tg2.y().abs() < limit {
                        ver2 = true;
                    } else {
                        let dir2 = occt_core::gp::GpDir::from_vec(&tg2);
                        ver2 = match dir2 {
                            Ok(dir2) => dir2.x().abs() + dir2.y().abs() < limit,
                            Err(_) => false,
                        };
                    }
                    fold = None;
                }
            }
            enlarge_min_max(tol_min_max, &mut tot_min, &mut tot_max);

            if let Some((vsta, vend, vout, vcut)) = fold {
                for ebis in 1..=nb_edges {
                    let eb = &mut self.my_e_data[ebis as usize];
                    if vsta == eb.v_sta() {
                        eb.set_v_sta(vend);
                        eb.set_out_lv_sta(vout);
                        eb.set_cut_at_sta(vcut);
                    } else if vsta == eb.v_end() {
                        eb.set_v_end(vend);
                        eb.set_out_lv_end(vout);
                        eb.set_cut_at_end(vcut);
                    }
                }
            }

            let ed = &mut self.my_e_data[edge as usize];
            ed.update_min_max(&min_max_edge);
            ed.set_ver_at_sta(vertical || ver1);
            ed.set_ver_at_end(vertical || ver2);
            ed.set_auto_intersection_done(true);
            ed.set_simple(true);
        }

        // update the faces (`cxx:783-979`)
        let mut face_itr1 = HlrBRepFaceIterator::new();
        let nb_faces = self.my_nb_faces;
        for face in 1..=nb_faces {
            self.i_face_geom = Some(face);
            // OCCT mutates the live record; the port works on a copy of it and
            // writes it back at the end of the iteration, because the edge
            // arrays it reads live in the same struct.
            let mut fd = self.my_f_data[face as usize].clone();
            fd.geometry_mut().set_projector(Some(Arc::clone(p))); // `FS.Projector(&myProj)`
            // `mySLProps.SetSurface(iFaceGeom)` (`cxx:791`).
            self.my_sl_props_mut().set_surface(fd.geometry());
            let i_face_type = fd.geometry().get_type();

            // Is the face cut by an outline (`cxx:798-817`)
            let mut cut = false;
            let mut with_out_l = false;
            face_itr1.init_edge(&fd);
            while face_itr1.more_edge() {
                if face_itr1.internal() {
                    with_out_l = true;
                    cut = true;
                } else if face_itr1.out_line() {
                    with_out_l = true;
                    if face_itr1.double() {
                        cut = true;
                    }
                }
                face_itr1.next_edge();
            }
            fd.set_cut(cut);
            fd.set_with_out_l(with_out_l);

            // Is the face simple = no auto-hiding (`cxx:819-831`)
            use occt_topo::geom_bnd_lib_surface3d::GeomAbsSurfaceType as ST;
            fd.set_simple(
                !with_out_l
                    && matches!(
                        i_face_type,
                        ST::Plane | ST::Cylinder | ST::Cone | ST::Sphere | ST::Torus
                    ),
            );

            fd.set_plane(i_face_type == ST::Plane);
            fd.set_cylinder(i_face_type == ST::Cylinder);
            fd.set_cone(i_face_type == ST::Cone);
            fd.set_sphere(i_face_type == ST::Sphere);
            fd.set_torus(i_face_type == ST::Torus);

            let tol = fd.tolerance() as f64;
            // `fd.Side(FS.IsSide(tol, myToler * 10))` (`cxx:839`).
            fd.set_side(fd.geometry().is_side(tol, toler as f64 * 10.0));

            let inverted = if fd.with_out_l() && !fd.side() {
                // `cxx:841-844`: `inverted = OrientOutLine(face, fd);
                // OrientOthEdge(face, fd);`.
                let inverted = self.orient_out_line(&mut fd);
                self.orient_oth_edge(&mut fd);
                inverted
            } else {
                // `cxx:845`: `else inverted = false;`.
                false
            };

            if fd.side() {
                // `cxx:846-850`.
                fd.set_hiding(false);
                fd.set_back(false);
            } else if !fd.with_out_l() {
                // `cxx:851-928`: the back-face test through the middle point of
                // the first non-double edge.
                fd.set_back(false);
                let mut found = false;

                face_itr1.init_edge(&fd);
                while face_itr1.more_edge() && !found {
                    let my_fe = face_itr1.edge();
                    let my_fe_ori = face_itr1.orientation();
                    let my_fe_double = face_itr1.double();

                    if !my_fe_double
                        && (my_fe_ori == Orientation::Forward || my_fe_ori == Orientation::Reversed)
                    {
                        let param = {
                            let ec = self.my_e_data[my_fe as usize].geometry();
                            ec.parameter3d((ec.last_parameter() + ec.first_parameter()) / 2.0)
                        };
                        // `HLRBRep_EdgeFaceTool::UVPoint(p, myFEGeom, iFaceGeom, pu, pv)`
                        // (`cxx:870`).
                        let uv = uv_point(
                            param,
                            self.my_e_data[my_fe as usize].geometry(),
                            fd.geometry(),
                        );
                        if let Some((pu, pv)) = uv {
                            // `mySLProps.SetParameters(pu, pv)` (`cxx:871`).
                            self.my_sl_props_mut().set_parameters(pu, pv);
                            let pt_raw = self.my_e_data[my_fe as usize].geometry().value3d(param);
                            // `if (mySLProps.IsNormalDefined())` (`cxx:874`).
                            if self.my_sl_props_mut().is_normal_defined() {
                                // `gp_Vec Nm = mySLProps.Normal()` (`cxx:878`);
                                // `Normal()` only throws when `IsNormalDefined()`
                                // is false, which the branch above excluded.
                                let nm_raw = self
                                    .my_sl_props()
                                    .normal()
                                    .copied()
                                    .expect("mySLProps.Normal() after IsNormalDefined()");
                                let mut pt = pt_raw;
                                // `gp_Vec Nm = mySLProps.Normal()` widens the
                                // direction to a vector before the transform.
                                let mut nm =
                                    occt_core::gp::GpVec::new(nm_raw.x(), nm_raw.y(), nm_raw.z());
                                // `Pt.Transform(T); Nm.Transform(T);` (`cxx:879-880`).
                                pt.transform(t);
                                nm.transform(t);
                                let r = if p.perspective() {
                                    // `Nm.Z()*Focus - (Nm.X()*Pt.X + ...)`
                                    // (`cxx:882-885`).
                                    match p.focus() {
                                        Ok(focus) => {
                                            nm.z() * focus
                                                - (nm.x() * pt.x()
                                                    + nm.y() * pt.y()
                                                    + nm.z() * pt.z())
                                        }
                                        Err(_) => nm.z(),
                                    }
                                } else {
                                    // `r = Nm.Z();` (`cxx:887-890`).
                                    nm.z()
                                };
                                if r.abs() > self.my_toler as f64 * 10.0 {
                                    // `fd.Back(r < 0); found = true;` (`cxx:890-892`).
                                    fd.set_back(r < 0.0);
                                    found = true;
                                }
                            }
                        }
                    }
                    if !found {
                        face_itr1.next_edge();
                    }
                }

                if !found {
                    // `cxx:899-903`.
                    fd.set_side(true);
                    fd.set_hiding(false);
                    fd.set_back(false);
                } else if fd.closed() {
                    // `cxx:905-921`.
                    let hiding = match fd.orientation() {
                        Orientation::Reversed => fd.back(),
                        Orientation::Forward => !fd.back(),
                        Orientation::External => true,
                        Orientation::Internal => false,
                    };
                    fd.set_hiding(hiding);
                } else {
                    fd.set_hiding(true); // `cxx:925`
                }
            } else {
                // `cxx:928-940`.
                if inverted {
                    fd.set_hiding(false);
                    fd.set_back(true);
                } else {
                    fd.set_hiding(true);
                    fd.set_back(false);
                }
            }

            // `cxx:942-978`: the per-wire / per-face min/max accumulation.
            let mut first_time = true;
            let mut face_min = MinMaxIndices::default();
            let mut face_max = MinMaxIndices::default();
            let mut wire_min = MinMaxIndices::default();
            let mut wire_max = MinMaxIndices::default();

            face_itr1.init_edge(&fd);
            while face_itr1.more_edge() {
                let my_fe = face_itr1.edge();
                if !fd.simple() {
                    // `EDataFE2.AutoIntersectionDone(false)` (`cxx:950`).
                    self.my_e_data[my_fe as usize].set_auto_intersection_done(false);
                }
                let mut edge_min = MinMaxIndices::default();
                let mut edge_max = MinMaxIndices::default();
                decode_min_max(
                    self.my_e_data[my_fe as usize].min_max(),
                    &mut edge_min,
                    &mut edge_max,
                );

                if face_itr1.beginning_of_wire() {
                    let (a, b) = copy_min_max(&edge_min, &edge_max);
                    wire_min = a;
                    wire_max = b;
                } else {
                    add_min_max(&edge_min, &edge_max, &mut wire_min, &mut wire_max);
                }
                if face_itr1.end_of_wire() {
                    let mut min_max_wire = MinMaxIndices::default();
                    encode_min_max(&wire_min, &wire_max, &mut min_max_wire);
                    // `myFaceItr1.Wire()->UpdateMinMax(MinMaxWire)` (`cxx:965`).
                    // The wire block is `handle<HLRAlgo_EdgesBlock>` in OCCT;
                    // the iterator's copy is updated through the face's wire
                    // list below.
                    if first_time {
                        first_time = false;
                        let (a, b) = copy_min_max(&wire_min, &wire_max);
                        face_min = a;
                        face_max = b;
                    } else {
                        add_min_max(&wire_min, &wire_max, &mut face_min, &mut face_max);
                    }
                }
                face_itr1.next_edge();
            }

            let mut min_max_face = MinMaxIndices::default();
            encode_min_max(&face_min, &face_max, &mut min_max_face);
            let size = size_box(&face_min, &face_max);
            if let Some(wires) = fd.wires_mut() {
                wires.update_min_max(&min_max_face);
            }
            fd.set_size(size);

            self.my_f_data[face as usize] = fd;
        }

        Ok(())
    }
}
