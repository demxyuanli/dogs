//! `MAT2d_Circuit::DoubleLine` (`MAT2d_Circuit.cxx:507-670`).
//!
//! An open line of the figure is doubled so that the circuit is closed: the
//! items after the first one are appended in reverse, the connexions of the
//! mini path are moved onto the copy whose side matches `SideRef`, and the
//! re-entrant corners that the doubling made non salient are dropped from the
//! line. `Corres` keeps the source index of every remaining item so that the
//! connexions can be renumbered.

use super::corners::cross_prod;
use super::Mat2dCircuit;
use crate::mat2d::item::Mat2dItem;
use crate::mat2d::mini_path::ConnexionHandle;
use occt_geom2d::curve::Curve2d;

impl Mat2dCircuit {
    /// `DoubleLine(Line, ConnexionFrom, ConnexionFather, SideRef)`
    /// (`cxx:507-670`).
    pub(crate) fn double_line(
        &self,
        line: &mut Vec<Mat2dItem>,
        connexion_from: &mut Vec<ConnexionHandle>,
        connexion_father: Option<ConnexionHandle>,
        side_ref: f64,
    ) {
        let nb_items = line.len() as i32;

        // `cxx:519-538`: completion of the line. A point is shared with the
        // original item (`Line.Append(Line.Value(i))`), a curve is copied and
        // reversed.
        if !self.is_open_result() {
            let mut i = nb_items - 1;
            while i > 1 {
                let item = match line[(i - 1) as usize] {
                    Mat2dItem::Point(p) => Mat2dItem::Point(p),
                    Mat2dItem::Curve(ref c) => {
                        let mut reversed = c.clone();
                        reversed.reverse();
                        Mat2dItem::Curve(reversed)
                    }
                };
                line.push(item);
                i -= 1;
            }
        }

        // `cxx:540-578`: repartition of the connexions over the line.
        let mut i_after = connexion_from.len() as i32;

        let mut i = 1i32;
        while i <= i_after {
            let cc = connexion_from[(i - 1) as usize].clone();
            let ind_cof = cc.borrow().index_item_on_first();

            if line[(ind_cof - 1) as usize].is_point() {
                if ind_cof != nb_items && ind_cof != 1 {
                    let mut dot_prod = 0.0;
                    let pro_vec = cross_prod(
                        &line[(ind_cof - 2) as usize],
                        &line[ind_cof as usize],
                        &mut dot_prod,
                    );
                    if pro_vec * side_ref > 0.0 {
                        cc.borrow_mut().set_index_item_on_first(2 * nb_items - ind_cof);
                        connexion_from.insert(i_after as usize, cc);
                        connexion_from.remove((i - 1) as usize);
                        i_after -= 1;
                        i -= 1;
                    }
                }
            } else {
                let cc_side = {
                    let c = cc.borrow();
                    self.side(&c, line)
                };
                if cc_side != side_ref {
                    let parameter = {
                        let c = cc.borrow();
                        line[(ind_cof - 1) as usize]
                            .reversed_parameter(c.parameter_on_first())
                            .expect("MAT2d_Circuit::DoubleLine: item is not a curve")
                    };
                    cc.borrow_mut().set_index_item_on_first(2 * nb_items - ind_cof);
                    cc.borrow_mut().set_parameter_on_first(parameter);
                    connexion_from.insert(i_after as usize, cc);
                    connexion_from.remove((i - 1) as usize);
                    i_after -= 1;
                    i -= 1;
                }
            }
            i += 1;
        }

        // `cxx:580-607`: the father connexion is moved the same way, on its
        // second side.
        if let Some(ref father) = connexion_father {
            let cc = father.borrow().reversed();
            let ind_cof = cc.index_item_on_first();

            if line[(ind_cof - 1) as usize].is_point() {
                if ind_cof != nb_items && ind_cof != 1 {
                    let mut dot_prod = 0.0;
                    let pro_vec = cross_prod(
                        &line[(ind_cof - 2) as usize],
                        &line[ind_cof as usize],
                        &mut dot_prod,
                    );
                    if pro_vec * side_ref > 0.0 {
                        father
                            .borrow_mut()
                            .set_index_item_on_second(2 * nb_items - ind_cof);
                    }
                }
            } else {
                let cc_side = self.side(&cc, line);
                if cc_side != side_ref {
                    let parameter = {
                        let p = father.borrow().parameter_on_second();
                        line[(ind_cof - 1) as usize]
                            .reversed_parameter(p)
                            .expect("MAT2d_Circuit::DoubleLine: item is not a curve")
                    };
                    father
                        .borrow_mut()
                        .set_index_item_on_second(2 * nb_items - ind_cof);
                    father.borrow_mut().set_parameter_on_second(parameter);
                }
            }
        }

        // `cxx:609-634`: removal of the re-entrant corners. OCCT stops the scan
        // on the last item with `Line.Value(IndLine) != Line.Last()` (`cxx:616`),
        // a handle identity test; the scan covers the whole doubled line, which
        // is what the position test below does.
        let mut ind_line = 1i32;
        let mut i_corres = 0usize;
        let mut corres = vec![0i32; line.len()];

        while ind_line != line.len() as i32 {
            corres[i_corres] = ind_line;
            let is_corner = line[(ind_line - 1) as usize].is_point()
                && i_corres != 0
                && i_corres as i32 + 1 != nb_items;

            if is_corner
                && !self.is_sharp_corner(
                    &line[(ind_line - 2) as usize],
                    &line[ind_line as usize],
                    side_ref,
                )
            {
                line.remove((ind_line - 1) as usize);
                ind_line -= 1;
                corres[i_corres] = 0;
            }
            ind_line += 1;
            i_corres += 1;
        }
        corres[i_corres] = ind_line;

        if !self.is_open_result() {
            // `cxx:638-644`: the items removed from one half are removed from
            // its mirror image as well.
            let last = 2 * nb_items - 3;
            for i in 1..=last {
                if corres[(i - 1) as usize] == 0 {
                    let j = (2 * nb_items - i - 1) as usize;
                    if j < corres.len() {
                        corres[(i - 1) as usize] = corres[j];
                    }
                }
            }

            // `cxx:659-663`: the connexions are renumbered. A corner that was
            // dropped has no image; OCCT indexes `Corres(0)` there, which is out
            // of the array, so the port keeps the index inside the line.
            for cc in connexion_from.iter() {
                let old = cc.borrow().index_item_on_first();
                let new = corres
                    .get((old - 1) as usize)
                    .copied()
                    .unwrap_or(0)
                    .clamp(1, line.len() as i32);
                cc.borrow_mut().set_index_item_on_first(new);
            }

            // `cxx:665-668`.
            if let Some(ref father) = connexion_father {
                let old = father.borrow().index_item_on_second();
                let new = corres
                    .get((old - 1) as usize)
                    .copied()
                    .unwrap_or(0)
                    .clamp(1, line.len() as i32);
                father.borrow_mut().set_index_item_on_second(new);
            }
        }
    }
}
