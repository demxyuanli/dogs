//! `FDS_reduceONFACEinterferences` — the `(ON, FACE)` filter of
//! `TopOpeBRepDS_EIR`.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EIR.cxx:49-105` (OCCT 8.0.0).
//!
//! When the list holds a point interference whose transition is `(ON, FACE)`,
//! every other point interference of the list is dropped (`cxx:75-103`).
//! `BDS` and `SIX` are unused by the C++ body (both parameters are unnamed in
//! the definition, `cxx:51-52`).

use crate::abs::{ShapeType, State};

use super::super::interference_iterator::InterferenceList;
use super::super::kind::Kind;
use super::super::process_interferences_tool::fun_has_state_shape;

/// `FDS_reduceONFACEinterferences(LI, BDS, SIX)` (`cxx:49-105`).
pub fn fds_reduce_onface_interferences(li: &mut InterferenceList) {
    // set hasONFACE = True if LI contains interfs with (ON,FACE) transition(s).
    let mut has_onface = false;
    for h in li.iter() {
        let i = h.read().expect("interference lock");
        if i.geometry_type() == Kind::Point {
            has_onface = fun_has_state_shape(&i.transition(), State::On, ShapeType::Face);
            if has_onface {
                break;
            }
        }
    }

    if !has_onface {
        return;
    }

    // LI has (ON,FACE) : remove all other interf (POINT,(not(ON,FACE)))
    let mut n = 0;
    while n < li.len() {
        let (is_point, on_face) = {
            let i = li[n].read().expect("interference lock");
            (
                i.geometry_type() == Kind::Point,
                fun_has_state_shape(&i.transition(), State::On, ShapeType::Face),
            )
        };
        if is_point {
            if !on_face {
                li.remove(n);
            } else {
                n += 1;
            }
        } else {
            n += 1;
        }
    }
}
