
use super::super::data_structure::DataStructure;
use super::super::interference::Interference;
use super::super::kind::Kind;

use super::*;

/// `FDS_SIisGIofIofSBAofTofI(BDS, SI, I)` (`cxx:396-462`): the face that `I`
/// names (before side first, then after side) carries an interference whose
/// geometry is the edge `SI`.
pub fn fds_si_is_gi_of_i_of_sba_of_t_of_i(
    bds: &DataStructure,
    si: i32,
    i: &Interference,
) -> bool {
    if si == 0 {
        return false;
    }
    let (sb1, ib1, sa1, ia1, _gt1, _g1, _st1, _s1) = fds_idata(i);

    let face_index = if sb1 == crate::abs::ShapeType::Face {
        ib1
    } else if sa1 == crate::abs::ShapeType::Face {
        ia1
    } else {
        return false;
    };

    for ib in bds.shape_interferences_at(face_index, true) {
        let (gtb, gb) = {
            let other = ib.read().expect("interference lock");
            let (_sbb, _ibb, _sab, _iab, gtb, gb, _stb, _sb) = fds_idata(&other);
            (gtb, gb)
        };
        if gtb == Kind::Edge && gb == si {
            // the face carries an interference whose geometry is the edge SI
            return true;
        }
    }
    false
}
