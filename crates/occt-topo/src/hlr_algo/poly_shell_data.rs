//! All the `HLRAlgo_PolyData` of one shell plus its edge list. Source:
//! `HLRAlgo_PolyShellData`
//! (`ModelingAlgorithms/TKHLR/HLRAlgo/HLRAlgo_PolyShellData.hxx:26-63`,
//! `HLRAlgo_PolyShellData.cxx:26-62`).

use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt;

use super::bi_point::HlrAlgoBiPoint;
use super::poly_data::HlrAlgoPolyData;

/// `HLRAlgo_PolyShellData::ShellIndices` (`hxx:31-34`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PolyShellIndices {
    pub min: i32,
    pub max: i32,
}

/// All the `PolyData` of a shell. Source: `HLRAlgo_PolyShellData`.
///
/// `HidingPolyData()` stores `handle<HLRAlgo_PolyData>`s in OCCT (`hxx:47`) and
/// therefore aliases elements of `PolyData()`. The port owns the poly data in
/// `poly_data` and keeps the *indices* of the hiding ones, which preserves the
/// element identity without reference counting.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct HlrAlgoPolyShellData {
    indices: PolyShellIndices,
    poly_data: Vec<HlrAlgoPolyData>,
    hiding_poly_data: Vec<usize>,
    edges: Vec<HlrAlgoBiPoint>,
}

impl HlrAlgoPolyShellData {
    /// `HLRAlgo_PolyShellData(const int nbFace)` (`cxx:28-30`): the poly data
    /// array is allocated but its slots stay null until filled.
    pub fn new(nb_face: usize) -> Self {
        Self {
            indices: PolyShellIndices::default(),
            poly_data: vec![HlrAlgoPolyData::new(); nb_face],
            hiding_poly_data: Vec::new(),
            edges: Vec::new(),
        }
    }

    /// `Hiding()` (`hxx:45`): `!myHPolHi.IsEmpty()`.
    pub fn hiding(&self) -> bool {
        !self.hiding_poly_data.is_empty()
    }

    /// `PolyData()` (`hxx:47`).
    pub fn poly_data(&self) -> &[HlrAlgoPolyData] {
        &self.poly_data
    }

    /// `PolyData()` mutable form (the OCCT `NCollection_Array1&` accessor).
    pub fn change_poly_data(&mut self) -> &mut Vec<HlrAlgoPolyData> {
        &mut self.poly_data
    }

    /// `HidingPolyData()` (`hxx:49`): the hiding poly data, resolved from the
    /// stored indices.
    pub fn hiding_poly_data(&self) -> impl Iterator<Item = &HlrAlgoPolyData> {
        self.hiding_poly_data.iter().map(|&i| &self.poly_data[i])
    }

    /// Number of entries of `HidingPolyData()`.
    pub fn nb_hiding_poly_data(&self) -> usize {
        self.hiding_poly_data.len()
    }

    /// `HidingPolyData().SetValue(i, pd)` (`PolyAlgo.cxx:452`).
    pub fn push_hiding_poly_data(&mut self, index: usize) {
        self.hiding_poly_data.push(index);
    }

    /// `Edges()` (`hxx:49`).
    pub fn edges(&self) -> &[HlrAlgoBiPoint] {
        &self.edges
    }

    /// `Edges()` mutable form.
    pub fn change_edges(&mut self) -> &mut Vec<HlrAlgoBiPoint> {
        &mut self.edges
    }

    /// `Indices()` (`hxx:51`).
    pub fn indices(&self) -> &PolyShellIndices {
        &self.indices
    }

    /// `Indices()` mutable form.
    pub fn change_indices(&mut self) -> &mut PolyShellIndices {
        &mut self.indices
    }

    /// `UpdateGlobalMinMax(theBox)` (`cxx:34-46`): every edge endpoint, then
    /// every poly data.
    pub fn update_global_min_max(&self, the_box: &mut BndBox) {
        for e in &self.edges {
            let p = e.points();
            the_box.add_point(&GpPnt::new(p.pntp1.x(), p.pntp1.y(), p.pntp1.z()));
            the_box.add_point(&GpPnt::new(p.pntp2.x(), p.pntp2.y(), p.pntp2.z()));
        }
        for pd in &self.poly_data {
            pd.update_global_min_max(the_box);
        }
    }

    /// `UpdateHiding(nbHiding)` (`cxx:52-62`): `nbHiding > 0` resizes the
    /// hiding array without copying, `0` moves in an empty array.
    pub fn update_hiding(&mut self, nb_hiding: usize) {
        self.hiding_poly_data = if nb_hiding > 0 {
            Vec::with_capacity(nb_hiding)
        } else {
            Vec::new()
        };
    }
}
