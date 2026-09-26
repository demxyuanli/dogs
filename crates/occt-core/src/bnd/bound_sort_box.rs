//! Port of Bnd_BoundSortBox and its private helper Bnd_VoxelGrid
//! (src/FoundationClasses/TKMath/Bnd/Bnd_BoundSortBox.hxx, 150 lines;
//!  src/FoundationClasses/TKMath/Bnd/Bnd_BoundSortBox.cxx, 624 lines), OCCT 8.0.0 (tag V8_0_0).
//!
//! Source: D:\source\OCCT-src\src\FoundationClasses\TKMath\Bnd\Bnd_BoundSortBox.{hxx,cxx}
//!
//! Coverage (OCCT line numbers):
//! * getBnd_VoxelGridResolution          cxx:102-121
//! * Bnd_VoxelGrid (class declaration)   cxx:129-210
//!   - ctor                              cxx:216-231
//!   - AddBox                            cxx:235-247
//!   - GetSliceX / GetSliceY / GetSliceZ cxx:251-271
//!   - AppendSliceX / Y / Z              cxx:275-307
//! * Bnd_BoundSortBox()                  cxx:311-320
//! * Initialize(setOfBoxes)              cxx:324-348
//! * Initialize(enclosing, setOfBoxes)   cxx:352-368
//! * Initialize(enclosing, nbBoxes)      cxx:372-390
//! * Add(box, index)                     cxx:394-406
//! * Compare(const Bnd_Box&)             cxx:410-504
//! * Compare(const gp_Pln&)              cxx:508-521
//! * calculateCoefficients               cxx:525-532
//! * resetVoxelGrid                      cxx:536-542
//! * sortBoxes                           cxx:546-552
//! * getBoundingVoxels                   cxx:556-596
//! * addBox                              cxx:600-624
//! * members myEnclosingBox .. myVoxelGrid  hxx:139-147
//!
//! Deviations from OCCT (none change the set returned by Compare):
//! * NCollection_DynamicArray / NCollection_IncAllocator
//!   (cxx:206-209, 222-230, 541) are replaced by Vec; the allocator
//!   increment is only a capacity hint.
//! * NCollection_HArray1<Bnd_Box> (hxx:140) is modelled as a 1-based
//!   Vec<BndBox> whose slot 0 is a void placeholder, so the indices returned
//!   by compare() are exactly OCCT's 1-based indices (cxx:396, 424, 494).
//!   OCCT also allows a caller-supplied array whose Lower is not 1; the port
//!   always uses 1..=len.
//! * Standard_MultiplyDefined / Standard_NullValue / Standard_OutOfRange
//!   (cxx:374, 396) become Rust assert! / slice-index panics.
//! * Initialize(setOfBoxes) does not reset myEnclosingBox before accumulating
//!   (cxx:328-335); the port reproduces that OCCT behaviour faithfully.
//!
//! UNPORTED / known divergence in this module:
//! * theBox.Get(...) (cxx:562) and Bnd_Box::IsOut(gp_Pln)
//!   (Bnd_Box.cxx:709-757) use this crate's BndBox::get, which returns
//!   +/-f64::INFINITY for open directions where OCCT returns
//!   +/-THE_BND_PRECISION_INFINITE = +/-1e100 (Bnd_Box.cxx:31, 238-274).
//!   For finite boxes the values are identical; for open boxes the crate's
//!   pre-existing representation governs. Likewise myEnclosingBox.CornerMin()
//!   (cxx:559): BndBox::corner_min ignores the open flags, where
//!   Bnd_Box::CornerMin (Bnd_Box.cxx:278-287) substitutes the infinite.
//! * the .hxx exposes no Boxes() / EnclosingBox() / Extent() getters in OCCT
//!   8.0.0; only the private members at hxx:139-147. Read-only accessors are
//!   provided at the end of this file and tagged with those member lines.

use crate::bnd::box3d::BndBox;
use crate::gp::{GpPln, GpPnt};

/// getBnd_VoxelGridResolution (cxx:102-121).
///
/// Resolution depends only on the number of boxes; it is 8, 16, 32, 64 or
/// 128 for 0..=100, 101..=1000, 1001..=10000, 10001..=40000 and > 40000 boxes.
fn get_bnd_voxel_grid_resolution(the_boxes_count: usize) -> usize {
    if the_boxes_count > 40000 {
        return 128;
    }
    if the_boxes_count > 10000 {
        return 64;
    }
    if the_boxes_count > 1000 {
        return 32;
    }
    if the_boxes_count > 100 {
        return 16;
    }
    8
}

/// C++ static_cast<int>(double): truncation toward zero.
///
/// Rust's f64 as i64 saturates on overflow and maps NaN to 0, where C++ is
/// undefined; every caller then clamps into [0, myResolution - 1]
/// (cxx:582-593), so the saturated value is still in range.
#[inline]
fn static_cast_int(v: f64) -> i64 {
    v as i64
}

/// std::clamp(index, 0, myResolution - 1) (cxx:582-593).
#[inline]
fn clamp_voxel_index(index: i64, the_resolution: usize) -> usize {
    index.clamp(0, the_resolution as i64 - 1) as usize
}

/// Bnd_VoxelGrid (cxx:129-210): three per-axis slice tables; slice i of each
/// table lists every box whose voxel extent covers i on that axis.
#[derive(Debug, Clone)]
struct BndVoxelGrid {
    /// mySlicesX (hxx-style member at cxx:207), one entry per X voxel.
    my_slices_x: Vec<Vec<usize>>,
    /// mySlicesY (cxx:208).
    my_slices_y: Vec<Vec<usize>>,
    /// mySlicesZ (cxx:209).
    my_slices_z: Vec<Vec<usize>>,
}

impl BndVoxelGrid {
    /// Bnd_VoxelGrid(theResolution, theExpectedBoxCount) ctor (cxx:216-231).
    /// theExpectedBoxCount only pre-sizes the per-slice vectors in OCCT
    /// (cxx:223-230); Vec grows on demand, so it is ignored.
    fn new(the_resolution: usize, _the_expected_box_count: usize) -> Self {
        Self {
            my_slices_x: vec![Vec::new(); the_resolution],
            my_slices_y: vec![Vec::new(); the_resolution],
            my_slices_z: vec![Vec::new(); the_resolution],
        }
    }

    /// AddBox(theBoxIndex, theVoxelBox) (cxx:235-247);
    /// theVoxelBox is [minX, minY, minZ, maxX, maxY, maxZ].
    fn add_box(&mut self, the_box_index: usize, the_voxel_box: [usize; 6]) {
        let [a_min_voxel_x, a_min_voxel_y, a_min_voxel_z, a_max_voxel_x, a_max_voxel_y, a_max_voxel_z] =
            the_voxel_box;
        append_slice(&mut self.my_slices_x, a_min_voxel_x, a_max_voxel_x, the_box_index);
        append_slice(&mut self.my_slices_y, a_min_voxel_y, a_max_voxel_y, the_box_index);
        append_slice(&mut self.my_slices_z, a_min_voxel_z, a_max_voxel_z, the_box_index);
    }

    /// GetSliceX(theVoxelIndex) (cxx:251-255): None stands for the null
    /// pointer OCCT returns for an empty slice.
    fn get_slice_x(&self, the_voxel_index: usize) -> Option<&[usize]> {
        let a_slice = &self.my_slices_x[the_voxel_index];
        if a_slice.is_empty() { None } else { Some(a_slice) }
    }

    /// GetSliceY(theVoxelIndex) (cxx:259-263).
    fn get_slice_y(&self, the_voxel_index: usize) -> Option<&[usize]> {
        let a_slice = &self.my_slices_y[the_voxel_index];
        if a_slice.is_empty() { None } else { Some(a_slice) }
    }

    /// GetSliceZ(theVoxelIndex) (cxx:267-271).
    fn get_slice_z(&self, the_voxel_index: usize) -> Option<&[usize]> {
        let a_slice = &self.my_slices_z[the_voxel_index];
        if a_slice.is_empty() { None } else { Some(a_slice) }
    }
}

/// AppendSliceX / AppendSliceY / AppendSliceZ (cxx:275-307): append the box
/// index to every voxel in [theVoxelIndexMin, theVoxelIndexMax].
fn append_slice(the_slices: &mut [Vec<usize>], the_index_min: usize, the_index_max: usize, the_box_index: usize) {
    for i in the_index_min..=the_index_max {
        the_slices[i].push(the_box_index);
    }
}

/// Bnd_BoundSortBox (hxx:51-148). Sorts a set of bounding boxes into a voxel
/// grid so that Compare returns the boxes touched by a query box/plane.
#[derive(Debug, Clone, Default)]
pub struct BoundSortBox {
    /// myEnclosingBox (hxx:139).
    my_enclosing_box: BndBox,
    /// myBoxes (hxx:140), 1-based: my_boxes[i] is OCCT index i;
    /// slot 0 is an unused void placeholder. Length == Upper + 1.
    my_boxes: Vec<BndBox>,
    /// myCoeffX (hxx:141).
    my_coeff_x: f64,
    /// myCoeffY (hxx:142).
    my_coeff_y: f64,
    /// myCoeffZ (hxx:143).
    my_coeff_z: f64,
    /// myResolution (hxx:144).
    my_resolution: usize,
    /// myLastResult (hxx:145); returned by Compare.
    my_last_result: Vec<usize>,
    /// myLargeBoxes (hxx:146): boxes too large to voxelize, tested directly.
    my_large_boxes: Vec<usize>,
    /// myVoxelGrid (hxx:147); None before Initialize/resetVoxelGrid, as in
    /// OCCT's null handle.
    my_voxel_grid: Option<BndVoxelGrid>,
}

impl BoundSortBox {
    /// Bnd_BoundSortBox() (cxx:311-320).
    pub fn new() -> Self {
        Self::default()
    }

    /// Initialize(const Handle(NCollection_HArray1<Bnd_Box>)&) (cxx:324-348).
    ///
    /// The enclosing box is the union of the non-void boxes. Note that, as in
    /// OCCT, this does not reset myEnclosingBox first (cxx:328-335), so calling
    /// it twice accumulates onto the previous enclosing box.
    pub fn initialize(&mut self, the_set_of_boxes: &[BndBox]) {
        self.my_boxes = with_placeholder(the_set_of_boxes);

        for a_box in the_set_of_boxes {
            // if (!aBox.IsVoid()) myEnclosingBox.Add(aBox) (cxx:330-335).
            if !a_box.is_void() {
                self.my_enclosing_box.add_box(a_box);
            }
        }

        // myResolution = getBnd_VoxelGridResolution(myBoxes->Length()) (cxx:337).
        self.my_resolution = get_bnd_voxel_grid_resolution(self.my_boxes.len() - 1);

        if self.my_enclosing_box.is_void() {
            return; // cxx:339-342
        }

        self.calculate_coefficients(); // cxx:344
        self.reset_voxel_grid();       // cxx:345
        self.sort_boxes();             // cxx:347
    }

    /// Initialize(const Bnd_Box& theEnclosingBox,
    ///            const Handle(NCollection_HArray1<Bnd_Box>)&) (cxx:352-368).
    pub fn initialize_with_enclosing(&mut self, the_enclosing_box: &BndBox, the_set_of_boxes: &[BndBox]) {
        self.my_boxes = with_placeholder(the_set_of_boxes);      // cxx:355
        self.my_enclosing_box = *the_enclosing_box;              // cxx:356
        self.my_resolution = get_bnd_voxel_grid_resolution(self.my_boxes.len() - 1); // cxx:357

        if self.my_enclosing_box.is_void() {
            return; // cxx:359-362
        }

        self.calculate_coefficients(); // cxx:364
        self.reset_voxel_grid();       // cxx:365
        self.sort_boxes();             // cxx:367
    }

    /// Initialize(const Bnd_Box& theEnclosingBox, const int theNbBoxes)
    /// (cxx:372-390). Boxes are then supplied with add().
    pub fn initialize_empty(&mut self, the_enclosing_box: &BndBox, the_nb_boxes: usize) {
        // Standard_NullValue_Raise_if(theNbBoxes <= 0, ...) (cxx:374).
        assert!(the_nb_boxes > 0, "Unexpected: theNbBoxes <= 0");
        // new NCollection_HArray1<Bnd_Box>(1, theNbBoxes) + Init(emptyBox)
        // (cxx:375-379): slot 0 is the placeholder, 1..=theNbBoxes are void.
        self.my_boxes = vec![BndBox::new(); the_nb_boxes + 1];
        self.my_enclosing_box = *the_enclosing_box; // cxx:380
        self.my_resolution = get_bnd_voxel_grid_resolution(the_nb_boxes); // cxx:381

        if self.my_enclosing_box.is_void() {
            return; // cxx:383-386
        }

        self.calculate_coefficients(); // cxx:388
        self.reset_voxel_grid();       // cxx:389
    }

    /// Add(const Bnd_Box& theBox, const int theIndex) (cxx:394-406).
    ///
    /// theIndex is 1-based and must be in [1, nbBoxes] from initialize_empty.
    /// Panics like Standard_MultiplyDefined if the slot already holds a box
    /// (cxx:396-397).
    pub fn add(&mut self, the_box: &BndBox, the_index: usize) {
        assert!(
            the_index >= 1 && the_index < self.my_boxes.len(),
            "Standard_OutOfRange: index {the_index} not in [1, {}]",
            self.my_boxes.len().saturating_sub(1)
        );
        // Standard_MultiplyDefined_Raise_if(!(myBoxes->Value(theIndex).IsVoid()), ...)
        // (cxx:396-397).
        assert!(self.my_boxes[the_index].is_void(), " This box is already defined !");
        if the_box.is_void() {
            return; // cxx:398-401
        }
        self.my_boxes[the_index] = *the_box; // cxx:403
        self.add_box(the_box, the_index);    // cxx:405
    }

    /// Compare(const Bnd_Box& theBox) (cxx:410-504): the indices of all stored
    /// boxes touched by theBox (i.e. not IsOut(theBox)), in OCCT's order:
    /// large boxes first (cxx:422-429), then boxes found through the Z slices
    /// (cxx:479-501).
    pub fn compare(&mut self, the_box: &BndBox) -> &[usize] {
        self.my_last_result.clear(); // cxx:413

        // cxx:415-418.
        if the_box.is_void() || the_box.is_out_box(&self.my_enclosing_box) {
            return &self.my_last_result;
        }

        // Processing the large boxes (cxx:422-429).
        for i in 0..self.my_large_boxes.len() {
            let a_box_index = self.my_large_boxes[i];
            if !self.my_boxes[a_box_index].is_out_box(the_box) {
                self.my_last_result.push(a_box_index);
            }
        }

        // Obtaining the box voxel coordinates (cxx:432-433).
        let a_voxels = self.get_bounding_voxels(the_box);
        let [a_min_voxel_x, a_min_voxel_y, a_min_voxel_z, a_max_voxel_x, a_max_voxel_y, a_max_voxel_z] =
            a_voxels;

        // Bit mask per box index (cxx:442-446).
        let mut a_result_indices = vec![0u8; self.my_boxes.len()]; // Resize(Upper+1, 0) cxx:443
        const AN_OCCUPIED_X: u8 = 0b01;
        const AN_OCCUPIED_Y: u8 = 0b10;
        const AN_OCCUPIED_XY: u8 = 0b11;

        if let Some(a_voxel_grid) = self.my_voxel_grid.as_ref() {
            // Checking the voxels along X-axis (cxx:449-461).
            for a_voxel_x in a_min_voxel_x..=a_max_voxel_x {
                if let Some(a_box_indices) = a_voxel_grid.get_slice_x(a_voxel_x) {
                    for &a_box_index in a_box_indices {
                        a_result_indices[a_box_index] |= AN_OCCUPIED_X;
                    }
                }
            }

            // Checking the voxels along Y-axis (cxx:464-476).
            for a_voxel_y in a_min_voxel_y..=a_max_voxel_y {
                if let Some(a_box_indices) = a_voxel_grid.get_slice_y(a_voxel_y) {
                    for &a_box_index in a_box_indices {
                        a_result_indices[a_box_index] |= AN_OCCUPIED_Y;
                    }
                }
            }

            // Checking the voxels along Z-axis (cxx:479-501).
            for a_voxel_z in a_min_voxel_z..=a_max_voxel_z {
                if let Some(a_box_indices) = a_voxel_grid.get_slice_z(a_voxel_z) {
                    for &a_box_index in a_box_indices {
                        if a_result_indices[a_box_index] == AN_OCCUPIED_XY {
                            a_result_indices[a_box_index] = 0; // cxx:492
                            if !self.my_boxes[a_box_index].is_out_box(the_box) {
                                self.my_last_result.push(a_box_index); // cxx:497
                            }
                        }
                    }
                }
            }
        }

        &self.my_last_result // cxx:503
    }

    /// Compare(const gp_Pln& thePlane) (cxx:508-521): the indices of all stored
    /// boxes not IsOut(thePlane).
    pub fn compare_plane(&mut self, the_plane: &GpPln) -> &[usize] {
        self.my_last_result.clear(); // cxx:511
        for a_box_index in 1..self.my_boxes.len() {
            let a_box = &self.my_boxes[a_box_index];
            if !is_out_plane(a_box, the_plane) {
                self.my_last_result.push(a_box_index); // cxx:517
            }
        }
        &self.my_last_result // cxx:520
    }

    /// calculateCoefficients() (cxx:525-532).
    fn calculate_coefficients(&mut self) {
        // BndBox::get returns the tuple (xmin, xmax, ymin, ymax, zmin, zmax)
        // (box3d.rs:139), not the out-parameter order of Bnd_Box::Get
        // (xmin, ymin, zmin, xmax, ymax, zmax) at cxx:528; rebind accordingly.
        let (a_xmin, a_xmax, a_ymin, a_ymax, a_zmin, a_zmax) =
            self.my_enclosing_box.get().expect("Bnd_Box is void"); // cxx:528
        let a_resolution = self.my_resolution as f64;
        self.my_coeff_x = if a_xmax - a_xmin == 0.0 { 0.0 } else { a_resolution / (a_xmax - a_xmin) }; // cxx:529
        self.my_coeff_y = if a_ymax - a_ymin == 0.0 { 0.0 } else { a_resolution / (a_ymax - a_ymin) }; // cxx:530
        self.my_coeff_z = if a_zmax - a_zmin == 0.0 { 0.0 } else { a_resolution / (a_zmax - a_zmin) }; // cxx:531
    }

    /// resetVoxelGrid() (cxx:536-542).
    fn reset_voxel_grid(&mut self) {
        self.my_voxel_grid = Some(BndVoxelGrid::new(self.my_resolution, self.my_boxes.len() - 1)); // cxx:538
        self.my_large_boxes.clear(); // cxx:539
        // myLargeBoxes.SetIncrement(max(Length / 16, 16)) (cxx:541) is a
        // capacity hint for NCollection_DynamicArray; Vec grows as needed.
    }

    /// sortBoxes() (cxx:546-552).
    fn sort_boxes(&mut self) {
        for a_box_index in 1..self.my_boxes.len() {
            let a_box = self.my_boxes[a_box_index];
            self.add_box(&a_box, a_box_index); // cxx:550
        }
    }

    /// getBoundingVoxels(const Bnd_Box&) (cxx:556-596):
    /// [minX, minY, minZ, maxX, maxY, maxZ] voxel indices.
    fn get_bounding_voxels(&self, the_box: &BndBox) -> [usize; 6] {
        // Start point of the voxel grid (cxx:559).
        let a_grid_start: GpPnt = self.my_enclosing_box.corner_min();

        // Same tuple-order note as in calculate_coefficients (box3d.rs:139).
        let (a_xmin, a_xmax, a_ymin, a_ymax, a_zmin, a_zmax) =
            the_box.get().expect("Bnd_Box is void"); // cxx:562
        let a_resolution = self.my_resolution;

        // cxx:582-593. The -1 / +1 safety margin is OCCT's.
        let a_xmin_index = clamp_voxel_index(
            static_cast_int((a_xmin - a_grid_start.x()) * self.my_coeff_x) - 1,
            a_resolution,
        );
        let a_ymin_index = clamp_voxel_index(
            static_cast_int((a_ymin - a_grid_start.y()) * self.my_coeff_y) - 1,
            a_resolution,
        );
        let a_zmin_index = clamp_voxel_index(
            static_cast_int((a_zmin - a_grid_start.z()) * self.my_coeff_z) - 1,
            a_resolution,
        );
        let a_xmax_index = clamp_voxel_index(
            static_cast_int((a_xmax - a_grid_start.x()) * self.my_coeff_x) + 1,
            a_resolution,
        );
        let a_ymax_index = clamp_voxel_index(
            static_cast_int((a_ymax - a_grid_start.y()) * self.my_coeff_y) + 1,
            a_resolution,
        );
        let a_zmax_index = clamp_voxel_index(
            static_cast_int((a_zmax - a_grid_start.z()) * self.my_coeff_z) + 1,
            a_resolution,
        );

        [
            a_xmin_index, a_ymin_index, a_zmin_index, a_xmax_index, a_ymax_index, a_zmax_index,
        ] // cxx:595
    }

    /// addBox(const Bnd_Box& theBox, const int theIndex) (cxx:600-624).
    fn add_box(&mut self, the_box: &BndBox, the_index: usize) {
        if the_box.is_void() {
            return; // cxx:602-605
        }

        let a_voxel_box = self.get_bounding_voxels(the_box); // cxx:607-608
        let [a_min_voxel_x, a_min_voxel_y, a_min_voxel_z, a_max_voxel_x, a_max_voxel_y, a_max_voxel_z] =
            a_voxel_box;
        let a_box_min_side = (a_max_voxel_x - a_min_voxel_x)
            .min(a_max_voxel_y - a_min_voxel_y)
            .min(a_max_voxel_z - a_min_voxel_z); // cxx:610-611

        if a_box_min_side * 4 > self.my_resolution {
            self.my_large_boxes.push(the_index); // cxx:617
        } else if let Some(a_voxel_grid) = self.my_voxel_grid.as_mut() {
            a_voxel_grid.add_box(the_index, a_voxel_box); // cxx:621-622
        }
    }

    // ---- read-only accessors for the hxx:139-147 members ----

    /// myBoxes, values 1..=Upper (hxx:140).
    pub fn boxes(&self) -> &[BndBox] {
        if self.my_boxes.is_empty() {
            &[]
        } else {
            &self.my_boxes[1..]
        }
    }

    /// myEnclosingBox (hxx:139).
    pub fn enclosing_box(&self) -> &BndBox {
        &self.my_enclosing_box
    }

    /// myCoeffX (hxx:141).
    pub fn coeff_x(&self) -> f64 {
        self.my_coeff_x
    }

    /// myCoeffY (hxx:142).
    pub fn coeff_y(&self) -> f64 {
        self.my_coeff_y
    }

    /// myCoeffZ (hxx:143).
    pub fn coeff_z(&self) -> f64 {
        self.my_coeff_z
    }

    /// myResolution (hxx:144).
    pub fn resolution(&self) -> usize {
        self.my_resolution
    }

    /// myLastResult (hxx:145).
    pub fn last_result(&self) -> &[usize] {
        &self.my_last_result
    }

    /// myLargeBoxes (hxx:146).
    pub fn large_boxes(&self) -> &[usize] {
        &self.my_large_boxes
    }
}

/// Build the 1-based box vector: slot 0 is the void placeholder, so that
/// my_boxes[i] carries OCCT index i.
fn with_placeholder(the_set_of_boxes: &[BndBox]) -> Vec<BndBox> {
    let mut a_boxes = Vec::with_capacity(the_set_of_boxes.len() + 1);
    a_boxes.push(BndBox::new());
    a_boxes.extend_from_slice(the_set_of_boxes);
    a_boxes
}

/// Bnd_Box::IsOut(const gp_Pln&) (Bnd_Box.cxx:709-757): true when the box lies
/// strictly on one side of the plane.
fn is_out_plane(the_box: &BndBox, the_plane: &GpPln) -> bool {
    if the_box.is_whole() {
        return false; // Bnd_Box.cxx:711-714
    }
    if the_box.is_void() {
        return true; // Bnd_Box.cxx:715-718
    }
    let (a, b, c, d) = the_plane.coefficients(); // Bnd_Box.cxx:721-722
    // get() order is (xmin, xmax, ymin, ymax, zmin, zmax) (box3d.rs:139).
    let (a_xmin, a_xmax, a_ymin, a_ymax, a_zmin, a_zmax) = the_box.get().unwrap(); // Bnd_Box.cxx:723-728
    let f = |x: f64, y: f64, z: f64| a * x + b * y + c * z + d;
    let plus = f(a_xmin, a_ymin, a_zmin) > 0.0; // Bnd_Box.cxx:729-730
    if plus != (f(a_xmin, a_ymin, a_zmax) > 0.0) {
        return false; // Bnd_Box.cxx:731-734
    }
    if plus != (f(a_xmin, a_ymax, a_zmin) > 0.0) {
        return false; // Bnd_Box.cxx:735-738
    }
    if plus != (f(a_xmin, a_ymax, a_zmax) > 0.0) {
        return false; // Bnd_Box.cxx:739-742
    }
    if plus != (f(a_xmax, a_ymin, a_zmin) > 0.0) {
        return false; // Bnd_Box.cxx:743-746
    }
    if plus != (f(a_xmax, a_ymin, a_zmax) > 0.0) {
        return false; // Bnd_Box.cxx:747-750
    }
    if plus != (f(a_xmax, a_ymax, a_zmin) > 0.0) {
        return false; // Bnd_Box.cxx:751-754
    }
    plus == (f(a_xmax, a_ymax, a_zmax) > 0.0) // Bnd_Box.cxx:755
}
