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
fn get_bnd_voxel_grid_resolution(boxes_count: usize) -> usize {
    if boxes_count > 40000 {
        return 128;
    }
    if boxes_count > 10000 {
        return 64;
    }
    if boxes_count > 1000 {
        return 32;
    }
    if boxes_count > 100 {
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
fn clamp_voxel_index(index: i64, resolution: usize) -> usize {
    if resolution == 0 {
        return 0;
    }
    index.clamp(0, resolution as i64 - 1) as usize
}

/// Bnd_VoxelGrid (cxx:129-210): three per-axis slice tables; slice i of each
/// table lists every box whose voxel extent covers i on that axis.
#[derive(Debug, Clone)]
struct BndVoxelGrid {
    /// mySlicesX (hxx-style member at cxx:207), one entry per X voxel.
    slices_x: Vec<Vec<usize>>,
    /// mySlicesY (cxx:208).
    slices_y: Vec<Vec<usize>>,
    /// mySlicesZ (cxx:209).
    slices_z: Vec<Vec<usize>>,
}

impl BndVoxelGrid {
    /// Bnd_VoxelGrid(theResolution, theExpectedBoxCount) ctor (cxx:216-231).
    /// theExpectedBoxCount only pre-sizes the per-slice vectors in OCCT
    /// (cxx:223-230); Vec grows on demand, so it is ignored.
    fn new(resolution: usize, _the_expected_box_count: usize) -> Self {
        Self {
            slices_x: vec![Vec::new(); resolution],
            slices_y: vec![Vec::new(); resolution],
            slices_z: vec![Vec::new(); resolution],
        }
    }

    /// AddBox(theBoxIndex, theVoxelBox) (cxx:235-247);
    /// theVoxelBox is [minX, minY, minZ, maxX, maxY, maxZ].
    fn add_box(&mut self, box_index: usize, voxel_box: [usize; 6]) {
        let [min_voxel_x, min_voxel_y, min_voxel_z, max_voxel_x, max_voxel_y, max_voxel_z] =
            voxel_box;
        append_slice(&mut self.slices_x, min_voxel_x, max_voxel_x, box_index);
        append_slice(&mut self.slices_y, min_voxel_y, max_voxel_y, box_index);
        append_slice(&mut self.slices_z, min_voxel_z, max_voxel_z, box_index);
    }

    /// GetSliceX(theVoxelIndex) (cxx:251-255): None stands for the null
    /// pointer OCCT returns for an empty slice.
    fn get_slice_x(&self, voxel_index: usize) -> Option<&[usize]> {
        let slice = &self.slices_x[voxel_index];
        if slice.is_empty() { None } else { Some(slice) }
    }

    /// GetSliceY(theVoxelIndex) (cxx:259-263).
    fn get_slice_y(&self, voxel_index: usize) -> Option<&[usize]> {
        let slice = &self.slices_y[voxel_index];
        if slice.is_empty() { None } else { Some(slice) }
    }

    /// GetSliceZ(theVoxelIndex) (cxx:267-271).
    fn get_slice_z(&self, voxel_index: usize) -> Option<&[usize]> {
        let slice = &self.slices_z[voxel_index];
        if slice.is_empty() { None } else { Some(slice) }
    }
}

/// AppendSliceX / AppendSliceY / AppendSliceZ (cxx:275-307): append the box
/// index to every voxel in [theVoxelIndexMin, theVoxelIndexMax].
fn append_slice(slices: &mut [Vec<usize>], index_min: usize, index_max: usize, box_index: usize) {
    for i in index_min..=index_max {
        slices[i].push(box_index);
    }
}

/// Bnd_BoundSortBox (hxx:51-148). Sorts a set of bounding boxes into a voxel
/// grid so that Compare returns the boxes touched by a query box/plane.
#[derive(Debug, Clone, Default)]
pub struct BoundSortBox {
    /// myEnclosingBox (hxx:139).
    enclosing_box: BndBox,
    /// myBoxes (hxx:140), 1-based: boxes[i] is OCCT index i;
    /// slot 0 is an unused void placeholder. Length == Upper + 1.
    boxes: Vec<BndBox>,
    /// myCoeffX (hxx:141).
    coeff_x: f64,
    /// myCoeffY (hxx:142).
    coeff_y: f64,
    /// myCoeffZ (hxx:143).
    coeff_z: f64,
    /// myResolution (hxx:144).
    resolution: usize,
    /// myLastResult (hxx:145); returned by Compare.
    last_result: Vec<usize>,
    /// myLargeBoxes (hxx:146): boxes too large to voxelize, tested directly.
    large_boxes: Vec<usize>,
    /// myVoxelGrid (hxx:147); None before Initialize/resetVoxelGrid, as in
    /// OCCT's null handle.
    voxel_grid: Option<BndVoxelGrid>,
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
    pub fn initialize(&mut self, set_of_boxes: &[BndBox]) {
        self.boxes = with_placeholder(set_of_boxes);

        for bnd_box in set_of_boxes {
            // if (!aBox.IsVoid()) myEnclosingBox.Add(aBox) (cxx:330-335).
            if !bnd_box.is_void() {
                self.enclosing_box.add_box(bnd_box);
            }
        }

        // myResolution = getBnd_VoxelGridResolution(myBoxes->Length()) (cxx:337).
        self.resolution = get_bnd_voxel_grid_resolution(self.boxes.len() - 1);

        if self.enclosing_box.is_void() {
            return; // cxx:339-342
        }

        self.calculate_coefficients(); // cxx:344
        self.reset_voxel_grid();       // cxx:345
        self.sort_boxes();             // cxx:347
    }

    /// Initialize(const Bnd_Box& theEnclosingBox,
    ///            const Handle(NCollection_HArray1<Bnd_Box>)&) (cxx:352-368).
    pub fn initialize_with_enclosing(&mut self, the_enclosing_box: &BndBox, set_of_boxes: &[BndBox]) {
        self.boxes = with_placeholder(set_of_boxes);      // cxx:355
        self.enclosing_box = *the_enclosing_box;              // cxx:356
        self.resolution = get_bnd_voxel_grid_resolution(self.boxes.len() - 1); // cxx:357

        if self.enclosing_box.is_void() {
            return; // cxx:359-362
        }

        self.calculate_coefficients(); // cxx:364
        self.reset_voxel_grid();       // cxx:365
        self.sort_boxes();             // cxx:367
    }

    /// Initialize(const Bnd_Box& theEnclosingBox, const int theNbBoxes)
    /// (cxx:372-390). Boxes are then supplied with add().
    pub fn initialize_empty(&mut self, the_enclosing_box: &BndBox, nb_boxes: usize) {
        // Standard_NullValue_Raise_if(theNbBoxes <= 0, ...) (cxx:374).
        assert!(nb_boxes > 0, "Unexpected: theNbBoxes <= 0");
        // new NCollection_HArray1<Bnd_Box>(1, theNbBoxes) + Init(emptyBox)
        // (cxx:375-379): slot 0 is the placeholder, 1..=theNbBoxes are void.
        self.boxes = vec![BndBox::new(); nb_boxes + 1];
        self.enclosing_box = *the_enclosing_box; // cxx:380
        self.resolution = get_bnd_voxel_grid_resolution(nb_boxes); // cxx:381

        if self.enclosing_box.is_void() {
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
    pub fn add(&mut self, bnd_box: &BndBox, index: usize) {
        assert!(
            index >= 1 && index < self.boxes.len(),
            "Standard_OutOfRange: index {index} not in [1, {}]",
            self.boxes.len().saturating_sub(1)
        );
        // Standard_MultiplyDefined_Raise_if(!(myBoxes->Value(theIndex).IsVoid()), ...)
        // (cxx:396-397).
        assert!(self.boxes[index].is_void(), " This box is already defined !");
        if bnd_box.is_void() {
            return; // cxx:398-401
        }
        self.boxes[index] = *bnd_box; // cxx:403
        self.add_box(bnd_box, index);    // cxx:405
    }

    /// Compare(const Bnd_Box& theBox) (cxx:410-504): the indices of all stored
    /// boxes touched by theBox (i.e. not IsOut(theBox)), in OCCT's order:
    /// large boxes first (cxx:422-429), then boxes found through the Z slices
    /// (cxx:479-501).
    pub fn compare(&mut self, bnd_box: &BndBox) -> &[usize] {
        self.last_result.clear(); // cxx:413

        // cxx:415-418.
        if bnd_box.is_void() || bnd_box.is_out_box(&self.enclosing_box) {
            return &self.last_result;
        }

        // Processing the large boxes (cxx:422-429).
        for i in 0..self.large_boxes.len() {
            let box_index = self.large_boxes[i];
            if !self.boxes[box_index].is_out_box(bnd_box) {
                self.last_result.push(box_index);
            }
        }

        // Obtaining the box voxel coordinates (cxx:432-433).
        let voxels = self.get_bounding_voxels(bnd_box);
        let [min_voxel_x, min_voxel_y, min_voxel_z, max_voxel_x, max_voxel_y, max_voxel_z] =
            voxels;

        // Bit mask per box index (cxx:442-446).
        let mut result_indices = vec![0u8; self.boxes.len()]; // Resize(Upper+1, 0) cxx:443
        const AN_OCCUPIED_X: u8 = 0b01;
        const AN_OCCUPIED_Y: u8 = 0b10;
        const AN_OCCUPIED_XY: u8 = 0b11;

        if let Some(a_voxel_grid) = self.voxel_grid.as_ref() {
            // Checking the voxels along X-axis (cxx:449-461).
            for voxel_x in min_voxel_x..=max_voxel_x {
                if let Some(box_indices) = a_voxel_grid.get_slice_x(voxel_x) {
                    for &box_index in box_indices {
                        result_indices[box_index] |= AN_OCCUPIED_X;
                    }
                }
            }

            // Checking the voxels along Y-axis (cxx:464-476).
            for voxel_y in min_voxel_y..=max_voxel_y {
                if let Some(box_indices) = a_voxel_grid.get_slice_y(voxel_y) {
                    for &box_index in box_indices {
                        result_indices[box_index] |= AN_OCCUPIED_Y;
                    }
                }
            }

            // Checking the voxels along Z-axis (cxx:479-501).
            for voxel_z in min_voxel_z..=max_voxel_z {
                if let Some(box_indices) = a_voxel_grid.get_slice_z(voxel_z) {
                    for &box_index in box_indices {
                        if result_indices[box_index] == AN_OCCUPIED_XY {
                            result_indices[box_index] = 0; // cxx:492
                            if !self.boxes[box_index].is_out_box(bnd_box) {
                                self.last_result.push(box_index); // cxx:497
                            }
                        }
                    }
                }
            }
        }

        &self.last_result // cxx:503
    }

    /// Compare(const gp_Pln& thePlane) (cxx:508-521): the indices of all stored
    /// boxes not IsOut(thePlane).
    pub fn compare_plane(&mut self, plane: &GpPln) -> &[usize] {
        self.last_result.clear(); // cxx:511
        for box_index in 1..self.boxes.len() {
            let bnd_box = &self.boxes[box_index];
            if !is_out_plane(bnd_box, plane) {
                self.last_result.push(box_index); // cxx:517
            }
        }
        &self.last_result // cxx:520
    }

    /// calculateCoefficients() (cxx:525-532).
    fn calculate_coefficients(&mut self) {
        // BndBox::get returns the tuple (xmin, xmax, ymin, ymax, zmin, zmax)
        // (box3d.rs:139), not the out-parameter order of Bnd_Box::Get
        // (xmin, ymin, zmin, xmax, ymax, zmax) at cxx:528; rebind accordingly.
        let (xmin, xmax, ymin, ymax, zmin, zmax) =
            self.enclosing_box.get().expect("Bnd_Box is void"); // cxx:528
        let a_resolution = self.resolution as f64;
        self.coeff_x = if xmax - xmin == 0.0 { 0.0 } else { a_resolution / (xmax - xmin) }; // cxx:529
        self.coeff_y = if ymax - ymin == 0.0 { 0.0 } else { a_resolution / (ymax - ymin) }; // cxx:530
        self.coeff_z = if zmax - zmin == 0.0 { 0.0 } else { a_resolution / (zmax - zmin) }; // cxx:531
    }

    /// resetVoxelGrid() (cxx:536-542).
    fn reset_voxel_grid(&mut self) {
        self.voxel_grid = Some(BndVoxelGrid::new(self.resolution, self.boxes.len() - 1)); // cxx:538
        self.large_boxes.clear(); // cxx:539
        // myLargeBoxes.SetIncrement(max(Length / 16, 16)) (cxx:541) is a
        // capacity hint for NCollection_DynamicArray; Vec grows as needed.
    }

    /// sortBoxes() (cxx:546-552).
    fn sort_boxes(&mut self) {
        for box_index in 1..self.boxes.len() {
            let bnd_box = self.boxes[box_index];
            self.add_box(&bnd_box, box_index); // cxx:550
        }
    }

    /// getBoundingVoxels(const Bnd_Box&) (cxx:556-596):
    /// [minX, minY, minZ, maxX, maxY, maxZ] voxel indices.
    fn get_bounding_voxels(&self, bnd_box: &BndBox) -> [usize; 6] {
        // Start point of the voxel grid (cxx:559).
        let grid_start: GpPnt = self.enclosing_box.corner_min();

        // Same tuple-order note as in calculate_coefficients (box3d.rs:139).
        let (xmin, xmax, ymin, ymax, zmin, zmax) =
            bnd_box.get().expect("Bnd_Box is void"); // cxx:562
        let a_resolution = self.resolution;

        // cxx:582-593. The -1 / +1 safety margin is OCCT's. `static_cast<int>`
        // of a coordinate far outside the grid saturates, and the following
        // `+/- 1` is then `std::clamp`ed into `[0, myResolution - 1]`. Saturating
        // arithmetic keeps that clamp reachable: a wrapping `i64` subtract
        // overflows before `clamp_voxel_index` runs.
        let xmin_index = clamp_voxel_index(
            static_cast_int((xmin - grid_start.x()) * self.coeff_x).saturating_sub(1),
            a_resolution,
        );
        let ymin_index = clamp_voxel_index(
            static_cast_int((ymin - grid_start.y()) * self.coeff_y).saturating_sub(1),
            a_resolution,
        );
        let zmin_index = clamp_voxel_index(
            static_cast_int((zmin - grid_start.z()) * self.coeff_z).saturating_sub(1),
            a_resolution,
        );
        let xmax_index = clamp_voxel_index(
            static_cast_int((xmax - grid_start.x()) * self.coeff_x).saturating_add(1),
            a_resolution,
        );
        let ymax_index = clamp_voxel_index(
            static_cast_int((ymax - grid_start.y()) * self.coeff_y).saturating_add(1),
            a_resolution,
        );
        let zmax_index = clamp_voxel_index(
            static_cast_int((zmax - grid_start.z()) * self.coeff_z).saturating_add(1),
            a_resolution,
        );

        [
            xmin_index, ymin_index, zmin_index, xmax_index, ymax_index, zmax_index,
        ] // cxx:595
    }

    /// addBox(const Bnd_Box& theBox, const int theIndex) (cxx:600-624).
    fn add_box(&mut self, bnd_box: &BndBox, index: usize) {
        if bnd_box.is_void() {
            return; // cxx:602-605
        }

        let voxel_box = self.get_bounding_voxels(bnd_box); // cxx:607-608
        let [min_voxel_x, min_voxel_y, min_voxel_z, max_voxel_x, max_voxel_y, max_voxel_z] =
            voxel_box;
        let box_min_side = (max_voxel_x - min_voxel_x)
            .min(max_voxel_y - min_voxel_y)
            .min(max_voxel_z - min_voxel_z); // cxx:610-611

        if box_min_side * 4 > self.resolution {
            self.large_boxes.push(index); // cxx:617
        } else if let Some(a_voxel_grid) = self.voxel_grid.as_mut() {
            a_voxel_grid.add_box(index, voxel_box); // cxx:621-622
        }
    }

    // ---- read-only accessors for the hxx:139-147 members ----

    /// myBoxes, values 1..=Upper (hxx:140).
    pub fn boxes(&self) -> &[BndBox] {
        if self.boxes.is_empty() {
            &[]
        } else {
            &self.boxes[1..]
        }
    }

    /// myEnclosingBox (hxx:139).
    pub fn enclosing_box(&self) -> &BndBox {
        &self.enclosing_box
    }

    /// myCoeffX (hxx:141).
    pub fn coeff_x(&self) -> f64 {
        self.coeff_x
    }

    /// myCoeffY (hxx:142).
    pub fn coeff_y(&self) -> f64 {
        self.coeff_y
    }

    /// myCoeffZ (hxx:143).
    pub fn coeff_z(&self) -> f64 {
        self.coeff_z
    }

    /// myResolution (hxx:144).
    pub fn resolution(&self) -> usize {
        self.resolution
    }

    /// myLastResult (hxx:145).
    pub fn last_result(&self) -> &[usize] {
        &self.last_result
    }

    /// myLargeBoxes (hxx:146).
    pub fn large_boxes(&self) -> &[usize] {
        &self.large_boxes
    }
}

/// Build the 1-based box vector: slot 0 is the void placeholder, so that
/// boxes[i] carries OCCT index i.
fn with_placeholder(set_of_boxes: &[BndBox]) -> Vec<BndBox> {
    let mut boxes = Vec::with_capacity(set_of_boxes.len() + 1);
    boxes.push(BndBox::new());
    boxes.extend_from_slice(set_of_boxes);
    boxes
}

/// Bnd_Box::IsOut(const gp_Pln&) (Bnd_Box.cxx:709-757): true when the box lies
/// strictly on one side of the plane.
fn is_out_plane(bnd_box: &BndBox, plane: &GpPln) -> bool {
    if bnd_box.is_whole() {
        return false; // Bnd_Box.cxx:711-714
    }
    if bnd_box.is_void() {
        return true; // Bnd_Box.cxx:715-718
    }
    let (a, b, c, d) = plane.coefficients(); // Bnd_Box.cxx:721-722
    // get() order is (xmin, xmax, ymin, ymax, zmin, zmax) (box3d.rs:139).
    let (xmin, xmax, ymin, ymax, zmin, zmax) = bnd_box.get().unwrap(); // Bnd_Box.cxx:723-728
    let f = |x: f64, y: f64, z: f64| a * x + b * y + c * z + d;
    let plus = f(xmin, ymin, zmin) > 0.0; // Bnd_Box.cxx:729-730
    if plus != (f(xmin, ymin, zmax) > 0.0) {
        return false; // Bnd_Box.cxx:731-734
    }
    if plus != (f(xmin, ymax, zmin) > 0.0) {
        return false; // Bnd_Box.cxx:735-738
    }
    if plus != (f(xmin, ymax, zmax) > 0.0) {
        return false; // Bnd_Box.cxx:739-742
    }
    if plus != (f(xmax, ymin, zmin) > 0.0) {
        return false; // Bnd_Box.cxx:743-746
    }
    if plus != (f(xmax, ymin, zmax) > 0.0) {
        return false; // Bnd_Box.cxx:747-750
    }
    if plus != (f(xmax, ymax, zmin) > 0.0) {
        return false; // Bnd_Box.cxx:751-754
    }
    plus == (f(xmax, ymax, zmax) > 0.0) // Bnd_Box.cxx:755
}
