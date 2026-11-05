//! `ChFi3d_IsInFront` — do two surf data face each other across a common face?
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:906-1216`
//! (OCCT 8.0.0).
//!
//! The body is the same block four times, once per (side of `cd1`, side of
//! `cd2`) pair that shares a face index: the face is re-read from the DS, the
//! two stripe orientations are compared with `ChFi3d::SameSide`, the 2d
//! intersection of the two p-curves is taken (`ChFi3d_IntTraces`) and the
//! answer of the first arm that produced one is kept. Arms two to four apply
//! OCCT's `restore` rule (`cxx:1032-1054` / `:1105-1128` / `:1183-1206`): an
//! intersection found later must not move either parameter backwards along its
//! own `sens`, otherwise the earlier answer is put back.
//!
//! Translation notes:
//! - `p1` / `p2` are `double&` out-parameters in OCCT and are **in/out**: a
//!   failing `ChFi3d_IntTraces` leaves them as they were. [`IsInFrontState`]
//!   carries them in and out so that the caller's value survives a failure,
//! - `face` / `sameside` / `jf1` / `jf2` / `visavis` / `ok` are the other
//!   out-parameters; `ok` is the `bool` the OCCT function returns,
//! - `Standard_NullObject` (`cxx:944`, `:952`, `:996`, `:1004`, `:1072`,
//!   `:1080`, `:1148`, `:1156`) is the `Err` of the port. OCCT has already
//!   written `jf1` / `jf2` / `face` when it raises; the port answers `Err`
//!   before the state is read, like every other `chfi3d` helper,
//! - `gp_Pnt2d P2d;` (`cxx:961` etc.) is only filled when `Check2dDistance`;
//!   `BRep_Tool::Parameters(Vtx, face)` (`cxx:965`) raises
//!   `Standard_NoSuchObject` (`BRep_Tool.cxx:1702`) and becomes the `Err` of
//!   [`BRepTool::parameters_of_vertex_on_face`],
//! - `ChFi3d_IntTraces` is [`crate::chfi3d::builder_0_helpers::ch_fi3d_int_traces`];
//!   it answers `Option<(p1, p2)>`.

use occt_core::gp::GpPnt2d;

use crate::brep_tool::BRepTool;
use crate::chfi3d::side::same_side;
use crate::chfi_ds::{Stripe, SurfData};
use crate::shape::{Face, TopoShape};
use crate::topopebrep::ds::HDataStructure;

use super::ch_fi3d_int_traces;

/// The locals and out-parameters of `ChFi3d_IsInFront` that outlive one arm.
///
/// Field names follow the OCCT locals: `p1` / `p2` / `face` / `sameside` /
/// `jf1` / `jf2` / `visavis` are the out-parameters, `u1` / `u2` / `ss` / `j1`
/// / `j2` / `ff` are the saved answer of the previous arm and `ok` is the
/// returned flag.
#[derive(Debug, Clone)]
pub struct IsInFrontState {
    /// `p1`
    pub p1: f64,
    /// `p2`
    pub p2: f64,
    /// `face`
    pub face: Option<Face>,
    /// `sameside`
    pub sameside: bool,
    /// `jf1`
    pub jf1: i32,
    /// `jf2`
    pub jf2: i32,
    /// `visavis`
    pub visavis: bool,
    /// `ok`
    pub ok: bool,
    /// `u1`
    pub u1: f64,
    /// `u2`
    pub u2: f64,
    /// `ss`
    pub ss: bool,
    /// `j1`
    pub j1: i32,
    /// `j2`
    pub j2: i32,
    /// `ff`
    pub ff: Option<Face>,
}

impl Default for IsInFrontState {
    fn default() -> Self {
        Self {
            p1: 0.0,
            p2: 0.0,
            face: None,
            sameside: false,
            jf1: 0,
            jf2: 0,
            visavis: false,
            ok: false,
            u1: 0.0,
            u2: 0.0,
            ss: false,
            j1: 0,
            j2: 0,
            ff: None,
        }
    }
}

impl IsInFrontState {
    /// The state a caller enters with: `p1` / `p2` as OCCT passes them in,
    /// every other field zeroed (`cxx:932-937`).
    pub fn new(p1: f64, p2: f64) -> Self {
        Self {
            p1,
            p2,
            ..Self::default()
        }
    }
}

/// One of the four `if (fd1->IndexOfSx() == fd2->IndexOfSy())` arms
/// (`cxx:989-990`, `:1063-1064`, `:1141-1142` and the `:936-937` one, which is
/// the only arm without the `restore` rule).
///
/// `restore` is `false` for the `S1`/`S1` arm.
#[allow(clippy::too_many_arguments)]
fn arm(
    ds: &HDataStructure,
    cd1: &Stripe,
    cd2: &Stripe,
    fd1: &SurfData,
    fd2: &SurfData,
    jf1: i32,
    jf2: i32,
    isf1: bool,
    isf2: bool,
    sens1: i32,
    sens2: i32,
    vtx: Option<&TopoShape>,
    check_2d_distance: bool,
    enlarge: bool,
    restore: bool,
    st: &mut IsInFrontState,
) -> Result<(), String> {
    // `face = TopoDS::Face(DStr.Shape(fd1->Index(jf1)));` + the null test.
    let face = ds
        .shape(fd1.index(jf1), false)
        .ok_or_else(|| {
            "ChFi3d_IsInFront : Trying to check orientation of NULL face".to_string()
        })?;
    let face = Face::wrap(face).ok_or_else(|| {
        "ChFi3d_IsInFront : Trying to check orientation of NULL face".to_string()
    })?;
    st.jf1 = jf1;
    st.jf2 = jf2;
    st.face = Some(face.clone());

    // `OrSave1 = cd1->Orientation(jf1); Or = OrFace1 = face.Orientation();
    //  OrSave2 = cd2->Orientation(jf2);`
    let or_save1 = cd1.orientation_on_surface(jf1);
    let or = face.0.orientation();
    let or_face1 = or;
    let or_save2 = cd2.orientation_on_surface(jf2);

    // `const TopoDS_Shape& shape2 = DStr.Shape(fd2->Index(jf2));` + null test.
    let shape2 = ds
        .shape(fd2.index(jf2), false)
        .ok_or_else(|| {
            "ChFi3d_IsInFront : Trying to check orientation of NULL shape".to_string()
        })?;
    let or_face2 = shape2.orientation();

    st.visavis = true;
    st.sameside = same_side(or, or_save1, or_save2, or_face1, or_face2);

    // The parameters of the other side are not used for orientation.
    let pref1 = fd1.interference(jf1).parameter(isf1);
    let pref2 = fd2.interference(jf2).parameter(isf2);

    // `if (Check2dDistance) P2d = BRep_Tool::Parameters(Vtx, face);`
    let p2d = if check_2d_distance {
        let vtx = vtx.ok_or_else(|| {
            "ChFi3d_IsInFront : BRep_Tool::Parameters with a null vertex".to_string()
        })?;
        BRepTool::parameters_of_vertex_on_face(vtx, &face)?
    } else {
        GpPnt2d::new(0.0, 0.0)
    };

    let Some((new_p1, new_p2)) = ch_fi3d_int_traces(
        fd1,
        pref1,
        jf1,
        sens1,
        fd2,
        pref2,
        jf2,
        sens2,
        &p2d,
        check_2d_distance,
        enlarge,
    ) else {
        // `else if (ok)` of the three `restore` arms: the re-initialization
        // added in case `p1`, ... take wrong values.
        if restore && st.ok {
            st.p1 = st.u1;
            st.p2 = st.u2;
            st.sameside = st.ss;
            st.jf1 = st.j1;
            st.jf2 = st.j2;
            st.face = st.ff.clone();
        }
        return Ok(());
    };

    st.p1 = new_p1;
    st.p2 = new_p2;

    if restore {
        // `bool restore = ok && ((j1 == jf1 && sens1 * (p1 - u1) > 0.) ||
        //                       (j2 == jf2 && sens2 * (p2 - u2) > 0.));`
        let restore_now = st.ok
            && ((st.j1 == jf1 && f64::from(sens1) * (st.p1 - st.u1) > 0.0)
                || (st.j2 == jf2 && f64::from(sens2) * (st.p2 - st.u2) > 0.0));
        st.ok = true;
        if restore_now {
            st.p1 = st.u1;
            st.p2 = st.u2;
            st.sameside = st.ss;
            st.jf1 = st.j1;
            st.jf2 = st.j2;
            st.face = st.ff.clone();
        } else {
            st.u1 = st.p1;
            st.u2 = st.p2;
            st.ss = st.sameside;
            st.j1 = st.jf1;
            st.j2 = st.jf2;
            st.ff = st.face.clone();
        }
    } else {
        // `cxx:980-985`: the first arm only records, it never restores.
        st.u1 = st.p1;
        st.u2 = st.p2;
        st.ss = st.sameside;
        st.j1 = st.jf1;
        st.j2 = st.jf2;
        st.ff = st.face.clone();
        st.ok = true;
    }
    Ok(())
}

/// `ChFi3d_IsInFront(DStr, cd1, cd2, i1, i2, sens1, sens2, p1, p2, face,
/// sameside, jf1, jf2, visavis, Vtx, Check2dDistance, enlarge)`
/// (`cxx:909-1217`).
///
/// `i1` / `i2` are 1-based indices into the two stripes' surf data.
#[allow(clippy::too_many_arguments)]
pub fn ch_fi3d_is_in_front(
    ds: &HDataStructure,
    cd1: &Stripe,
    cd2: &Stripe,
    i1: i32,
    i2: i32,
    sens1: i32,
    sens2: i32,
    st: &mut IsInFrontState,
    vtx: Option<&TopoShape>,
    check_2d_distance: bool,
    enlarge: bool,
) -> Result<(), String> {
    let isf1 = sens1 == 1;
    let isf2 = sens2 == 1;

    // `const occ::handle<ChFiDS_SurfData>& fd1 = cd1->SetOfSurfData()->Value(i1);`
    let fd1 = surf_data(cd1, i1)?;
    let fd2 = surf_data(cd2, i2)?;

    // `visavis = false;` (`cxx:934`).
    st.visavis = false;

    // `cxx:936-987`.
    if fd1.index_of_s1() == fd2.index_of_s1() {
        arm(
            ds,
            cd1,
            cd2,
            &fd1,
            &fd2,
            1,
            1,
            isf1,
            isf2,
            sens1,
            sens2,
            vtx,
            check_2d_distance,
            enlarge,
            false,
            st,
        )?;
    }
    // `cxx:989-1061`.
    if fd1.index_of_s2() == fd2.index_of_s1() {
        arm(
            ds,
            cd1,
            cd2,
            &fd1,
            &fd2,
            2,
            1,
            isf1,
            isf2,
            sens1,
            sens2,
            vtx,
            check_2d_distance,
            enlarge,
            true,
            st,
        )?;
    }
    // `cxx:1063-1140`.
    if fd1.index_of_s1() == fd2.index_of_s2() {
        arm(
            ds,
            cd1,
            cd2,
            &fd1,
            &fd2,
            1,
            2,
            isf1,
            isf2,
            sens1,
            sens2,
            vtx,
            check_2d_distance,
            enlarge,
            true,
            st,
        )?;
    }
    // `cxx:1141-1215`.
    if fd1.index_of_s2() == fd2.index_of_s2() {
        arm(
            ds,
            cd1,
            cd2,
            &fd1,
            &fd2,
            2,
            2,
            isf1,
            isf2,
            sens1,
            sens2,
            vtx,
            check_2d_distance,
            enlarge,
            true,
            st,
        )?;
    }
    Ok(())
}

/// `cd->SetOfSurfData()->Value(i)`, `i` 1-based.
pub(crate) fn surf_data(cd: &Stripe, i: i32) -> Result<std::sync::Arc<SurfData>, String> {
    if i < 1 {
        return Err(format!("ChFi3d: surf data index {i} is out of range"));
    }
    cd.set_of_surf_data()
        .get(i as usize - 1)
        .cloned()
        .ok_or_else(|| format!("ChFi3d: surf data index {i} is out of range"))
}
