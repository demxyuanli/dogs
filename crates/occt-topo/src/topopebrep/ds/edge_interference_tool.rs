//! `TopOpeBRepDS_EdgeInterferenceTool` — complex transition on an edge.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EdgeInterferenceTool.hxx:31-55`,
//! `TopOpeBRepDS_EdgeInterferenceTool.cxx:30-174` (OCCT 8.0.0).
//!
//! OCCT keeps the edge's own orientation in `myEdgeOrientation` and the support
//! index written into the resulting interference in `myEdgeOriented`; both are
//! written by `Init` (`cxx:54-73`) and read by `Add` / `Transition`. An
//! `INTERNAL` / `EXTERNAL` edge never reaches `TopTrans_CurveTransition`
//! (`cxx:60-65`, `:86-90`, `:136-141`) and `Transition` answers fixed states
//! for it (`cxx:164-169`).
//!
//! Translation notes:
//! - `::Parameter(I)` (`cxx:35-51`) throws `Standard_ProgramError` for any
//!   interference that is neither an `EdgeVertexInterference` nor a
//!   `CurvePointInterference`; that is the `Err` arm here.
//! - `TopoDS_Iterator(E, false)` (`cxx:98`) is the raw child list of the edge
//!   (`cumOri = false`), i.e. `TShape::children` as stored.
//! - `TopOpeBRepTool_ShapeTool::EdgeData(E, par, T, N, C)` is
//!   [`edge_data`] applied to `BRepAdaptor_Curve(E)`; the returned `tol` is
//!   `Precision::Angular()` (`ShapeTool.cxx:625`) and `N` is only defined when
//!   `|C|` exceeds the (always `RealLast()`) `tolm` of `cxx:626-630`.

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};

use super::interference::Interference;
use crate::top_trans_curve_transition::CurveTransition;
use crate::topopebrep::tool::adaptor::CurveAdaptor;
use crate::topopebrep::tool::shape_tool::{edge_data, EdgeDataResult};
use occt_core::gp::GpDir;

/// `static double Parameter(I)` (`cxx:35-51`).
pub fn parameter_of(i: &Interference) -> Result<f64, String> {
    match i.data() {
        super::interference::InterferenceData::EdgeVertex(_)
        | super::interference::InterferenceData::CurvePoint { .. } => i
            .parameter()
            .ok_or_else(|| "TopOpeBRepDS_EdgeInterferenceTool1".to_string()),
        _ => Err("TopOpeBRepDS_EdgeInterferenceTool1".to_string()),
    }
}

/// The `gp_Dir N` a caller of `EdgeData` passes in (`cxx:71`, `:119`, `:149`).
///
/// `gp_Dir()` is `(1, 0, 0)` (`gp_Dir.hxx:52-55`) and `EdgeData` only writes
/// `N` when `|C| > tolm` with `tolm = max(Angular, max(Epsilon(0.), RealLast()))`
/// (`ShapeTool.cxx:626-632`), i.e. never. `N` is therefore always `(1, 0, 0)`
/// here, exactly as in OCCT 8.0.0 (the `xpu150399` workaround).
fn normal_of(_r: &EdgeDataResult) -> GpDir {
    GpDir::new(1.0, 0.0, 0.0).unwrap_or_else(|_| GpDir::default_dir())
}

/// `TopOpeBRepDS_EdgeInterferenceTool` (`hxx:31-55`).
#[derive(Debug, Clone)]
pub struct EdgeInterferenceTool {
    /// `myEdgeOrientation`
    edge_orientation: Orientation,
    /// `myEdgeOriented`
    edge_oriented: i32,
    /// `myTool`
    tool: CurveTransition,
}

impl EdgeInterferenceTool {
    /// `TopOpeBRepDS_EdgeInterferenceTool()` (`cxx:32`): `= default`, so
    /// `myEdgeOrientation` is `FORWARD` (the enum's 0 value), `myEdgeOriented`
    /// is 0 and `myTool` is an empty curve transition.
    pub fn new() -> Self {
        Self {
            edge_orientation: Orientation::Forward,
            edge_oriented: 0,
            tool: CurveTransition::new(),
        }
    }

    /// `Init(E, I)` (`cxx:53-74`).
    pub fn init(&mut self, e: &TopoShape, i: &Interference) -> Result<(), String> {
        self.edge_orientation = e.orientation();
        self.edge_oriented = i.support();

        // NYI 971219 : the orientation of the mother edge is not taken into
        // account.
        if self.edge_orientation == Orientation::Internal
            || self.edge_orientation == Orientation::External
        {
            return Ok(());
        }

        // par = parametre sur l'arete croisee
        let par = parameter_of(i)?;
        let r = edge_data(&CurveAdaptor::of_edge(&Edge(e.clone())), par)?;
        let n = normal_of(&r);
        self.tool.reset(&r.tangent, &n, r.curvature);
        Ok(())
    }

    /// `Add(E, V, I)` (`cxx:77-124`: the `TopoDS_Shape& V` overload).
    ///
    /// `VofE` walks `TopoDS_Iterator(E, false)`; the last visited child's
    /// orientation is the `oritan` handed to `TopTrans_CurveTransition::Compare`
    /// (`cxx:100-107`, `:122`), which is the child that matched `V` because the
    /// loop `break`s there.
    pub fn add_vertex(&mut self, e: &TopoShape, v: &TopoShape, i: &Interference) -> Result<(), String> {
        let eori = e.orientation();
        if eori == Orientation::Internal || eori == Orientation::External {
            return Ok(());
        }

        // premiere interference sur arete orientee : Init
        if self.edge_orientation == Orientation::Internal
            || self.edge_orientation == Orientation::External
        {
            return self.init(e, i);
        }

        // V est un sommet de E ?
        let children: Vec<TopoShape> = e
            .tshape
            .read()
            .expect("poisoned TShape lock")
            .children
            .clone();
        let mut v_of_e = false;
        let mut oritan = Orientation::Forward;
        for s in &children {
            if crate::topo_tools::is_same(s, v) {
                v_of_e = true;
                oritan = s.orientation();
                break;
            }
        }

        if !v_of_e {
            return Ok(());
        }
        // V est un sommet de E
        let par = BRepTool::parameter_of_vertex(e, v);
        let r = edge_data(&CurveAdaptor::of_edge(&Edge(e.clone())), par)?;
        let n = normal_of(&r);
        let oriloc = i
            .transition()
            .orientation(crate::abs::State::In, crate::abs::ShapeType::Face);
        self.tool
            .compare(r.tolerance, &r.tangent, &n, r.curvature, oriloc, oritan);
        Ok(())
    }

    /// `Add(E, P, I)` (`cxx:127-155`): the `TopOpeBRepDS_Point` overload, whose
    /// point argument OCCT ignores (`cxx:130`).
    pub fn add_point(&mut self, e: &TopoShape, i: &Interference) -> Result<(), String> {
        let eori = e.orientation();
        if eori == Orientation::Internal || eori == Orientation::External {
            return Ok(());
        }

        // premiere interference sur arete orientee : Init
        if self.edge_orientation == Orientation::Internal
            || self.edge_orientation == Orientation::External
        {
            return self.init(e, i);
        }

        let par = parameter_of(i)?;
        let r = edge_data(&CurveAdaptor::of_edge(&Edge(e.clone())), par)?;
        let n = normal_of(&r);
        let oriloc = i
            .transition()
            .orientation(crate::abs::State::In, crate::abs::ShapeType::Face);
        let oritan = Orientation::Internal;
        self.tool
            .compare(r.tolerance, &r.tangent, &n, r.curvature, oriloc, oritan);
        Ok(())
    }

    /// `Transition(I)` (`cxx:158-174`).
    pub fn transition(&self, i: &mut Interference) {
        if self.edge_orientation == Orientation::Internal {
            let t = i.change_transition();
            t.set_state_before(crate::abs::State::In);
            t.set_state_after(crate::abs::State::In);
        } else if self.edge_orientation == Orientation::External {
            let t = i.change_transition();
            t.set_state_before(crate::abs::State::Out);
            t.set_state_after(crate::abs::State::Out);
        } else {
            i.set_support(self.edge_oriented);
            let stb = self.tool.state_before();
            let sta = self.tool.state_after();
            let t = i.change_transition();
            t.set_state_before(stb);
            t.set_state_after(sta);
        }
    }
}

