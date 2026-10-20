//! `ShapeFix_Face::FixPeriodicDegenerated` (`ShapeFix_Face.cxx:3101-3259`).

use super::prelude::*;
use super::*;

impl ShapeFixFace {
    /// `ShapeFix_Face::FixPeriodicDegenerated` (`ShapeFix_Face.cxx:3101-3259`):
    /// a conical face whose single wire belts the surface by one period gets a
    /// degenerated apex edge and a second (apex) wire, so that the following
    /// `FixMissingSeam` sees the 2-wire configuration and stitches the seam.
    /// Without it a cone read from STEP keeps 1 edge and 0 V-span, and the
    /// face's discrete range is empty (`IMeshData_Failure`).
    ///
    /// `Precision()` in this routine is `ShapeFix_Root::Precision()`
    /// (`ShapeFix_Root.lxx:34`) = `Precision::Confusion()` by default
    /// (`ShapeFix_Root.cxx:26`; the STEP operator does not `SetPrecision`).
    pub fn fix_periodic_degenerated(&mut self) -> bool {
        // cxx:3105-3112.
        let mut face = match &self.face {
            Some(f) => f.clone(),
            None => return false,
        };
        {
            let applied = self.context.apply(&face.0);
            if applied.is_face() {
                face = Face(applied);
                self.face = Some(face.clone());
            }
        }

        // cxx:3116-3129: `TopoDS_Iterator(myFace, false)` — oriented wires only.
        let children: Vec<TopoShape> = {
            let ts = face.0.tshape.read().unwrap();
            ts.children.clone()
        };
        let mut wire_seq: Vec<Wire> = Vec::new();
        for child in children {
            let o = child.orientation();
            if child.shape_type() != ShapeType::Wire
                || (o != Orientation::Forward && o != Orientation::Reversed)
            {
                continue;
            }
            wire_seq.push(Wire(child));
        }

        // cxx:3131-3139: only a single wire on a conical surface is checked.
        let surf = match BRepTool::face_surface(&face) {
            Some(s) => s,
            None => return false,
        };
        let cone = match surf.gp_cone() {
            Some(c) => c,
            None => return false,
        };
        if wire_seq.len() != 1 {
            return false;
        }
        let mut sole_wire = wire_seq.remove(0);

        // cxx:3146-3161: does the wire belt the cone by one period?
        let lp = match is_periodic_conical_loop(&sole_wire, &face, CONFUSION) {
            Some(v) => v,
            None => return false,
        };

        // cxx:3167-3180: the base circle the cone was built from
        // (`VIso(0.0)` -> `Geom_Circle`, whose radius is `Radius`).
        let base_r = match surf.v_iso_curve(0.0).and_then(|c| c.circle_radius()) {
            Some(r) => r,
            None => return false,
        };
        let semi_angle = cone.semi_angle();
        if semi_angle.abs() <= CONFUSION {
            return false; // cxx:3177-3180: bad surface
        }

        // cxx:3183-3184: the V parameter of the apex.
        let apex_v = -(base_r / semi_angle.sin());

        // cxx:3187: `BRepBuilderAPI_MakeVertex(aConeSurf->Apex())`.
        let builder = TopoBuilder::new();
        let apex_pnt = cone.apex();

        // cxx:3198-3201: reject when the apex V sits on or between the wire's
        // V bounds — the 2D support line would not be consistent with the wire.
        if (apex_v - lp.min_v).abs() <= CONFUSION
            || (apex_v - lp.max_v).abs() <= CONFUSION
            || (apex_v < lp.max_v && apex_v > lp.min_v)
        {
            return false;
        }

        // cxx:3205-3223: the 2D apex support line plus the wire flip that puts
        // the apex before the wire along the seam direction. Exactly one of the
        // two branches holds, the third V relationship was rejected above.
        let (line_pnt, line_dir, flip_wire) = if apex_v < lp.min_v {
            (
                GpPnt2d::new(lp.min_u, apex_v),
                GpDir2d::new(1.0, 0.0).expect("dir"),
                !lp.is_u_decrease,
            )
        } else if apex_v > lp.max_v {
            (
                GpPnt2d::new(lp.max_u, apex_v),
                GpDir2d::new(-1.0, 0.0).expect("dir"),
                lp.is_u_decrease,
            )
        } else {
            return false;
        };

        // cxx:3226-3232: `UpdateEdge` (pcurve on `myFace` with the
        // `Precision()` tolerance), the degenerate vertex twice, the
        // `Degenerated` flag and `Range(E, First, Last)` (the edge's own range,
        // which `BRep_TEdge` shares between the 3D curve and the pcurves).
        let reg = GeometryRegistry::global();
        let face_key = GeometryRegistry::shape_key(&face.0);
        let line: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(line_pnt, line_dir));
        let mut apex_edge = builder.make_shape(ShapeType::Edge);
        reg.set_edge_pcurve(&apex_edge, face_key, line);
        reg.set_edge_tolerance(&apex_edge, CONFUSION);
        reg.set_degenerated(&apex_edge, true);
        // cxx:3232 `Range(E, First, Last)`: `BRep_TEdge` carries a single range
        // shared by the 3D curve and every CurveOnSurface representation, so
        // the apex edge's pcurve range becomes `[0, dU]` as well.
        reg.set_edge_range(&apex_edge, 0.0, (lp.max_u - lp.min_u).abs());
        reg.set_pcurve_range(&apex_edge, face_key, 0.0, (lp.max_u - lp.min_u).abs());
        let mut apex = builder.make_vertex(apex_pnt, CONFUSION);
        apex.0.set_orientation(Orientation::Forward);
        builder.add(&mut apex_edge, &apex.0);
        apex.0.set_orientation(Orientation::Reversed);
        builder.add(&mut apex_edge, &apex.0);
        let apex_wire = builder.make_wire(&[Edge(apex_edge)]);

        // cxx:3239-3253: the new face over the old surface, carrying the old
        // face's location and orientation (`EmptyCopied`) with the two wires.
        if flip_wire {
            sole_wire.0.reverse();
        }
        let mut new_face = builder.make_face(surf.clone(), &[sole_wire, apex_wire]);
        new_face.0.set_location(face.0.location());
        new_face.0.set_orientation(face.0.orientation());

        // cxx:3255-3257.
        self.result = Some(new_face.0.clone());
        self.context.replace(&face.0, &new_face.0);
        true
    }
}
