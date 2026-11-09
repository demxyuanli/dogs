use super::prelude::*;
use super::statics::*;


/// Axis-aligned box primitive. Source: `BRepPrimAPI_MakeBox`.
#[derive(Debug, Clone)]
pub struct BRepPrimBox {
    pub solid: Solid,
    pub bbox: BndBox,
}


impl BRepPrimBox {
    /// Axis-aligned box from the origin with dimensions `dx × dy × dz`.
    /// Builds the full boundary: 8 vertices, 12 line edges, 6 planar faces.
    pub fn make_box(dx: f64, dy: f64, dz: f64) -> Self {
        assert!(
            dx > 0.0 && dy > 0.0 && dz > 0.0,
            "BRepPrimBox::make_box: dimensions must be positive"
        );
        let b = TopoBuilder::new();
        let corners = box_corners(dx, dy, dz);
        let solid = build_box(&b, &corners);
        let bbox = BndBox::from_corners(&GpPnt::zero(), &GpPnt::new(dx, dy, dz));
        Self { solid, bbox }
    }

    /// Box spanning two corner points (axis-aligned). Same `BRepPrim_GWedge`
    /// shared-TShape construction as [`Self::make_box`], translated to `p1`/`p2`.
    pub fn make_box_corner(p1: &GpPnt, p2: &GpPnt) -> Self {
        let min = GpPnt::new(p1.x().min(p2.x()), p1.y().min(p2.y()), p1.z().min(p2.z()));
        let max = GpPnt::new(p1.x().max(p2.x()), p1.y().max(p2.y()), p1.z().max(p2.z()));
        let dx = max.x() - min.x();
        let dy = max.y() - min.y();
        let dz = max.z() - min.z();
        assert!(
            dx > 0.0 && dy > 0.0 && dz > 0.0,
            "make_box_corner: box must have positive extent"
        );
        let b = TopoBuilder::new();
        let corners = [
            GpPnt::new(min.x(), min.y(), min.z()),
            GpPnt::new(max.x(), min.y(), min.z()),
            GpPnt::new(max.x(), max.y(), min.z()),
            GpPnt::new(min.x(), max.y(), min.z()),
            GpPnt::new(min.x(), min.y(), max.z()),
            GpPnt::new(max.x(), min.y(), max.z()),
            GpPnt::new(max.x(), max.y(), max.z()),
            GpPnt::new(min.x(), max.y(), max.z()),
        ];
        let solid = build_box(&b, &corners);
        Self {
            solid,
            bbox: BndBox::from_corners(&min, &max),
        }
    }

    pub fn volume(&self) -> f64 {
        let (x0, x1, y0, y1, z0, z1) = self.bbox.get().unwrap();
        (x1 - x0) * (y1 - y0) * (z1 - z0)
    }

    pub fn surface_area(&self) -> f64 {
        let (x0, x1, y0, y1, z0, z1) = self.bbox.get().unwrap();
        let dx = x1 - x0;
        let dy = y1 - y0;
        let dz = z1 - z0;
        2.0 * (dx * dy + dy * dz + dx * dz)
    }

    pub fn center(&self) -> GpPnt {
        let (x0, x1, y0, y1, z0, z1) = self.bbox.get().unwrap();
        GpPnt::new((x0 + x1) / 2.0, (y0 + y1) / 2.0, (z0 + z1) / 2.0)
    }

    /// Number of distinct vertices/edges/faces in the real boundary tree
    /// (uniquified by `TShape` identity, matching OCCT's `TopExp_Explorer` +
    /// `TopTools_MapOfShape`).
    pub fn counts(&self) -> (usize, usize, usize) {
        fn unique(shape: &crate::shape::TopoShape, t: ShapeType) -> usize {
            let mut seen = std::collections::HashSet::new();
            let mut ex = Explorer::new(shape, t);
            while ex.more() {
                seen.insert(Arc::as_ptr(&ex.current().tshape));
                ex.next();
            }
            seen.len()
        }
        (unique(&self.solid.0, ShapeType::Vertex), unique(&self.solid.0, ShapeType::Edge), unique(&self.solid.0, ShapeType::Face))
    }
}
