//! Phase 5 module: render_svg — SVG wireframe and shaded rendering.
//!
//! **UNPORTED (audit A14)**: OCCT does not write SVG — its visualization is
//! TKV3d (`V3d_View`, `AIS_*`, OpenGl). This module is a port-local raster/vector
//! renderer inspired by `V3d_View`'s orthographic projection, not a translation.
//!
//! Ports OCCT's `V3d_View` orthographic projection to a 2-D SVG scene. A
//! shape is tessellated (`BRepMesh_IncrementalMesh`), projected perpendicular
//! to a view direction (`HLRBRep`), painter-sorted, and emitted as either
//! filled `<polygon>`s (depth-shaded grayscale or per-shape color) or a
//! `<polyline>` wireframe of the visible edges. Back faces are culled by the
//! projection stage, so only front-facing (and silhouette) geometry is drawn.

use occt_core::gp::GpVec;
use occt_core::quantity::Color;

use crate::hlr;
use crate::hlr::ProjectedMesh;

use crate::model::BRepModel;
use crate::shape::TopoShape;

/// Options controlling an SVG render.
#[derive(Debug, Clone)]
pub struct SvgRenderOptions {
    /// SVG canvas width in user units.
    pub width: f64,
    /// SVG canvas height in user units.
    pub height: f64,
    /// View direction for the orthographic projection.
    pub view_dir: GpVec,
    /// Emit shaded filled polygons (`true`) or a wireframe (`false`).
    pub fill: bool,
    /// Background color string (used for the `background` attribute).
    pub background: String,
    /// Stroke width for polygon/polyline outlines.
    pub stroke_width: f64,
}

impl Default for SvgRenderOptions {
    fn default() -> Self {
        Self {
            width: 800.0,
            height: 600.0,
            view_dir: GpVec::new(0.0, 0.0, 1.0),
            fill: true,
            background: "#ffffff".into(),
            stroke_width: 1.0,
        }
    }
}

/// Render `shape` to an SVG document.
///
/// The shape's mesh is projected along `options.view_dir`; filled output
/// painter-sorts the front-facing triangles into `<polygon>` elements, while
/// `fill=false` emits only the visible edges as `<polyline>`s. The projected
/// geometry is scaled and centered to fit the `width` × `height` canvas.
pub fn render_svg(shape: &TopoShape, options: &SvgRenderOptions, deflection: f64) -> String {
    let mesh = crate::shape_mesh::mesh_shape(shape, deflection);
    let pm = hlr::orthographic_project(&mesh, &options.view_dir);
    if options.fill {
        render_filled(&pm, None, options)
    } else {
        render_wireframe(&pm, options)
    }
}

/// Render every shape of a model as one merged projected scene.
///
/// Each shape is projected independently; the resulting triangles are
/// concatenated into a single mesh (vertex indices rebased) and drawn in
/// painter order. When a shape carries a [`Color`], its triangles use that
/// color; otherwise they fall back to depth-based grayscale.
pub fn render_svg_scene(model: &BRepModel, options: &SvgRenderOptions, deflection: f64) -> String {
    if options.fill {
        let mut points = Vec::new();
        let mut triangles = Vec::new();
        let mut colors: Vec<Option<Color>> = Vec::new();
        let mut dmin = f64::INFINITY;
        let mut dmax = f64::NEG_INFINITY;
        for ms in &model.shapes {
            let mesh = crate::shape_mesh::mesh_shape(&ms.shape, deflection);
            let pm = hlr::orthographic_project(&mesh, &options.view_dir);
            let base = points.len();
            points.extend(pm.points.iter().copied());
            for t in &pm.triangles {
                dmin = dmin.min(t.depth);
                dmax = dmax.max(t.depth);
                triangles.push(hlr::ProjectedTriangle {
                    a: base + t.a,
                    b: base + t.b,
                    c: base + t.c,
                    depth: t.depth,
                    normal_z: t.normal_z,
                });
                colors.push(ms.color);
            }
        }
        let merged = ProjectedMesh { points, triangles };
        render_filled_scene(&merged, &colors, dmin, dmax, options)
    } else {
        // Wireframe scene: emit each shape's visible edges sequentially.
        let mut out = String::new();
        out.push_str(&header(options));
        for ms in &model.shapes {
            let mesh = crate::shape_mesh::mesh_shape(&ms.shape, deflection);
            let pm = hlr::orthographic_project(&mesh, &options.view_dir);
            out.push_str(&polyline_body(&pm, options));
        }
        out.push_str("</svg>\n");
        out
    }
}

/// Count `<polygon` elements in an SVG document (test helper).
pub fn svg_polygon_count(svg: &str) -> usize {
    svg.matches("<polygon").count()
}

// ---------------------------------------------------------------------------
// Internals
// ---------------------------------------------------------------------------

fn header(options: &SvgRenderOptions) -> String {
    format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{}\" height=\"{}\" viewBox=\"0 0 {} {}\" background=\"{}\">\n",
        options.width, options.height, options.width, options.height, options.background
    )
}

/// Painter-sorted filled rendering of a single projection.
fn render_filled(pm: &ProjectedMesh, color: Option<Color>, options: &SvgRenderOptions) -> String {
    let f = fit(&pm.points, options.width, options.height);
    let order = hlr::painter_sort(pm);
    let (dmin, dmax) = depth_range(pm);
    let mut out = String::new();
    out.push_str(&header(options));
    for &ti in &order {
        let t = &pm.triangles[ti];
        let (xa, ya) = map(&pm.points[t.a], f, options.width, options.height);
        let (xb, yb) = map(&pm.points[t.b], f, options.width, options.height);
        let (xc, yc) = map(&pm.points[t.c], f, options.width, options.height);
        if signed_area(xa, ya, xb, yb, xc, yc) <= 1e-12 {
            continue; // edge-on triangle → zero projected area
        }
        let fill = color.map(color_to_rgb).unwrap_or_else(|| grayscale(t.depth, dmin, dmax));
        out.push_str(&format!(
            "<polygon points=\"{:.3},{:.3} {:.3},{:.3} {:.3},{:.3}\" fill=\"{fill}\" stroke=\"black\" stroke-width=\"{}\" />\n",
            xa, ya, xb, yb, xc, yc, options.stroke_width
        ));
    }
    out.push_str("</svg>\n");
    out
}

/// Filled rendering of a merged multi-shape projection with per-triangle color.
fn render_filled_scene(
    pm: &ProjectedMesh,
    colors: &[Option<Color>],
    dmin: f64,
    dmax: f64,
    options: &SvgRenderOptions,
) -> String {
    let f = fit(&pm.points, options.width, options.height);
    let order = hlr::painter_sort(pm);
    let mut out = String::new();
    out.push_str(&header(options));
    for &ti in &order {
        let t = &pm.triangles[ti];
        let (xa, ya) = map(&pm.points[t.a], f, options.width, options.height);
        let (xb, yb) = map(&pm.points[t.b], f, options.width, options.height);
        let (xc, yc) = map(&pm.points[t.c], f, options.width, options.height);
        if signed_area(xa, ya, xb, yb, xc, yc) <= 1e-12 {
            continue;
        }
        let fill = colors
            .get(ti)
            .and_then(|c| *c)
            .map(color_to_rgb)
            .unwrap_or_else(|| grayscale(t.depth, dmin, dmax));
        out.push_str(&format!(
            "<polygon points=\"{:.3},{:.3} {:.3},{:.3} {:.3},{:.3}\" fill=\"{fill}\" stroke=\"black\" stroke-width=\"{}\" />\n",
            xa, ya, xb, yb, xc, yc, options.stroke_width
        ));
    }
    out.push_str("</svg>\n");
    out
}

/// Wireframe `<polyline>`s of the visible edges of one projection.
fn polyline_body(pm: &ProjectedMesh, options: &SvgRenderOptions) -> String {
    let f = fit(&pm.points, options.width, options.height);
    let edges = hlr::visible_edges(pm, 1e-9);
    let mut out = String::new();
    for ((x1, y1), (x2, y2)) in edges {
        let (mx1, my1) = map_xy(x1, y1, f, options.width, options.height);
        let (mx2, my2) = map_xy(x2, y2, f, options.width, options.height);
        out.push_str(&format!(
            "<polyline points=\"{:.3},{:.3} {:.3},{:.3}\" fill=\"none\" stroke=\"black\" stroke-width=\"{}\" />\n",
            mx1, my1, mx2, my2, options.stroke_width
        ));
    }
    out
}

/// Wireframe rendering of a single projection.
fn render_wireframe(pm: &ProjectedMesh, options: &SvgRenderOptions) -> String {
    let mut out = String::new();
    out.push_str(&header(options));
    out.push_str(&polyline_body(pm, options));
    out.push_str("</svg>\n");
    out
}

/// Projected-depth range of all triangles.
fn depth_range(pm: &ProjectedMesh) -> (f64, f64) {
    let mut dmin = f64::INFINITY;
    let mut dmax = f64::NEG_INFINITY;
    for t in &pm.triangles {
        dmin = dmin.min(t.depth);
        dmax = dmax.max(t.depth);
    }
    (dmin, dmax)
}

/// Signed doubled area of a triangle in view-plane coordinates.
fn signed_area(xa: f64, ya: f64, xb: f64, yb: f64, xc: f64, yc: f64) -> f64 {
    (xb - xa) * (yc - ya) - (yb - ya) * (xc - xa)
}

/// Fit transform: `(minx, miny, sx, sy)` mapping projected coords to the
/// padded canvas (5% margin on each side). Degenerate/empty extents fall back
/// to a unit window so the output stays finite.
fn fit(points: &[hlr::ProjectedPoint], width: f64, height: f64) -> (f64, f64, f64, f64) {
    let mut minx = f64::INFINITY;
    let mut maxx = f64::NEG_INFINITY;
    let mut miny = f64::INFINITY;
    let mut maxy = f64::NEG_INFINITY;
    for p in points {
        minx = minx.min(p.x);
        maxx = maxx.max(p.x);
        miny = miny.min(p.y);
        maxy = maxy.max(p.y);
    }
    if !minx.is_finite() || !(maxx > minx) {
        minx = -1.0;
        maxx = 1.0;
    }
    if !miny.is_finite() || !(maxy > miny) {
        miny = -1.0;
        maxy = 1.0;
    }
    let sx = width * 0.9 / (maxx - minx);
    let sy = height * 0.9 / (maxy - miny);
    (minx, miny, sx, sy)
}

/// Map a projected point to SVG canvas coordinates (Y flipped for SVG).
fn map(p: &hlr::ProjectedPoint, f: (f64, f64, f64, f64), width: f64, height: f64) -> (f64, f64) {
    map_xy(p.x, p.y, f, width, height)
}

fn map_xy(x: f64, y: f64, f: (f64, f64, f64, f64), width: f64, height: f64) -> (f64, f64) {
    let (minx, miny, sx, sy) = f;
    let mx = width * 0.05 + (x - minx) * sx;
    let my = height * 0.05 + (y - miny) * sy;
    (mx, height - my)
}

/// `Color` → `rgb(r,g,b)` SVG fill string.
fn color_to_rgb(c: Color) -> String {
    format!(
        "rgb({},{},{})",
        (c.r * 255.0).round() as i32,
        (c.g * 255.0).round() as i32,
        (c.b * 255.0).round() as i32
    )
}

/// Depth-based grayscale fill: nearer triangles are lighter, farther darker.
fn grayscale(depth: f64, dmin: f64, dmax: f64) -> String {
    let t = if dmax > dmin { (depth - dmin) / (dmax - dmin) } else { 0.5 };
    let g = (200.0 - 150.0 * t.clamp(0.0, 1.0)).round() as i32;
    format!("rgb({g},{g},{g})")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::BRepPrimBox;
    use crate::shape::TopoShape;

    fn unit_box() -> TopoShape {
        BRepPrimBox::make_box(1.0, 1.0, 1.0).solid.0
    }

    /// Every `points="..."` attribute in the document parses to finite
    /// numbers within the canvas.
    fn assert_points_within_viewbox(svg: &str, w: f64, h: f64) {
        for attr in svg.split("<polygon points=\"") {
            let Some(end) = attr.find('"') else { continue };
            let body = &attr[..end];
            let mut nums: Vec<f64> = Vec::new();
            for part in body.split(|c: char| c == ',' || c.is_whitespace()) {
                if let Ok(v) = part.parse::<f64>() {
                    nums.push(v);
                }
            }
            assert!(
                nums.iter().all(|v| v.is_finite()),
                "non-finite coordinate in {body}"
            );
            for (i, v) in nums.iter().enumerate() {
                let limit = if i % 2 == 0 { w } else { h };
                assert!(
                    (0.0..=limit).contains(v),
                    "coordinate {v} outside viewBox (limit {limit}) in {body}"
                );
            }
        }
    }

    #[test]
    fn box_filled_render_front_faces_only() {
        let svg = render_svg(&unit_box(), &SvgRenderOptions::default(), 0.25);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        // Only front-facing triangles are drawn: exactly the triangles of the
        // near face (edge-on side faces project to zero area and are skipped).
        let mesh = crate::shape_mesh::mesh_shape(&unit_box(), 0.25);
        let pm = hlr::orthographic_project(&mesh, &GpVec::new(0.0, 0.0, 1.0));
        let front = pm.triangles.iter().filter(|t| t.normal_z < 0.0).count();
        assert!(front > 0, "box has front-facing triangles");
        assert_eq!(svg_polygon_count(&svg), front, "polygons:\n{svg}");
        assert_points_within_viewbox(&svg, 800.0, 600.0);
    }

    #[test]
    fn box_wireframe_emits_polyline() {
        let opts = SvgRenderOptions { fill: false, ..Default::default() };
        let svg = render_svg(&unit_box(), &opts, 0.25);
        assert!(svg.contains("<polyline"));
        assert!(!svg.contains("<polygon"));
    }

    #[test]
    fn scene_merges_two_shapes() {
        let mut model = BRepModel::new();
        model.add("A", unit_box());
        model.add("B", unit_box());
        let svg = render_svg_scene(&model, &SvgRenderOptions::default(), 0.25);
        let mesh = crate::shape_mesh::mesh_shape(&unit_box(), 0.25);
        let pm = hlr::orthographic_project(&mesh, &GpVec::new(0.0, 0.0, 1.0));
        let front = pm.triangles.iter().filter(|t| t.normal_z < 0.0).count();
        // Two boxes, each contributing its front-facing triangles.
        assert_eq!(svg_polygon_count(&svg), 2 * front, "polygons:\n{svg}");
        assert_points_within_viewbox(&svg, 800.0, 600.0);
    }

    #[test]
    fn scene_uses_shape_color() {
        let mut model = BRepModel::new();
        model.add_with_color("Red", unit_box(), Color::RED);
        let svg = render_svg_scene(&model, &SvgRenderOptions::default(), 0.25);
        assert!(
            svg.contains("fill=\"rgb(255,0,0)\""),
            "expected red fill, got:\n{svg}"
        );
    }

    #[test]
    fn svg_polygon_count_counts_elements() {
        let svg = "<svg></svg><polygon a/><polygon b/>";
        assert_eq!(svg_polygon_count(svg), 2);
    }

    #[test]
    fn empty_model_yields_valid_svg() {
        let svg = render_svg_scene(&BRepModel::new(), &SvgRenderOptions::default(), 0.25);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert_eq!(svg_polygon_count(&svg), 0);
    }
}
