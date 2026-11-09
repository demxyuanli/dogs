use super::entities::{BoundaryCurve, EntityData};

fn push_boundary(v: &mut Vec<usize>, surface: &Option<usize>, curves: &[BoundaryCurve]) {
    if let Some(s) = surface {
        v.push(*s);
    }
    for curve in curves {
        v.push(curve.model);
        v.extend(curve.parameters.iter().copied());
    }
}

pub(in crate::iges_read) fn pointer_targets(data: &EntityData) -> Vec<usize> {
    let mut v = Vec::new();
    match data {
        EntityData::Group { entities } => v.extend(entities.iter().copied()),
        EntityData::SingleParent { parent, children } => {
            if let Some(p) = parent {
                v.push(*p);
            }
            for c in children.iter().flatten() {
                v.push(*c);
            }
        }
        EntityData::SubfigureDef { entities, .. } => {
            for entity in entities {
                if let Some(index) = entity {
                    v.push(*index);
                }
            }
        }
        EntityData::SingularSubfigure { subfigure, .. } => v.push(*subfigure),
        EntityData::Trimmed {
            surface,
            outer,
            inners,
        } => {
            v.push(*surface);
            if let Some(o) = outer {
                v.push(*o);
            }
            for i in inners.iter().flatten() {
                v.push(*i);
            }
        }
        EntityData::Plane { bounding, .. } => {
            if let Some(b) = bounding {
                v.push(*b);
            }
        }
        EntityData::CurveOnSurface {
            surface,
            curve_uv,
            curve_3d,
            ..
        } => {
            if let Some(s) = surface {
                v.push(*s);
            }
            if let Some(u) = curve_uv {
                v.push(*u);
            }
            if let Some(c) = curve_3d {
                v.push(*c);
            }
        }
        EntityData::Composite { curves } => v.extend(curves.iter().copied()),
        EntityData::Tabulated { directrix, .. } => v.push(*directrix),
        EntityData::OffsetSurface { surface, .. } => v.push(*surface),
        EntityData::Bounded {
            surface, boundaries, ..
        } => {
            v.push(*surface);
            v.extend(boundaries.iter().copied());
        }
        EntityData::Boundary {
            surface, curves, ..
        } => {
            push_boundary(&mut v, surface, curves);
        }
        EntityData::RuledSurface { first, second, .. } => {
            v.push(*first);
            v.push(*second);
        }
        EntityData::Revolution { axis, generatrix, .. } => {
            v.push(*axis);
            v.push(*generatrix);
        }
        EntityData::PlaneSurface {
            location,
            normal,
            reference,
        } => {
            v.push(*location);
            v.push(*normal);
            if let Some(r) = reference {
                v.push(*r);
            }
        }
        EntityData::CylindricalSurface {
            location,
            axis,
            reference,
            ..
        }
        | EntityData::ConicalSurface {
            location,
            axis,
            reference,
            ..
        } => {
            v.push(*location);
            v.push(*axis);
            if let Some(r) = reference {
                v.push(*r);
            }
        }
        EntityData::SphericalSurface {
            center,
            axis,
            reference,
            ..
        } => {
            v.push(*center);
            if let Some(a) = axis {
                v.push(*a);
            }
            if let Some(r) = reference {
                v.push(*r);
            }
        }
        EntityData::ToroidalSurface {
            center,
            axis,
            reference,
            ..
        } => {
            v.push(*center);
            v.push(*axis);
            if let Some(r) = reference {
                v.push(*r);
            }
        }
        EntityData::Line { .. }
        | EntityData::CircularArc { .. }
        | EntityData::ConicArc { .. }
        | EntityData::Copious { .. }
        | EntityData::SplineCurve { .. }
        | EntityData::SplineSurface { .. }
        | EntityData::TransformMatrix { .. }
        | EntityData::BSplineCurve { .. }
        | EntityData::BSplineSurface { .. } => {}
        EntityData::OffsetCurve { base, function, .. } => {
            v.push(*base);
            if let Some(f) = function {
                v.push(*f);
            }
        }
        EntityData::Point { .. }
        | EntityData::Direction { .. }
        | EntityData::Raw { .. } => {}
        // `IGESSolid_*::OwnShared` (`IGESSolid_ToolVertexList.cxx:...`,
        // `EdgeList.cxx:...`, `Loop.cxx:...`, `Face.cxx:...`, `Shell.cxx:...`,
        // `ManifoldSolid.cxx:165-179`).
        EntityData::VertexList { .. } => {}
        EntityData::EdgeList { edges } => {
            for edge in edges {
                if let Some(c) = edge.curvature {
                    v.push(c);
                }
                if let Some(l) = edge.start_vertex_list {
                    v.push(l);
                }
                if let Some(l) = edge.end_vertex_list {
                    v.push(l);
                }
            }
        }
        EntityData::SolidLoop { edges } => {
            for edge in edges {
                if let Some(e) = edge.edge {
                    v.push(e);
                }
                for c in edge.param_curves.iter().flatten() {
                    v.push(*c);
                }
            }
        }
        EntityData::SolidFace { surface, loops, .. } => {
            if let Some(s) = surface {
                v.push(*s);
            }
            for l in loops.iter().flatten() {
                v.push(*l);
            }
        }
        EntityData::SolidShell { faces, .. } => {
            for f in faces.iter().flatten() {
                v.push(*f);
            }
        }
        EntityData::ManifoldSolid {
            shell, void_shells, ..
        } => {
            if let Some(s) = shell {
                v.push(*s);
            }
            for s in void_shells.iter().flatten() {
                v.push(*s);
            }
        }
    }
    v
}
