use super::prelude::*;
use super::statics::*;


/// `BRepFill::Face` (`BRepFill.cxx:211-422`).
pub fn face(edge1: &Edge, edge2: &Edge) -> Result<Face, String> {
    let (same_loc, location) = shared_location(edge1, edge2);
    let (c1, v1f, v1l) = prepare_curve(edge1, same_loc)?;
    let (c2, v2f, v2l) = prepare_curve(edge2, same_loc)?;
    let closed = is_same(&v1f.0, &v1l.0) && is_same(&v2f.0, &v2l.0);
    let surf = occt_geom::geom_fill::generator_surface(&[c1, c2])?;
    let (u0, u1) = surf.u_range();
    let (v0, v1) = surf.v_range();
    let mut face = TopoBuilder::new().make_face(surf.clone(), &[]);
    GeometryRegistry::global().set_face_tolerance(&face.0, CONFUSION);

    let mut edge3 = iso_edge(surf.as_ref(), u0, v0, v1, &v1f, &v2f)?;
    let edge4 = if closed {
        edge3.clone()
    } else {
        iso_edge(surf.as_ref(), u1, v0, v1, &v1l, &v2l)?
    };
    edge3.0.reverse();

    let mut wire = Wire::new();
    let builder = TopoBuilder::new();
    builder.add_edge(&mut wire, edge1);
    builder.add_edge(&mut wire, &edge4);
    let mut edge2_rev = edge2.clone();
    edge2_rev.0.reverse();
    builder.add_edge(&mut wire, &edge2_rev);
    builder.add_edge(&mut wire, &edge3);
    wire.0.set_closed(true);
    builder.add_wire(&mut face, &wire);

    let reg = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(&face.0);
    bind_profile_pcurve(&reg, edge1, face_key, 0.0, v0, u0, u1);
    bind_profile_pcurve(&reg, edge2, face_key, 0.0, v1, u0, u1);
    if closed {
        reg.set_edge_pcurves(
            &edge3.0,
            face_key,
            vec![uv_line(u1, 0.0, 0.0, 1.0), uv_line(u0, 0.0, 0.0, 1.0)],
        );
    } else {
        reg.set_edge_pcurve(&edge3.0, face_key, uv_line(u0, 0.0, 0.0, 1.0));
        reg.set_edge_pcurve(&edge4.0, face_key, uv_line(u1, 0.0, 0.0, 1.0));
    }
    clear_same_flags(edge1);
    clear_same_flags(edge2);
    clear_same_flags(&edge3);
    clear_same_flags(&edge4);
    same_parameter_shape(&face.0, 1.0e-5, false);
    if same_loc {
        face.0.move_location(&location);
    }
    Ok(face)
}


/// `BRepFill::Shell` (`BRepFill.cxx:426-667`).
pub fn shell(wire1: &Wire, wire2: &Wire) -> Result<Shell, String> {
    let mut shell = Shell::new();
    let edges1 = edges_of_wire(wire1);
    let edges2 = edges_of_wire(wire2);
    let closed = wire1.0.is_closed() && wire2.0.is_closed();
    let mut edge4: Option<Edge> = None;
    let mut couture: Option<Edge> = None;
    let n = edges1.len().min(edges2.len());
    for i in 0..n {
        let edge1 = &edges1[i];
        let edge2 = &edges2[i];
        let periodic = edge_is_closed(edge1) && edge_is_closed(edge2);
        let (same_loc, location) = shared_location(edge1, edge2);
        let (c1, v1f, v1l) = prepare_curve(edge1, same_loc)?;
        let (c2, v2f, v2l) = prepare_curve(edge2, same_loc)?;
        let surf = occt_geom::geom_fill::generator_surface(&[c1, c2])?;
        let (u0, u1) = surf.u_range();
        let (v0, v1) = surf.v_range();
        let mut face = TopoBuilder::new().make_face(surf.clone(), &[]);
        GeometryRegistry::global().set_face_tolerance(&face.0, CONFUSION);

        let edge3 = if i == 0 {
            let built = iso_edge(surf.as_ref(), u0, v0, v1, &v1f, &v2f)?;
            if closed {
                couture = Some(built.clone());
            }
            let mut rev = built;
            rev.0.reverse();
            rev
        } else {
            let mut rev = edge4.clone().ok_or_else(|| "BRepFill::Shell".to_string())?;
            rev.0.reverse();
            rev
        };
        let edge4_now = if closed && i + 1 == n {
            couture.clone().ok_or_else(|| "BRepFill::Shell".to_string())?
        } else {
            iso_edge(surf.as_ref(), u1, v0, v1, &v1l, &v2l)?
        };

        let mut wire = Wire::new();
        let builder = TopoBuilder::new();
        builder.add_edge(&mut wire, edge1);
        builder.add_edge(&mut wire, &edge4_now);
        let mut edge2_rev = edge2.clone();
        edge2_rev.0.reverse();
        builder.add_edge(&mut wire, &edge2_rev);
        builder.add_edge(&mut wire, &edge3);
        wire.0.set_closed(true);
        builder.add_wire(&mut face, &wire);
        if same_loc {
            face.0.move_location(&location);
        }
        builder.add_face(&mut shell, &face);

        let reg = GeometryRegistry::global();
        let face_key = GeometryRegistry::shape_key(&face.0);
        bind_profile_pcurve(&reg, edge1, face_key, 0.0, v0, u0, u1);
        bind_profile_pcurve(&reg, edge2, face_key, 0.0, v1, u0, u1);
        if periodic {
            reg.set_edge_pcurves(
                &edge3.0,
                face_key,
                vec![uv_line(u1, 0.0, 0.0, 1.0), uv_line(u0, 0.0, 0.0, 1.0)],
            );
        } else {
            reg.set_edge_pcurve(&edge3.0, face_key, uv_line(u0, 0.0, 0.0, 1.0));
            reg.set_edge_pcurve(&edge4_now.0, face_key, uv_line(u1, 0.0, 0.0, 1.0));
        }
        clear_same_flags(edge1);
        clear_same_flags(edge2);
        clear_same_flags(&edge3);
        clear_same_flags(&edge4_now);
        edge4 = Some(edge4_now);
    }
    shell.0.set_closed(shell_is_closed(&shell.0));
    same_parameter_shape(&shell.0, 1.0e-5, false);
    Ok(shell)
}
