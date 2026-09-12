pub use super::prelude::*;
pub use super::*;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::mesh_box;
    use crate::primitives::{BRepPrimBox, BRepPrimSphere};

    fn temp_subdir(name: &str) -> std::path::PathBuf {
        let mut d = std::env::temp_dir();
        d.push(format!("rwmesh_{}_{}", std::process::id(), name));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    fn triangle_bin() -> Vec<u8> {
        let mut bin = Vec::new();
        for p in [GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.), GpPnt::new(0., 1., 0.)] {
            for c in [p.x() as f32, p.y() as f32, p.z() as f32] {
                bin.extend_from_slice(&c.to_le_bytes());
            }
        }
        for i in [0u32, 1, 2] {
            bin.extend_from_slice(&i.to_le_bytes());
        }
        bin
    }

    /// Single-triangle glTF JSON; `node_extra` is spliced into the node object
    /// (e.g. a `,"matrix":[...]` attribute).
    fn tri_gltf_json(uri: &str, node_extra: &str) -> String {
        format!(
            r#"{{"asset":{{"version":"2.0"}},"scene":0,"scenes":[{{"nodes":[0]}}],"nodes":[{{"mesh":0{node_extra}}}],"meshes":[{{"primitives":[{{"attributes":{{"POSITION":0}},"indices":1}}]}}],"buffers":[{{"uri":"{uri}","byteLength":48}}],"bufferViews":[{{"buffer":0,"byteOffset":0,"byteLength":36,"target":34962}},{{"buffer":0,"byteOffset":36,"byteLength":12,"target":34963}}],"accessors":[{{"bufferView":0,"componentType":5126,"count":3,"type":"VEC3"}},{{"bufferView":1,"componentType":5125,"count":3,"type":"SCALAR"}}]}}"#
        )
    }

    #[test]
    fn read_obj_roundtrip() {
        use occt_core::io::obj::{ObjFace, ObjMesh};
        let mesh = ObjMesh {
            vertices: vec![
                GpPnt::new(0., 0., 0.),
                GpPnt::new(1., 0., 0.),
                GpPnt::new(0., 1., 0.),
                GpPnt::new(1., 1., 0.),
            ],
            texcoords: vec![],
            normals: vec![],
            faces: vec![
                ObjFace { v: vec![0, 1, 2], vt: None, vn: None },
                ObjFace { v: vec![1, 3, 2], vt: None, vn: None },
            ],
        };
        let path = temp_subdir("obj").join("m.obj");
        std::fs::write(&path, occt_core::io::obj::write_obj(&mesh)).unwrap();
        let scene = read_obj_scene(&path.to_string_lossy()).unwrap();
        assert_eq!(mesh_vertex_count(&scene), 4);
        assert_eq!(mesh_triangle_count(&scene), 2);
    }

    #[test]
    fn read_ply_roundtrip() {
        let ply = PlyMesh {
            vertices: vec![
                GpPnt::new(0., 0., 0.),
                GpPnt::new(1., 0., 0.),
                GpPnt::new(1., 1., 0.),
                GpPnt::new(0., 1., 0.),
            ],
            faces: vec![vec![0, 1, 2], vec![0, 2, 3]],
        };
        let path = temp_subdir("ply").join("m.ply");
        std::fs::write(&path, occt_core::io::ply::write_ply(&ply)).unwrap();
        let scene = read_ply_scene(&path.to_string_lossy()).unwrap();
        assert_eq!(mesh_vertex_count(&scene), 4);
        assert_eq!(mesh_triangle_count(&scene), 2);
    }

    #[test]
    fn read_stl_ascii() {
        let text = "solid demo\n\
            facet normal 0 0 1\n  outer loop\n    vertex 0 0 0\n    vertex 1 0 0\n    vertex 0 1 0\n  endloop\nendfacet\n\
            endsolid demo\n";
        let path = temp_subdir("stl_a").join("m.stl");
        std::fs::write(&path, text).unwrap();
        let scene = read_stl_scene(&path.to_string_lossy()).unwrap();
        assert_eq!(mesh_triangle_count(&scene), 1);
        assert_eq!(mesh_vertex_count(&scene), 3);
    }

    #[test]
    fn read_stl_binary() {
        let stl = StlMesh {
            triangles: vec![
                [GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.), GpPnt::new(0., 1., 0.)],
                [GpPnt::new(1., 0., 0.), GpPnt::new(1., 1., 0.), GpPnt::new(0., 1., 0.)],
            ],
            normals: vec![],
        };
        let path = temp_subdir("stl_b").join("m.stl");
        std::fs::write(&path, occt_core::io::stl::write_binary_stl(&stl)).unwrap();
        let scene = read_stl_scene(&path.to_string_lossy()).unwrap();
        assert_eq!(mesh_triangle_count(&scene), 2);
        // Shared vertices deduplicated.
        assert_eq!(mesh_vertex_count(&scene), 4);
    }

    #[test]
    fn read_gltf_embedded_base64() {
        let bin = triangle_bin();
        let uri = format!("data:application/octet-stream;base64,{}", crate::gltf::base64_encode(&bin));
        let path = temp_subdir("gltf_b64").join("m.gltf");
        std::fs::write(&path, tri_gltf_json(&uri, "")).unwrap();
        let scene = read_gltf_scene(&path.to_string_lossy()).unwrap();
        assert_eq!(mesh_vertex_count(&scene), 3);
        assert_eq!(mesh_triangle_count(&scene), 1);
    }

    #[test]
    fn read_gltf_external_bin() {
        let dir = temp_subdir("gltf_ext");
        std::fs::write(dir.join("mesh.bin"), triangle_bin()).unwrap();
        let path = dir.join("m.gltf");
        std::fs::write(&path, tri_gltf_json("mesh.bin", "")).unwrap();
        let scene = read_gltf_scene(&path.to_string_lossy()).unwrap();
        assert_eq!(mesh_vertex_count(&scene), 3);
        assert_eq!(mesh_triangle_count(&scene), 1);
    }

    #[test]
    fn scene_to_compound_box() {
        let shape = BRepPrimBox::make_box(2.0, 2.0, 2.0).solid.0;
        let scene = mesh_from_shape(&shape, 0.5);
        let comp = scene_to_compound(&scene, 1e-7).expect("compound");
        assert!(crate::shape_mesh::count_faces(&comp) > 0, "no faces in compound");
        assert!(mesh_vertex_count(&scene) > 0);

        // Closedness: a clean watertight box mesh converts to a solid.
        let m = mesh_box((GpPnt::zero(), GpPnt::new(2.0, 2.0, 2.0)));
        let scene2 = MeshScene {
            name: "box".into(),
            nodes: vec![MeshNode {
                name: "box".into(),
                vertices: m.vertices,
                triangles: m.triangles.iter().map(|t| (t.n0, t.n1, t.n2)).collect(),
                transform: None,
                material: None,
                uv: None,
                texture_path: None,
            }],
            materials: Vec::new(),
        };
        let s = scene_to_shape(&scene2, 0, 1e-7).expect("box shape");
        assert!(s.is_solid(), "watertight box should produce a solid");
    }

    #[test]
    fn scene_bounds_ok() {
        let scene = MeshScene {
            name: "s".into(),
            nodes: vec![MeshNode {
                name: "n".into(),
                vertices: vec![
                    GpPnt::new(1., 2., 3.),
                    GpPnt::new(-1., -2., -3.),
                    GpPnt::new(4., 0., 5.),
                ],
                triangles: vec![(0, 1, 2)],
                transform: None,
                material: None,
                uv: None,
                texture_path: None,
            }],
            materials: Vec::new(),
        };
        let (lo, hi) = scene_bounds(&scene).expect("bounds");
        for v in &scene.nodes[0].vertices {
            assert!(lo.x() <= v.x() && v.x() <= hi.x(), "x in bounds");
            assert!(lo.y() <= v.y() && v.y() <= hi.y(), "y in bounds");
            assert!(lo.z() <= v.z() && v.z() <= hi.z(), "z in bounds");
        }
        assert!(lo.is_equal(&GpPnt::new(-1., -2., -3.)));
        assert!(hi.is_equal(&GpPnt::new(4., 2., 5.)));
        assert!(scene_bounds(&MeshScene { name: "e".into(), nodes: vec![], materials: Vec::new() }).is_none());
    }

    #[test]
    fn convert_obj_to_ply() {
        let src = temp_subdir("conv").join("m.obj");
        let dst = temp_subdir("conv").join("m.ply");
        std::fs::write(
            &src,
            "v 0 0 0\nv 1 0 0\nv 0 1 0\nv 1 1 0\nf 1 2 3\nf 2 4 3\n",
        )
        .unwrap();
        convert_mesh_format(&src.to_string_lossy(), &dst.to_string_lossy()).unwrap();
        let scene = read_ply_scene(&dst.to_string_lossy()).unwrap();
        assert_eq!(mesh_triangle_count(&scene), 2);
        assert_eq!(mesh_vertex_count(&scene), 4);
    }

    #[test]
    fn roundtrip_shape_mesh() {
        let shape = BRepPrimSphere::make_sphere(1.0).solid.0;
        let scene = mesh_from_shape(&shape, 0.2);
        let n_before = mesh_vertex_count(&scene);
        assert!(n_before > 100, "sphere mesh too coarse: {n_before}");

        let s = scene_to_shape(&scene, 0, 1e-7).expect("sphere shape");
        // Each input triangle becomes one face.
        assert_eq!(crate::shape_mesh::count_faces(&s), mesh_triangle_count(&scene));

        // Vertex count preserved by the conversion (mesh is a single closed grid).
        let tri = Triangulation::new(
            scene.nodes[0].vertices.clone(),
            scene.nodes[0].triangles.iter().map(|&(a, b, c)| Triangle::new(a, b, c)).collect(),
        );
        let b = triangulation_to_brep(&tri);
        // The sphere grid collapses its pole rows (nu points → 1 each), so the
        // deduped conversion keeps slightly fewer vertices — but it must stay
        // close to the source count.
        let loss = n_before.saturating_sub(b.vertices.len());
        assert!(
            loss <= n_before / 4,
            "vertex count drifted: before {n_before} after {}",
            b.vertices.len()
        );
    }

    #[test]
    fn node_transform_applied_gltf() {
        let bin = triangle_bin();
        let uri = format!("data:application/octet-stream;base64,{}", crate::gltf::base64_encode(&bin));
        let dir = temp_subdir("gltf_trs");
        let path = dir.join("m.gltf");
        let json = tri_gltf_json(&uri, r#","matrix":[1,0,0,0,0,1,0,0,0,0,1,0,5,0,0,1]"#);
        std::fs::write(&path, json).unwrap();
        let scene = read_gltf_scene(&path.to_string_lossy()).unwrap();
        let node = &scene.nodes[0];
        assert!((node.vertices[0].x() - 5.0).abs() < 1e-6, "v0.x = {}", node.vertices[0].x());
        assert!((node.vertices[1].x() - 6.0).abs() < 1e-6, "v1.x = {}", node.vertices[1].x());
        assert!((node.vertices[2].y() - 1.0).abs() < 1e-6, "v2.y = {}", node.vertices[2].y());
        assert!(node.transform.is_some(), "node transform recorded");
    }

    #[test]
    fn obj_groups_become_nodes() {
        let obj = "o partA\nv 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\no partB\nv 5 5 5\nv 6 5 5\nv 5 6 5\nf 4 5 6\n";
        let scene = parse_obj_scene(obj, "groups").unwrap();
        assert_eq!(scene.nodes.len(), 2, "two object groups -> two nodes");
        assert_eq!(mesh_triangle_count(&scene), 2);
        assert_eq!(scene.nodes[0].name, "partA");
        assert_eq!(scene.nodes[1].name, "partB");
        assert_eq!(mesh_vertex_count(&scene), 6);
    }

    #[test]
    fn obj_mtl_material_parsed() {
        let dir = temp_subdir("obj_mtl");
        std::fs::write(
            dir.join("colors.mtl"),
            "newmtl red\nKd 1 0 0\nKs 0.2 0.2 0.2\nd 0.5\n",
        )
        .unwrap();
        let obj = "mtllib colors.mtl\nusemtl red\no part\nv 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n";
        let path = dir.join("m.obj");
        std::fs::write(&path, obj).unwrap();
        let scene = read_obj_scene(&path.to_string_lossy()).unwrap();
        assert_eq!(scene.materials.len(), 1, "mtl loaded into library");
        assert_eq!(scene.materials[0].name, "red");
        let (r, g, b) = scene.materials[0].diffuse;
        assert!((r - 1.0).abs() < 1e-6 && g.abs() < 1e-6 && b.abs() < 1e-6, "Kd 1 0 0");
        assert!((scene.materials[0].opacity - 0.5).abs() < 1e-6, "d 0.5");
        let mat = scene.nodes[0].material.as_ref().expect("usemtl applied to node");
        let (r, g, b) = mat.diffuse;
        assert!((r - 1.0).abs() < 1e-6 && g.abs() < 1e-6 && b.abs() < 1e-6, "node diffuse red");
    }

    #[test]
    fn obj_uv_coordinates() {
        let obj = "v 0 0 0\nv 1 0 0\nv 0 1 0\nvt 0 0\nvt 1 0\nvt 0 1\nf 1/1 2/2 3/3\n";
        let scene = parse_obj_scene(obj, "uv").unwrap();
        let node = &scene.nodes[0];
        let uv = node.uv.as_ref().expect("uv present");
        assert_eq!(uv.len(), 3, "one uv per vertex");
        assert!((uv[0].x() - 0.0).abs() < 1e-9 && (uv[0].y() - 0.0).abs() < 1e-9);
        assert!((uv[1].x() - 1.0).abs() < 1e-9 && (uv[1].y() - 0.0).abs() < 1e-9);
        assert!((uv[2].x() - 0.0).abs() < 1e-9 && (uv[2].y() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn gltf_material_color() {
        let bin = triangle_bin();
        let uri = format!("data:application/octet-stream;base64,{}", crate::gltf::base64_encode(&bin));
        let json = format!(
            r#"{{"asset":{{"version":"2.0"}},"scene":0,"scenes":[{{"nodes":[0]}}],"nodes":[{{"mesh":0}}],"meshes":[{{"primitives":[{{"attributes":{{"POSITION":0}},"indices":1,"material":0}}]}}],"materials":[{{"name":"redmat","pbrMetallicRoughness":{{"baseColorFactor":[1,0,0,1]}}}}],"buffers":[{{"uri":"{uri}","byteLength":48}}],"bufferViews":[{{"buffer":0,"byteOffset":0,"byteLength":36,"target":34962}},{{"buffer":0,"byteOffset":36,"byteLength":12,"target":34963}}],"accessors":[{{"bufferView":0,"componentType":5126,"count":3,"type":"VEC3"}},{{"bufferView":1,"componentType":5125,"count":3,"type":"SCALAR"}}]}}"#
        );
        let path = temp_subdir("gltf_mat").join("m.gltf");
        std::fs::write(&path, json).unwrap();
        let scene = read_gltf_scene(&path.to_string_lossy()).unwrap();
        let mat = scene.nodes[0].material.as_ref().expect("gltf material attached");
        let (r, g, b) = mat.diffuse;
        assert!((r - 1.0).abs() < 1e-6 && g.abs() < 1e-6 && b.abs() < 1e-6, "baseColorFactor [1,0,0,1]");
        assert_eq!(scene.materials.len(), 1, "gltf materials copied to library");
        assert_eq!(scene.materials[0].name, "redmat");
    }

    fn box_vrml() -> &'static str {
        "#VRML V2.0 utf8\n\
         Shape {\n\
           geometry IndexedFaceSet {\n\
             coord Coordinate {\n\
               point [\n\
                 0 0 0, 1 0 0, 1 1 0, 0 1 0,\n\
                 0 0 1, 1 0 1, 1 1 1, 0 1 1\n\
               ]\n\
             }\n\
             coordIndex [\n\
               0,1,2,3,-1,  4,7,6,5,-1,  0,4,5,1,-1,\n\
               1,5,6,2,-1,  2,6,7,3,-1,  3,7,4,0,-1\n\
             ]\n\
           }\n\
         }\n"
    }

    #[test]
    fn vrml_simple_box() {
        let scene = parse_vrml_scene(box_vrml(), "cube").unwrap();
        assert_eq!(scene.nodes.len(), 1, "one shape -> one node");
        assert_eq!(mesh_triangle_count(&scene), 12, "six quads fan to twelve tris");
        assert_eq!(mesh_vertex_count(&scene), 8, "cube has eight points");
    }

    #[test]
    fn vrml_transform_nesting() {
        let vrml = "#VRML V2.0 utf8\n\
            Transform {\n\
              translation 2 0 0\n\
              children [\n\
                Shape {\n\
                  geometry IndexedFaceSet {\n\
                    coord Coordinate { point [ 0 0 0, 1 0 0, 0 1 0 ] }\n\
                    coordIndex [ 0,1,2,-1 ]\n\
                  }\n\
                }\n\
              ]\n\
            }\n";
        let scene = parse_vrml_scene(vrml, "moved").unwrap();
        let node = &scene.nodes[0];
        assert_eq!(node.vertices.len(), 3);
        // Vertices (0,0,0), (1,0,0), (0,1,0) shifted by translation [2,0,0].
        assert!((node.vertices[0].x() - 2.0).abs() < 1e-9, "v0.x = {}", node.vertices[0].x());
        assert!((node.vertices[1].x() - 3.0).abs() < 1e-9, "v1.x = {}", node.vertices[1].x());
        assert!((node.vertices[2].x() - 2.0).abs() < 1e-9, "v2.x = {}", node.vertices[2].x());
        assert!((node.vertices[2].y() - 1.0).abs() < 1e-9, "v2.y = {}", node.vertices[2].y());
        assert_eq!(mesh_triangle_count(&scene), 1);
    }

    #[test]
    fn vrml_material_diffuse() {
        let vrml = "#VRML V2.0 utf8\n\
            Shape {\n\
              appearance Appearance {\n\
                material Material { diffuseColor 0 0 1 }\n\
              }\n\
              geometry IndexedFaceSet {\n\
                coord Coordinate { point [ 0 0 0, 1 0 0, 0 1 0 ] }\n\
                coordIndex [ 0,1,2,-1 ]\n\
              }\n\
            }\n";
        let scene = parse_vrml_scene(vrml, "colored").unwrap();
        let mat = scene.nodes[0].material.as_ref().expect("diffuse material attached");
        let (r, g, b) = mat.diffuse;
        assert!(r.abs() < 1e-6 && g.abs() < 1e-6 && (b - 1.0).abs() < 1e-6, "diffuseColor 0 0 1");
    }

    #[test]
    fn convert_wrl_to_obj() {
        let mesh = crate::mesh::mesh_box((GpPnt::zero(), GpPnt::new(1.0, 1.0, 1.0)));
        let wrl_text = crate::vrml::write_vrml(&mesh, "box");
        let dir = temp_subdir("wrl2obj");
        let src = dir.join("m.wrl");
        let dst = dir.join("m.obj");
        std::fs::write(&src, wrl_text).unwrap();
        convert_mesh_format(&src.to_string_lossy(), &dst.to_string_lossy()).unwrap();
        let scene = read_obj_scene(&dst.to_string_lossy()).unwrap();
        assert_eq!(mesh_triangle_count(&scene), 12, "box has 12 tris after obj roundtrip");
        assert_eq!(mesh_vertex_count(&scene), 8);
    }

    #[test]
    fn material_lookup() {
        let scene = MeshScene {
            name: "s".into(),
            nodes: Vec::new(),
            materials: vec![Material {
                name: "brass".into(),
                diffuse: (0.8, 0.6, 0.2),
                ..Material::default()
            }],
        };
        assert!(material_from_name(&scene, "brass").is_some(), "known material found");
        let m = material_from_name(&scene, "brass").unwrap();
        assert!((m.diffuse.0 - 0.8).abs() < 1e-9);
        assert!(material_from_name(&scene, "nope").is_none(), "unknown material -> None");
    }

    #[test]
    fn mtl_file_roundtrip_uv() {
        let dir = temp_subdir("mtl_uv");
        std::fs::write(dir.join("m.mtl"), "newmtl blue\nKd 0 0 1\n").unwrap();
        let obj = "mtllib m.mtl\nusemtl blue\no part\n\
                   v 0 0 0\nv 1 0 0\nv 0 1 0\n\
                   vt 0.25 0.25\nvt 0.75 0.25\nvt 0.25 0.75\n\
                   f 1/1 2/2 3/3\n";
        let path = dir.join("m.obj");
        std::fs::write(&path, obj).unwrap();
        let scene = read_obj_scene(&path.to_string_lossy()).unwrap();
        let node = &scene.nodes[0];
        let mat = node.material.as_ref().expect("material preserved");
        assert!((mat.diffuse.2 - 1.0).abs() < 1e-6, "blue diffuse");
        let uv = node.uv.as_ref().expect("uv preserved alongside material");
        assert_eq!(uv.len(), 3);
        assert!((uv[1].x() - 0.75).abs() < 1e-9);
        assert_eq!(scene.materials.len(), 1);
    }

    #[test]
    fn no_material_is_none() {
        let obj = "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n";
        let scene = parse_obj_scene(obj, "plain").unwrap();
        assert!(scene.nodes[0].material.is_none(), "no usemtl -> no material");
        assert!(scene.nodes[0].uv.is_none(), "no vt -> no uv");
    }

    #[test]
    fn ppm_checkerboard_uv() {
        let scene = MeshScene {
            name: "s".into(),
            nodes: vec![MeshNode {
                name: "n".into(),
                vertices: vec![GpPnt::new(0., 0., 0.)],
                triangles: Vec::new(),
                transform: None,
                material: None,
                uv: Some(vec![GpPnt2d::new(0.1, 0.1)]),
                texture_path: None,
            }],
            materials: Vec::new(),
        };
        let ppm = scene_with_texture_to_ppm(&scene, 0, 2, 2).unwrap();
        assert!(ppm.starts_with("P3\n2 2\n255\n"), "P3 header with size");
        // Top-left quadrant (u<0.5, v<0.5) is red.
        assert!(ppm.contains("255 0 0\n"), "uv quadrant red present");
        // Without UVs the image falls back to material diffuse.
        let plain = MeshScene {
            name: "p".into(),
            nodes: vec![MeshNode {
                name: "n".into(),
                vertices: vec![GpPnt::new(0., 0., 0.)],
                triangles: Vec::new(),
                transform: None,
                material: None,
                uv: None,
                texture_path: None,
            }],
            materials: Vec::new(),
        };
        let ppm2 = scene_with_texture_to_ppm(&plain, 0, 1, 1).unwrap();
        assert!(ppm2.contains("179 179 179\n"), "default gray 0.7 -> 179");
        assert!(scene_with_texture_to_ppm(&plain, 5, 1, 1).is_err(), "bad index errors");
    }
}
