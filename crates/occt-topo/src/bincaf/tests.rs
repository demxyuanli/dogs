use super::prelude::*;
use super::*;
    use crate::builder::TopoBuilder;
    use crate::primitives::{BRepPrimBox, BRepPrimSphere};
    use crate::topo_tools_full::faces_of;

    fn attr(kind: &str, value: &str) -> XcafAttribute {
        XcafAttribute {
            kind: kind.into(),
            value: value.into(),
        }
    }

    #[test]
    fn roundtrip_empty() {
        let doc = BinXcaf {
            root: BinXcafEntry::default(),
        };
        let bytes = serialize_bincaf(&doc).expect("serialize");
        let got = deserialize_bincaf(&bytes).expect("deserialize");
        assert!(got.root.shape.is_none());
        assert!(got.root.attributes.is_empty());
        assert!(got.root.children.is_empty());
    }

    #[test]
    fn roundtrip_box_solid() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let doc = BinXcaf {
            root: BinXcafEntry {
                shape: Some(b.solid.0),
                ..Default::default()
            },
        };
        let bytes = serialize_bincaf(&doc).expect("serialize");
        let got = deserialize_bincaf(&bytes).expect("deserialize");
        let shape = got.root.shape.expect("shape");
        assert_eq!(shape.shape_type(), ShapeType::Solid);
        assert_eq!(faces_of(&shape).len(), 6, "face count");
    }

    #[test]
    fn roundtrip_attributes() {
        let doc = BinXcaf {
            root: BinXcafEntry {
                shape: None,
                attributes: vec![
                    attr("name", "Root"),
                    attr("color", "1.0,0.0,0.0"),
                    attr("layer", "L1"),
                ],
                children: vec![BinXcafEntry {
                    shape: None,
                    attributes: vec![attr("name", "Child")],
                    children: vec![],
                }],
            },
        };
        let bytes = serialize_bincaf(&doc).expect("serialize");
        let got = deserialize_bincaf(&bytes).expect("deserialize");
        assert_eq!(got.root.attributes.len(), 3);
        assert_eq!(got.root.attributes[0], attr("name", "Root"));
        assert_eq!(got.root.attributes[1], attr("color", "1.0,0.0,0.0"));
        assert_eq!(got.root.attributes[2], attr("layer", "L1"));
        assert_eq!(got.root.children.len(), 1);
        assert_eq!(got.root.children[0].attributes, vec![attr("name", "Child")]);
    }

    #[test]
    fn roundtrip_nested_tree() {
        let grandchild = BinXcafEntry {
            shape: None,
            attributes: vec![attr("name", "GC")],
            children: vec![],
        };
        let child = BinXcafEntry {
            shape: None,
            attributes: vec![],
            children: vec![grandchild],
        };
        let doc = BinXcaf {
            root: BinXcafEntry {
                shape: None,
                attributes: vec![],
                children: vec![child],
            },
        };
        let bytes = serialize_bincaf(&doc).expect("serialize");
        let got = deserialize_bincaf(&bytes).expect("deserialize");
        assert_eq!(got.root.children.len(), 1);
        assert_eq!(got.root.children[0].children.len(), 1);
        assert_eq!(got.root.children[0].children[0].attributes[0].value, "GC");
    }

    #[test]
    fn bad_magic_rejected() {
        assert!(deserialize_bincaf(b"NOTBINXCAF!!").is_err());
        assert!(deserialize_bincaf(&[0u8; 16]).is_err());
    }

    #[test]
    fn truncated_rejected() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let doc = BinXcaf {
            root: BinXcafEntry {
                shape: Some(b.solid.0),
                attributes: vec![attr("name", "Box")],
                children: vec![],
            },
        };
        let bytes = serialize_bincaf(&doc).expect("serialize");
        let cut = bytes.len() * 3 / 5;
        // Must not panic: either Err on a truncated read or a smaller valid doc.
        let _ = deserialize_bincaf(&bytes[..cut]);
    }

    #[test]
    fn roundtrip_sphere_solid() {
        let s = BRepPrimSphere::make_sphere(2.5);
        let doc = BinXcaf {
            root: BinXcafEntry {
                shape: Some(s.solid.0),
                ..Default::default()
            },
        };
        let bytes = serialize_bincaf(&doc).expect("serialize");
        let got = deserialize_bincaf(&bytes).expect("deserialize");
        let shape = got.root.shape.expect("shape");
        assert_eq!(shape.shape_type(), ShapeType::Solid);
        assert_eq!(faces_of(&shape).len(), 1, "sphere face count");
    }

    #[test]
    fn roundtrip_compound_solids() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let s = BRepPrimSphere::make_sphere(1.0);
        let bld = TopoBuilder::new();
        let comp = bld.make_compound_of(&[b.solid.0, s.solid.0]);
        let doc = BinXcaf {
            root: BinXcafEntry {
                shape: Some(comp.0),
                ..Default::default()
            },
        };
        let bytes = serialize_bincaf(&doc).expect("serialize");
        let got = deserialize_bincaf(&bytes).expect("deserialize");
        let shape = got.root.shape.expect("shape");
        assert_eq!(shape.shape_type(), ShapeType::Compound);
        let kids: Vec<TopoShape> = shape
            .tshape
            .read()
            .unwrap()
            .children
            .clone();
        assert_eq!(kids.len(), 2);
        assert_eq!(kids[0].shape_type(), ShapeType::Solid);
        assert_eq!(kids[1].shape_type(), ShapeType::Solid);
        assert_eq!(faces_of(&kids[0]).len(), 6);
        assert_eq!(faces_of(&kids[1]).len(), 1);
    }

    #[test]
    fn file_io_roundtrip() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let doc = BinXcaf {
            root: BinXcafEntry {
                shape: Some(b.solid.0),
                attributes: vec![attr("name", "Box")],
                children: vec![],
            },
        };
        let path = std::env::temp_dir().join("occt_bincaf_test.xbf");
        let p = path.to_str().unwrap();
        write_bincaf_file(&doc, p).expect("write file");
        let got = read_bincaf_file(p).expect("read file");
        assert_eq!(got.root.attributes[0].value, "Box");
        assert_eq!(got.root.shape.expect("shape").shape_type(), ShapeType::Solid);
        std::fs::remove_file(p).ok();
    }
