use super::prelude::*;
use super::*;
    use crate::primitives::{BRepPrimBox, BRepPrimSphere};
    use crate::topo_tools_full::faces_of;

    fn attr(kind: &str, value: &str) -> XmlAttribute {
        XmlAttribute {
            kind: kind.into(),
            value: value.into(),
        }
    }

    fn doc_with_attrs() -> XmlXcafDoc {
        XmlXcafDoc {
            version: "1.0".into(),
            root: XmlEntry {
                name: "Root".into(),
                shape: None,
                attributes: vec![
                    attr("name", "Root"),
                    attr("color", "1.0,0.0,0.0"),
                    attr("layer", "L1"),
                ],
                children: vec![XmlEntry {
                    name: "Child".into(),
                    shape: None,
                    attributes: vec![attr("name", "Child")],
                    children: vec![],
                }],
            },
        }
    }

    #[test]
    fn xml_header_and_escape() {
        let s = "a&b<c>d\"e'f";
        let esc = escape_xml(s);
        assert_eq!(esc, "a&amp;b&lt;c&gt;d&quot;e&apos;f");
        assert_eq!(unescape_xml(&esc).unwrap(), s);
        assert_eq!(escape_xml("1 < 2 && 3 > 0"), "1 &lt; 2 &amp;&amp; 3 &gt; 0");
    }

    #[test]
    fn to_xml_wellformed() {
        let doc = doc_with_attrs();
        let xml = to_xml(&doc).unwrap();
        assert!(xml.starts_with("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"), "header:\n{xml}");
        assert!(xml.contains("<xcaf version=\"1.0\">"));
        // <entry> and <shape> are written with matching open/close tags.
        for tag in ["entry", "shape"] {
            let open = xml.matches(&format!("<{tag}")).count();
            let close = xml.matches(&format!("</{tag}>")).count();
            assert_eq!(open, close, "unbalanced <{tag}>");
        }
        // <attribute> is self-closing: it must never appear with a close tag.
        assert!(xml.contains("<attribute"));
        assert_eq!(xml.matches("</attribute>").count(), 0);
    }

    #[test]
    fn roundtrip_attributes() {
        let doc = doc_with_attrs();
        let got = from_xml(&to_xml(&doc).unwrap()).unwrap();
        assert_eq!(got.version, "1.0");
        assert_eq!(got.root.name, "Root");
        assert_eq!(got.root.attributes, doc.root.attributes);
        assert_eq!(got.root.children.len(), 1);
        assert_eq!(got.root.children[0].attributes, vec![attr("name", "Child")]);
        assert!(got.root.shape.is_none());
    }

    #[test]
    fn roundtrip_nested_tree() {
        let doc = XmlXcafDoc {
            version: "1.0".into(),
            root: XmlEntry {
                name: "R".into(),
                shape: None,
                attributes: vec![],
                children: vec![
                    XmlEntry {
                        name: "C1".into(),
                        shape: None,
                        attributes: vec![attr("name", "c1")],
                        children: vec![XmlEntry {
                            name: "GC".into(),
                            shape: None,
                            attributes: vec![attr("name", "gc")],
                            children: vec![],
                        }],
                    },
                    XmlEntry {
                        name: "C2".into(),
                        shape: None,
                        attributes: vec![],
                        children: vec![],
                    },
                ],
            },
        };
        let got = from_xml(&to_xml(&doc).unwrap()).unwrap();
        assert_eq!(got.root.children.len(), 2);
        assert_eq!(got.root.children[0].name, "C1");
        assert_eq!(got.root.children[0].children.len(), 1);
        assert_eq!(got.root.children[0].children[0].name, "GC");
        assert_eq!(got.root.children[0].children[0].attributes[0].value, "gc");
        assert_eq!(got.root.children[1].name, "C2");
    }

    #[test]
    fn roundtrip_box_shape() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let doc = XmlXcafDoc {
            version: "1.0".into(),
            root: XmlEntry {
                name: "Box".into(),
                shape: Some(b.solid.0),
                attributes: vec![],
                children: vec![],
            },
        };
        let xml = to_xml(&doc).unwrap();
        assert!(xml.contains("<shape type=\"Solid\">"), "solid in:\n{xml}");
        let got = from_xml(&xml).unwrap();
        let shape = got.root.shape.expect("shape");
        assert_eq!(shape.shape_type(), ShapeType::Solid);
        assert_eq!(faces_of(&shape).len(), 6, "box face count");
    }

    #[test]
    fn roundtrip_sphere_shape() {
        let s = BRepPrimSphere::make_sphere(2.5);
        let doc = XmlXcafDoc {
            version: "1.0".into(),
            root: XmlEntry {
                name: "Sphere".into(),
                shape: Some(s.solid.0),
                attributes: vec![],
                children: vec![],
            },
        };
        let xml = to_xml(&doc).unwrap();
        assert!(xml.contains("kind=\"sphere\""), "sphere surface in:\n{xml}");
        let got = from_xml(&xml).unwrap();
        let shape = got.root.shape.expect("shape");
        assert_eq!(shape.shape_type(), ShapeType::Solid);
        assert_eq!(faces_of(&shape).len(), 1, "sphere face count");
    }

    #[test]
    fn malformed_xml_rejected() {
        assert!(from_xml("not xml at all").is_err());
        assert!(from_xml("<?xml version=\"1.0\"?><xcaf>").is_err());
        assert!(from_xml("<xcaf><entry name=\"a\"></entry>").is_err());
        assert!(from_xml("<xcaf><entry name=\"a\"><foo/></entry></xcaf>").is_err());
        assert!(from_xml("<xcaf><entry name=\"a\"/><entry name=\"b\"/></xcaf>").is_err());
        assert!(from_xml("<xcaf><entry name=\"a\"></entry></xcaf> junk").is_err());
        assert!(from_xml("<xcaf><entry name=\"a\"><entry></entry></xcf></xcaf>").is_err());
        assert!(unescape_xml("&bogus;").is_err());
        assert!(unescape_xml("a & b").is_err());
    }

    #[test]
    fn xml_escape_roundtrip() {
        let tricky = "Widget <A> & \"B\" 'C' > 5";
        let doc = XmlXcafDoc {
            version: "1.0".into(),
            root: XmlEntry {
                name: "Root".into(),
                shape: None,
                attributes: vec![attr("name", tricky)],
                children: vec![],
            },
        };
        let xml = to_xml(&doc).unwrap();
        assert!(xml.contains("&lt;A&gt; &amp; &quot;B&quot; &apos;C&apos; &gt; 5"), "escaped:\n{xml}");
        let got = from_xml(&xml).unwrap();
        assert_eq!(got.root.attributes[0].value, tricky);
    }

    #[test]
    fn file_roundtrip() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let doc = XmlXcafDoc {
            version: "1.0".into(),
            root: XmlEntry {
                name: "Box".into(),
                shape: Some(b.solid.0),
                attributes: vec![attr("name", "Box")],
                children: vec![],
            },
        };
        let path = std::env::temp_dir().join("occt_xmlcaf_test.xml");
        let p = path.to_str().unwrap();
        write_xml_file(&doc, p).unwrap();
        let got = read_xml_file(p).unwrap();
        assert_eq!(got.root.attributes[0].value, "Box");
        let shape = got.root.shape.expect("shape");
        assert_eq!(shape.shape_type(), ShapeType::Solid);
        assert_eq!(faces_of(&shape).len(), 6);
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn empty_doc_ok() {
        let doc = XmlXcafDoc::default();
        let got = from_xml(&to_xml(&doc).unwrap()).unwrap();
        assert_eq!(got.version, "1.0");
        assert_eq!(got.root.name, "");
        assert!(got.root.shape.is_none());
        assert!(got.root.attributes.is_empty());
        assert!(got.root.children.is_empty());
    }

    #[test]
    fn shape_element_helpers_roundtrip() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let shape = b.solid.0.clone();
        let xml = shape_to_xml_element(&shape);
        // `tag` is the opening tag; `inner` is everything up to the final
        // closing `</shape>` (the root shape's own close tag).
        let tag = xml.lines().next().unwrap_or("").to_string();
        let open_end = xml.find('>').unwrap();
        let close_start = xml.rfind("</shape>").unwrap();
        let inner = &xml[open_end + 1..close_start];
        let got = shape_from_xml_element(&tag, inner).unwrap();
        assert_eq!(got.shape_type(), ShapeType::Solid);
        assert_eq!(faces_of(&got).len(), 6);
    }
