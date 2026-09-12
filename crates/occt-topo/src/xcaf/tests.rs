use super::prelude::*;
use super::*;
    use crate::primitives::{BRepPrimBox, BRepPrimSphere};
    use occt_core::quantity::Color;

    fn two_shape_model() -> BRepModel {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let s = BRepPrimSphere::make_sphere(0.5);
        let mut model = BRepModel::new();
        model.add_with_color("RedBox", b.solid.0.clone(), Color::RED);
        model.add("GreyBall", s.solid.0.clone());
        model.find_mut("GreyBall").unwrap().layer = "L1".into();
        model
    }

    #[test]
    fn model_to_assembly_has_both_products() {
        let model = two_shape_model();
        let asm = model_to_step_assembly(&model);
        assert!(asm.find("RedBox").is_some());
        assert!(asm.find("GreyBall").is_some());
        assert_eq!(asm.colors.get("RedBox"), Some(&(1.0, 0.0, 0.0)));
        assert_eq!(asm.layers.get("GreyBall").map(|s| s.as_str()), Some("L1"));
        assert_eq!(asm.products.len(), 2);
    }

    #[test]
    fn write_contains_names() {
        let model = two_shape_model();
        let step = write_step_with_metadata(&model);
        assert!(step.contains("RedBox"), "step:\n{step}");
        assert!(step.contains("GreyBall"));
        assert!(step.contains(BLOCK_BEGIN));
    }

    #[test]
    fn extract_roundtrips_names_and_colors() {
        let model = two_shape_model();
        let step = write_step_with_metadata(&model);
        let asm = extract_metadata_from_step(&step);
        assert!(asm.find("RedBox").is_some());
        assert!(asm.find("GreyBall").is_some());
        assert_eq!(asm.colors.get("RedBox"), Some(&(1.0, 0.0, 0.0)));
        assert_eq!(asm.layers.get("GreyBall").map(|s| s.as_str()), Some("L1"));
        assert_eq!(asm.products.len(), 2);
    }

    #[test]
    fn full_roundtrip_succeeds() {
        let model = two_shape_model();
        model_with_metadata_step_roundtrip(&model).expect("roundtrip ok");
    }
