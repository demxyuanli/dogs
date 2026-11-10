use super::super::delaun_types::{DelaunTriangle, DelaunVertex, VertexState};
use super::{DelaunDataStructure, DelaunSelector};

fn v(u: f64, w: f64) -> DelaunVertex {
    DelaunVertex::new_parametric(u, w, VertexState::Free)
}

#[test]
fn add_node_dedupes_coincident() {
    let mut ds = DelaunDataStructure::new(10);
    let a = ds.add_node(v(1.0, 2.0));
    let b = ds.add_node(v(1.0 + 1e-8, 2.0)); // within Confusion (1e-7)
    let c = ds.add_node(v(5.0, 5.0));
    assert_eq!(a, b);
    assert_ne!(a, c);
    assert_eq!(ds.nb_nodes(), 2);
}

#[test]
fn add_link_detects_reverse_orientation() {
    let mut ds = DelaunDataStructure::new(10);
    let a = ds.add_node(v(0.0, 0.0));
    let b = ds.add_node(v(1.0, 0.0));
    let l1 = ds.add_link(a, b, VertexState::Free);
    assert!(l1 > 0);
    let l2 = ds.add_link(b, a, VertexState::Free); // same link, reversed
    assert_eq!(l1, -l2);
    assert_eq!(ds.links_connected_to(a), &[l1]);
    assert_eq!(ds.links_connected_to(b), &[l1]);
}

#[test]
fn element_link_registration_and_removal() {
    let mut ds = DelaunDataStructure::new(10);
    let a = ds.add_node(v(0.0, 0.0));
    let b = ds.add_node(v(1.0, 0.0));
    let c = ds.add_node(v(0.0, 1.0));
    let ab = ds.add_link(a, b, VertexState::Free);
    let bc = ds.add_link(b, c, VertexState::Free);
    let ca = ds.add_link(c, a, VertexState::Free);
    let tri = DelaunTriangle::new([ab, bc, ca], [a, b, c]);
    let id = ds.add_element(tri);
    assert_eq!(ds.elements_of_domain().len(), 1);
    let pair = ds.elements_connected_to(ab);
    assert_eq!(pair.first_index(), id);
    ds.remove_element(id);
    assert!(ds.elements_of_domain().is_empty());
    assert!(ds.elements_connected_to(ab).is_empty());
    assert_eq!(ds.element_movability(id), VertexState::Deleted);
}

#[test]
fn remove_link_requires_no_elements() {
    let mut ds = DelaunDataStructure::new(10);
    let a = ds.add_node(v(0.0, 0.0));
    let b = ds.add_node(v(1.0, 0.0));
    let ab = ds.add_link(a, b, VertexState::Free);
    assert_eq!(ds.links_of_domain().len(), 1);
    ds.remove_link(ab, false);
    assert!(ds.links_of_domain().is_empty());
    assert_eq!(ds.link_movability(ab), VertexState::Deleted);
}

#[test]
fn clear_domain_removes_elements_and_free_links() {
    let mut ds = DelaunDataStructure::new(10);
    let a = ds.add_node(v(0.0, 0.0));
    let b = ds.add_node(v(1.0, 0.0));
    let c = ds.add_node(v(0.0, 1.0));
    let ab = ds.add_link(a, b, VertexState::Free);
    let bc = ds.add_link(b, c, VertexState::Free);
    let ca = ds.add_link(c, a, VertexState::Free);
    ds.add_element(DelaunTriangle::new([ab, bc, ca], [a, b, c]));
    ds.clear_domain();
    assert!(ds.elements_of_domain().is_empty());
    assert!(ds.links_of_domain().is_empty());
}

#[test]
fn selector_collects_element_neighbors() {
    let mut ds = DelaunDataStructure::new(10);
    let a = ds.add_node(v(0.0, 0.0));
    let b = ds.add_node(v(2.0, 0.0));
    let c = ds.add_node(v(0.0, 2.0));
    let d = ds.add_node(v(2.0, 2.0));
    let ab = ds.add_link(a, b, VertexState::Free);
    let bc = ds.add_link(b, c, VertexState::Free);
    let _ca = ds.add_link(c, a, VertexState::Free);
    let cd = ds.add_link(c, d, VertexState::Free);
    let da = ds.add_link(d, a, VertexState::Free);
    let bd = ds.add_link(b, d, VertexState::Free);
    ds.add_element(DelaunTriangle::new([ab, bd, da], [a, b, d]));
    ds.add_element(DelaunTriangle::new([bc, cd, bd], [b, c, d]));
    let mut sel = DelaunSelector::new(&ds);
    sel.neighbors_of_node(a);
    // Triangle (b,c,d) does not touch node a; only the (a,b,d) triangle does.
    assert_eq!(sel.elements().len(), 1);
    let mut sel2 = DelaunSelector::new(&ds);
    sel2.neighbors_of_link(bd);
    assert_eq!(sel2.elements().len(), 2);
}
