use occt_core::gp::{GpAx3, GpCone, GpDir, GpPnt};
use occt_geom::{GeomCone, Surface};

fn main() {
    let h = 1.0f64;
    let alpha1 = 0.4f64;
    let alpha2 = 0.8f64;
    let r_shared = h * alpha1.tan();
    let h2 = h + r_shared / alpha2.tan();
    let c1_ax = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
    let c1 = GpCone::new(c1_ax, 0.0, alpha1).unwrap();
    println!("cone1 apex={:?}", c1.apex());
    let c2_ax = GpAx3::new(GpPnt::new(0.0, 0.0, h2), GpDir::new(0.0, 0.0, -1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
    let c2 = GpCone::new(c2_ax, 0.0, alpha2).unwrap();
    println!("cone2 apex={:?}", c2.apex());
    let s1 = GeomCone::new(c1);
    let s2 = GeomCone::new(c2);
    let c1i = occt_topo::fillet_curved::cone_from_surface(&s1);
    let c2i = occt_topo::fillet_curved::cone_from_surface(&s2);
    println!("cone1 info={c1i:?}");
    println!("cone2 info={c2i:?}");
    let k1 = occt_topo::fillet_curved::classify_surface_analytic(&s1);
    let k2 = occt_topo::fillet_curved::classify_surface_analytic(&s2);
    println!("k1={k1:?} k2={k2:?}");
}
