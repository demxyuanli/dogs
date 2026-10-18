//! TEMPORARY probe (T-100): does `ginter`'s generic arm actually handle a
//! Line/Conic pair? `ginter.rs` claims the dedicated `IntConicConic` overloads
//! are an optimisation, not new capability. This measures it end-to-end through
//! the public dispatcher -- the exact code path in question.
//!
//! Line: the x axis. Circle: centre (0,0), radius 5 => expect x = +/-5.
//! Before T-100 this pair fell into the empty UNPORTED arm, which returned no
//! intersection at all.
use occt_core::gp::{GpAx22d, GpCirc2d, GpDir2d, GpPnt2d};
use occt_core::intres2d::IntRes2dDomain;
use occt_geom2d::geom2d_int::Geom2dIntGInter;
use occt_geom2d::{Curve2d, Geom2dCircle, Geom2dLine};

fn main() {
    let line = Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, 0.0), GpDir2d::new(1.0, 0.0).unwrap());
    let circ = Geom2dCircle::new(GpCirc2d::new(GpAx22d::standard(), 5.0));

    // Domain of the circle: one full turn, given by its two end points.
    let p0 = circ.d0(0.0);
    let p2pi = circ.d0(2.0 * std::f64::consts::PI);
    let d2 = IntRes2dDomain::bounded(&p0, 0.0, 1e-9, &p2pi, 2.0 * std::f64::consts::PI, 1e-9);

    let mut inter = Geom2dIntGInter::new();
    inter.perform_with_d2(&line, &circ, &d2, 1e-9, 1e-9);

    println!("is_done={}", inter.is_done());
    let n = inter.nb_points();
    println!("nb_points={n}");
    let mut xs = Vec::new();
    for i in 1..=n {
        let p = inter.point(i);
        let t = p.param_on_second();
        let xy = circ.d0(t);
        println!("  point[{i}] t2={t:.9} xy=({:.9},{:.9})", xy.x(), xy.y());
        xs.push(xy.x());
    }
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    xs.dedup_by(|a, b| (*a - *b).abs() < 1e-6);
    println!("x values (deduped): {xs:?}");
    // The domain is the closed span [0, 2pi], so the point at t = 2pi coincides
    // with the one at t = 0 and appears twice; dedup above mirrors that.
    let ok = xs.len() == 2 && (xs[0] + 5.0).abs() < 1e-6 && (xs[1] - 5.0).abs() < 1e-6;
    println!("EXPECT +/-5 => {}", if ok { "PASS" } else { "FAIL" });
}
