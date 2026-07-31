//! Coordinate system conversions. Source: `Convert/`
//! Conversions between polar, cylindrical, spherical, and Cartesian coordinates.

/// Polar to Cartesian 2D. Returns (x, y).
pub fn polar_to_cartesian2d(rho: f64, theta: f64) -> (f64, f64) {
    (rho * theta.cos(), rho * theta.sin())
}

/// Cartesian 2D to polar. Returns (rho, theta).
pub fn cartesian2d_to_polar(x: f64, y: f64) -> (f64, f64) {
    ((x*x + y*y).sqrt(), y.atan2(x))
}

/// Cylindrical to Cartesian 3D. Returns (x, y, z).
pub fn cylindrical_to_cartesian(rho: f64, theta: f64, z: f64) -> (f64, f64, f64) {
    (rho * theta.cos(), rho * theta.sin(), z)
}

/// Cartesian 3D to cylindrical. Returns (rho, theta, z).
pub fn cartesian_to_cylindrical(x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    ((x*x + y*y).sqrt(), y.atan2(x), z)
}

/// Spherical to Cartesian 3D. Returns (x, y, z).
pub fn spherical_to_cartesian(r: f64, theta: f64, phi: f64) -> (f64, f64, f64) {
    let st = theta.sin(); let ct = theta.cos(); let sp = phi.sin(); let cp = phi.cos();
    (r * st * cp, r * st * sp, r * ct)
}

/// Cartesian 3D to spherical. Returns (r, theta, phi).
pub fn cartesian_to_spherical(x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    let r = (x*x + y*y + z*z).sqrt();
    let theta = if r > 0.0 { (z / r).acos() } else { 0.0 };
    let phi = y.atan2(x);
    (r, theta, phi)
}
