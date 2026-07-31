//! Map projection utilities: Mercator, equirectangular and Lambert azimuthal
//! equal-area projections plus a UTM zone helper.

use std::f64::consts::PI;

fn to_rad(deg: f64) -> f64 {
    deg * PI / 180.0
}

fn to_deg(rad: f64) -> f64 {
    rad * 180.0 / PI
}

/// Web-style Mercator projection: `x = R * lon`, `y = R * ln(tan(pi/4 + lat/2))`.
pub fn lonlat_to_mercator(lon_deg: f64, lat_deg: f64, radius: f64) -> (f64, f64) {
    let lon = to_rad(lon_deg);
    let lat = to_rad(lat_deg);
    let x = radius * lon;
    let y = radius * (PI / 4.0 + lat / 2.0).tan().ln();
    (x, y)
}

/// Inverse of [`lonlat_to_mercator`].
pub fn mercator_to_lonlat(x: f64, y: f64, radius: f64) -> (f64, f64) {
    let lon = x / radius;
    let lat = 2.0 * (y / radius).exp().atan() - PI / 2.0;
    (to_deg(lon), to_deg(lat))
}

/// Equirectangular (plate carrée) projection about `ref_lon`. The reference
/// latitude is fixed at 0, so `x = R * (lon - ref_lon)`, `y = R * lat`.
pub fn lonlat_to_equirectangular(
    lon_deg: f64,
    lat_deg: f64,
    radius: f64,
    ref_lon_deg: f64,
) -> (f64, f64) {
    let x = radius * to_rad(lon_deg - ref_lon_deg);
    let y = radius * to_rad(lat_deg);
    (x, y)
}

/// Lambert azimuthal equal-area projection, north polar aspect:
/// `k = 2 R sin((pi/2 - lat)/2)`, `theta = lon`, `x = k cos(theta)`,
/// `y = -k sin(theta)`.
pub fn lonlat_to_lambert_azimuthal(lon_deg: f64, lat_deg: f64, radius: f64) -> (f64, f64) {
    let lon = to_rad(lon_deg);
    let lat = to_rad(lat_deg);
    let k = 2.0 * radius * ((PI / 2.0 - lat) / 2.0).sin();
    let x = k * lon.cos();
    let y = -k * lon.sin();
    (x, y)
}

/// UTM longitude zone (`1..=60`) for a longitude in degrees.
pub fn utm_zone(lon_deg: f64) -> i32 {
    let zone = ((lon_deg + 180.0) / 6.0).floor() as i32 + 1;
    zone.clamp(1, 60)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_approx(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() < tol, "{a} != {b}");
    }

    #[test]
    fn mercator_roundtrip() {
        let (x, y) = lonlat_to_mercator(45.0, 30.0, 6378137.0);
        let (lon, lat) = mercator_to_lonlat(x, y, 6378137.0);
        assert_approx(lon, 45.0, 1e-9);
        assert_approx(lat, 30.0, 1e-9);
    }

    #[test]
    fn utm_zone_edges() {
        assert_eq!(utm_zone(-180.0), 1);
        assert_eq!(utm_zone(179.0), 60);
    }
}
