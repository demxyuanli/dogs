//! 2-D interpolation on a regular grid: bilinear and nearest-neighbor.
//! Source: `math_Interpolation` (2-D), bilinear form.

/// A regular grid of samples. `values[i][j]` holds the sample at
/// `(xs[i], ys[j])`.
pub struct Grid2d {
    pub xs: Vec<f64>,
    pub ys: Vec<f64>,
    pub values: Vec<Vec<f64>>,
}

impl Grid2d {
    /// Build a grid, validating strictly increasing axes and value dimensions.
    pub fn new(xs: Vec<f64>, ys: Vec<f64>, values: Vec<Vec<f64>>) -> Result<Self, String> {
        if xs.len() < 2 || ys.len() < 2 {
            return Err("Grid2d: need at least 2 samples along each axis".to_string());
        }
        if xs.windows(2).any(|w| w[1] <= w[0]) {
            return Err("Grid2d: xs must be strictly increasing".to_string());
        }
        if ys.windows(2).any(|w| w[1] <= w[0]) {
            return Err("Grid2d: ys must be strictly increasing".to_string());
        }
        if values.len() != xs.len() {
            return Err("Grid2d: values row count must equal xs length".to_string());
        }
        if values.iter().any(|row| row.len() != ys.len()) {
            return Err("Grid2d: every values row must equal ys length".to_string());
        }
        Ok(Self { xs, ys, values })
    }

    pub fn nx(&self) -> usize {
        self.xs.len()
    }

    pub fn ny(&self) -> usize {
        self.ys.len()
    }

    /// Sample value at grid node `(i, j)`.
    pub fn at(&self, i: usize, j: usize) -> f64 {
        self.values[i][j]
    }
}

/// Index of the cell `[axis[i], axis[i+1]]` containing `v` (assumed within
/// range). The upper grid line maps into the last cell.
fn cell_index(axis: &[f64], v: f64) -> usize {
    if v == axis[axis.len() - 1] {
        return axis.len() - 2;
    }
    let mut i = 0;
    while i + 1 < axis.len() && v > axis[i + 1] {
        i += 1;
    }
    i
}

/// Index of the axis sample nearest to `v`.
fn nearest_index(axis: &[f64], v: f64) -> usize {
    let mut best = 0;
    let mut best_d = (v - axis[0]).abs();
    for (k, &a) in axis.iter().enumerate().skip(1) {
        let d = (v - a).abs();
        if d < best_d {
            best_d = d;
            best = k;
        }
    }
    best
}

/// Bilinear interpolation on the cell containing `(x, y)`. Errors when the
/// point is outside the grid's bounding box. Exact grid-line hits are exact.
pub fn bilinear_interp(g: &Grid2d, x: f64, y: f64) -> Result<f64, String> {
    let nx = g.xs.len();
    let ny = g.ys.len();
    if x < g.xs[0] || x > g.xs[nx - 1] || y < g.ys[0] || y > g.ys[ny - 1] {
        return Err("bilinear_interp: point out of grid range".to_string());
    }
    let i = cell_index(&g.xs, x);
    let j = cell_index(&g.ys, y);
    let tx = (x - g.xs[i]) / (g.xs[i + 1] - g.xs[i]);
    let ty = (y - g.ys[j]) / (g.ys[j + 1] - g.ys[j]);
    let (v00, v10, v01, v11) = (
        g.values[i][j],
        g.values[i + 1][j],
        g.values[i][j + 1],
        g.values[i + 1][j + 1],
    );
    Ok(v00 * (1.0 - tx) * (1.0 - ty)
        + v10 * tx * (1.0 - ty)
        + v01 * (1.0 - tx) * ty
        + v11 * tx * ty)
}

/// Nearest-sample interpolation (closest grid node in both axes).
pub fn nearest_interp2d(g: &Grid2d, x: f64, y: f64) -> f64 {
    let i = nearest_index(&g.xs, x);
    let j = nearest_index(&g.ys, y);
    g.values[i][j]
}

/// Gradient `(∂f/∂x, ∂f/∂y)` of the bilinear surface at `(x, y)`. Errors when
/// the point is outside the grid.
pub fn bilinear_gradient(g: &Grid2d, x: f64, y: f64) -> Result<(f64, f64), String> {
    let nx = g.xs.len();
    let ny = g.ys.len();
    if x < g.xs[0] || x > g.xs[nx - 1] || y < g.ys[0] || y > g.ys[ny - 1] {
        return Err("bilinear_gradient: point out of grid range".to_string());
    }
    let i = cell_index(&g.xs, x);
    let j = cell_index(&g.ys, y);
    let tx = (x - g.xs[i]) / (g.xs[i + 1] - g.xs[i]);
    let ty = (y - g.ys[j]) / (g.ys[j + 1] - g.ys[j]);
    let (v00, v10, v01, v11) = (
        g.values[i][j],
        g.values[i + 1][j],
        g.values[i][j + 1],
        g.values[i + 1][j + 1],
    );
    let dx = g.xs[i + 1] - g.xs[i];
    let dy = g.ys[j + 1] - g.ys[j];
    let dfdx = ((1.0 - ty) * (v10 - v00) + ty * (v11 - v01)) / dx;
    let dfdy = ((1.0 - tx) * (v01 - v00) + tx * (v11 - v10)) / dy;
    Ok((dfdx, dfdy))
}

/// Sample `f` on a uniform `nx × ny` grid spanning `[x0, x1] × [y0, y1]`.
/// Axes include both endpoints; `nx`/`ny` below 2 are clamped to 2.
pub fn build_grid_uniform(
    x0: f64,
    x1: f64,
    nx: usize,
    y0: f64,
    y1: f64,
    ny: usize,
    f: &dyn Fn(f64, f64) -> f64,
) -> Grid2d {
    let nx = nx.max(2);
    let ny = ny.max(2);
    let xs: Vec<f64> = (0..nx).map(|i| x0 + i as f64 * (x1 - x0) / (nx - 1) as f64).collect();
    let ys: Vec<f64> = (0..ny).map(|j| y0 + j as f64 * (y1 - y0) / (ny - 1) as f64).collect();
    let values: Vec<Vec<f64>> = xs
        .iter()
        .map(|&x| ys.iter().map(|&y| f(x, y)).collect())
        .collect();
    Grid2d { xs, ys, values }
}

/// Maximum `|bilinear_interp − f|` over a `samples × samples` grid of query
/// points spanning the grid's bounding box.
pub fn grid_max_error(g: &Grid2d, f: &dyn Fn(f64, f64) -> f64, samples: usize) -> f64 {
    if samples == 0 {
        return 0.0;
    }
    let (x0, x1) = (g.xs[0], g.xs[g.xs.len() - 1]);
    let (y0, y1) = (g.ys[0], g.ys[g.ys.len() - 1]);
    let mut max = 0.0f64;
    for i in 0..samples {
        let tx = if samples == 1 { 0.0 } else { i as f64 / (samples - 1) as f64 };
        let x = x0 + tx * (x1 - x0);
        for j in 0..samples {
            let ty = if samples == 1 { 0.0 } else { j as f64 / (samples - 1) as f64 };
            let y = y0 + ty * (y1 - y0);
            let approx = bilinear_interp(g, x, y).unwrap_or(0.0);
            let err = (approx - f(x, y)).abs();
            if err > max {
                max = err;
            }
        }
    }
    max
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn linear_grid() -> Grid2d {
        // f(x, y) = x + 2y on xs = ys = [0, 1, 2].
        let xs = vec![0.0, 1.0, 2.0];
        let ys = vec![0.0, 1.0, 2.0];
        let values: Vec<Vec<f64>> = xs
            .iter()
            .map(|&x| ys.iter().map(|&y| x + 2.0 * y).collect())
            .collect();
        Grid2d::new(xs, ys, values).unwrap()
    }

    #[test]
    fn bilinear_reproduces_linear_function() {
        let g = linear_grid();
        for (x, y) in [(0.5, 0.5), (1.3, 0.7), (0.0, 1.0), (2.0, 2.0), (1.0, 2.0)] {
            let v = bilinear_interp(&g, x, y).unwrap();
            assert!((v - (x + 2.0 * y)).abs() < 1e-9, "at ({x}, {y}) got {v}");
        }
        // Exact grid-node hits.
        assert!((bilinear_interp(&g, 1.0, 1.0).unwrap() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn gradient_of_quadratic() {
        // f(x, y) = x² + y; bilinear gradient gives ∂f/∂x = 2x on cell centers.
        let xs = vec![0.0, 1.0, 2.0, 3.0];
        let ys = vec![0.0, 1.0, 2.0];
        let values: Vec<Vec<f64>> = xs
            .iter()
            .map(|&x| ys.iter().map(|&y| x * x + y).collect())
            .collect();
        let g = Grid2d::new(xs, ys, values).unwrap();
        let (dfdx, dfdy) = bilinear_gradient(&g, 1.5, 1.0).unwrap();
        assert!((dfdx - 3.0).abs() < 1e-12, "dfdx = {dfdx}");
        assert!((dfdy - 1.0).abs() < 1e-12, "dfdy = {dfdy}");
    }

    #[test]
    fn build_grid_uniform_approximates_smooth_function() {
        let g = build_grid_uniform(0.0, 2.0 * PI, 41, 0.0, 2.0 * PI, 41, &|x, y| {
            x.sin() * y.cos()
        });
        let err = grid_max_error(&g, &|x, y| x.sin() * y.cos(), 60);
        assert!(err < 0.02, "max error = {err}");
    }

    #[test]
    fn out_of_range_errors() {
        let g = linear_grid();
        assert!(bilinear_interp(&g, -0.1, 0.5).is_err());
        assert!(bilinear_interp(&g, 0.5, 2.5).is_err());
        assert!(bilinear_gradient(&g, 3.0, 1.0).is_err());
    }

    #[test]
    fn nearest_picks_closest_sample() {
        let xs = vec![0.0, 1.0, 10.0];
        let ys = vec![0.0, 1.0, 10.0];
        let values: Vec<Vec<f64>> = xs
            .iter()
            .map(|&x| ys.iter().map(|&y| x + 2.0 * y).collect())
            .collect();
        let g = Grid2d::new(xs, ys, values).unwrap();
        assert_eq!(nearest_interp2d(&g, 0.4, 0.3), 0.0); // (0, 0)
        assert_eq!(nearest_interp2d(&g, 7.0, 7.0), 30.0); // (10, 10)
        assert_eq!(nearest_interp2d(&g, 7.0, 0.2), 10.0); // (10, 0)
    }
}
