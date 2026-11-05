//! `BRepFill_PipeShell::Simulate` (`BRepFill_PipeShell.cxx:651-701`).

use crate::brep_fill_pipe_shell::PipeShell;
use crate::shape::TopoShape;

impl PipeShell {
    /// `Simulate(N, List)` (`cxx:651-701`): `N` sections spread by arc length
    /// along the path.
    ///
    /// `First` / `Last` keep the values of the previous iteration in OCCT, so
    /// they are carried here the same way rather than recomputed per `U`.
    pub fn simulate(&mut self, n: usize) -> Result<Vec<TopoShape>, String> {
        // Preparation.
        self.prepare()?;

        let section = self
            .section()
            .cloned()
            .ok_or_else(|| "BRepFill_PipeShell::Simulate: no section law".to_string())?;
        let location = self
            .location()
            .cloned()
            .ok_or_else(|| "BRepFill_PipeShell::Simulate: no location law".to_string())?;

        let nb_l = location.base().nb_law();
        let (first_s, last) = section
            .base().law(1)
            .ok_or_else(|| "BRepFill_PipeShell::Simulate: no section law 1".to_string())?
            .get_domain();
        let delta_s = last - first_s;

        // `cxx:665-671`: the parameters of digitalization.
        let (_length_f, length) = location.base().curvilinear_bounds(nb_l)?;
        let mut delta = length;
        if n > 1 {
            delta /= (n - 1) as f64;
        }

        // `cxx:676`: `CurvilinearBounds(1, First, Last)` initiates `Last`.
        let (mut first_edge, mut last_edge) = location.base().curvilinear_bounds(1)?;

        let mut list: Vec<TopoShape> = Vec::new();
        let mut u = 0.0_f64;
        let mut ii = 1usize;
        let mut finis = false;
        while !finis {
            if u >= length {
                u = length;
                finis = true;
            } else {
                if ii < nb_l {
                    let (first, last) = location.base().curvilinear_bounds(nb_l)?;
                    first_edge = first;
                    last_edge = last;
                }
                if u > last_edge {
                    u = (last_edge + first_edge) / 2.0; // The edge is not skipped
                }
                if u > first_edge {
                    ii += 1;
                }
            }

            let us = first_s + (u / length) * delta_s;
            let w = section
                .d0(us)
                .ok_or_else(|| "BRepFill_PipeShell::Simulate: no section shape".to_string())?;
            let moved = location
                .base()
                .d0(u, &w)?
                .ok_or_else(|| "BRepFill_PipeShell::Simulate: section off the path".to_string())?;
            list.push(moved);

            u += delta;
        }

        Ok(list)
    }
}

/// `Simulate` in free-function form (`cxx:651-701`).
pub fn simulate(pipe: &mut PipeShell, n: usize) -> Result<Vec<TopoShape>, String> {
    pipe.simulate(n)
}
