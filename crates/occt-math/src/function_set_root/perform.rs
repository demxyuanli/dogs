//! `math_FunctionSetRoot::Perform` / `PerformWithBounds`
//! (`math_FunctionSetRoot.cxx:796-1417`).

use super::direction::{bounds, search_direction, search_direction_constrained, DirFunction};
use super::minimize::{minimize_direction_2pt, minimize_direction_3pt};
use super::*;

impl MathFunctionSetRoot {

    /// `Perform(F, StartingPoint, StopOnDivergent = false)`
    /// (`math_FunctionSetRoot.cxx:787-792`).
    pub fn perform<F: MathFunctionSetWithDerivatives>(
        &mut self,
        f: &mut F,
        starting_point: &MathVector,
        stop_on_divergent: bool,
    ) {
        let inf = self.inf_bound.clone();
        let sup = self.sup_bound.clone();
        self.perform_with_bounds(f, starting_point, &inf, &sup, stop_on_divergent);
    }

    /// `Perform(F, StartingPoint, InfBound, SupBound, StopOnDivergent = false)`
    /// (`math_FunctionSetRoot.cxx:796-1417`).
    pub fn perform_with_bounds<F: MathFunctionSetWithDerivatives>(
        &mut self,
        f: &mut F,
        starting_point: &MathVector,
        inf_bound: &MathVector,
        sup_bound: &MathVector,
        stop_on_divergent: bool,
    ) {
        let ninc = f.nb_variables();
        let neq = f.nb_equations();

        if neq == 0
            || starting_point.len() != ninc
            || inf_bound.len() != ninc
            || sup_bound.len() != ninc
        {
            // math_FunctionSetRoot.cxx:804-808 throws Standard_DimensionError.
            panic!("math_FunctionSetRoot::Perform: dimension error");
        }

        let epsk = 1.0e-16; // EpsSqrt
        let eps = 1.0e-32; // Eps
        let eps2 = 1.0e-64; // Eps2
        let progres = 0.005; // Progres

        let mut change_direction = false;
        let mut sort = false;
        let mut is_new_sol = false;

        // InvLengthMax / aConstraints (math_FunctionSetRoot.cxx:817-824).
        let mut inv_length_max = MathVector::new(1, ninc);
        for i in 1..=ninc {
            let sup_bound_val = sup_bound.value(i).min(PRECISION_INFINITE);
            let an_inf_bound = inf_bound.value(i).max(-PRECISION_INFINITE);
            inv_length_max.set_value(
                i,
                1.0 / ((sup_bound_val - an_inf_bound) / 4.0).max(1.0e-9),
            );
        }
        let mut constraints = vec![0i32; ninc + 1];

        let mut dir = DirFunction::new(ninc, neq);

        self.done = false;
        self.sol = starting_point.clone();
        self.kount = 0;

        self.is_divergent = false;
        for i in 1..=ninc {
            self.is_divergent = self.is_divergent
                || self.sol.value(i) < inf_bound.value(i)
                || self.sol.value(i) > sup_bound.value(i);
        }
        if stop_on_divergent && self.is_divergent {
            return;
        }

        // Recentrage sur les bornes (math_FunctionSetRoot.cxx:846-856).
        for i in 1..=ninc {
            if self.sol.value(i) <= inf_bound.value(i) {
                self.sol.set_value(i, inf_bound.value(i));
            } else if self.sol.value(i) > sup_bound.value(i) {
                self.sol.set_value(i, sup_bound.value(i));
            }
        }

        // Premiere valeur de F et de son gradient.
        let (mut f2, mut gnr1) = match dir.value_full(
            f,
            &self.sol,
            &mut self.ff,
            &mut self.df,
            &mut self.gh,
        ) {
            Some(v) => v,
            None => {
                self.done = false;
                if !stop_on_divergent || !self.is_divergent {
                    self.state = f.state_number();
                }
                return;
            }
        };
        let mut ambda2 = gnr1;
        self.save[0] = f2.max(epsk);
        let tol_func = epsilon_of(f2);

        if f2 <= eps || gnr1 <= eps2 {
            self.done = false;
            if !stop_on_divergent || !self.is_divergent {
                self.done = true;
                self.state = f.state_number();
            }
            return;
        }

        let mut dy = 0.0f64;
        self.kount = 1;
        while self.kount <= self.itermax {
            let previous_minimum = f2;
            let oldgr = gnr1;
            self.previous_solution = self.sol.clone();
            self.sol_save = self.sol.clone();

            // math_FunctionSetRoot.cxx:895 — ChangeDirection is passed by value
            // to SearchDirection and is not updated by it.
            search_direction(
                &self.df,
                &self.gh,
                &self.ff,
                change_direction,
                &inv_length_max,
                &mut self.dh,
                &mut dy,
            );
            if dy.abs() <= eps {
                self.done = false;
                if !stop_on_divergent || !self.is_divergent {
                    self.done = true;
                    // modified by jgv, 31.08.2011: update F before GetStateNumber
                    f.value(&self.sol, &mut self.ff);
                    self.state = f.state_number();
                }
                return;
            }
            let mut ambda;
            if change_direction {
                ambda = ambda2 / dy.abs().sqrt();
                if ambda > 1.0 {
                    ambda = 1.0;
                }
            } else {
                ambda = 1.0;
                ambda2 = 0.5 * ambda / self.dh.norm();
            }

            for i in 1..=ninc {
                self.sol
                    .set_value(i, self.sol.value(i) + ambda * self.dh.value(i));
            }
            for i in 1..=ninc {
                self.is_divergent = self.is_divergent
                    || self.sol.value(i) < inf_bound.value(i)
                    || self.sol.value(i) > sup_bound.value(i);
            }
            if stop_on_divergent && self.is_divergent {
                return;
            }

            let br = bounds(
                inf_bound,
                sup_bound,
                &self.tol,
                &mut self.sol,
                &self.sol_save,
                &mut constraints,
                &mut self.delta,
            );
            sort = br.0;
            is_new_sol = br.1;

            self.dhsave = self.gh.clone();
            if is_new_sol {
                match dir.value_full(
                    f,
                    &self.sol,
                    &mut self.ff,
                    &mut self.df,
                    &mut self.gh,
                ) {
                    Some(v) => {
                        f2 = v.0;
                        gnr1 = v.1;
                    }
                    None => {
                        self.done = false;
                        if !stop_on_divergent || !self.is_divergent {
                            self.state = f.state_number();
                        }
                        return;
                    }
                }
            }

            if f2 <= eps || gnr1 <= eps2 {
                self.done = false;
                if !stop_on_divergent || !self.is_divergent {
                    self.done = true;
                    f.value(&self.sol, &mut self.ff);
                    self.state = f.state_number();
                }
                return;
            }

            if sort || (f2 / previous_minimum > progres) {
                dy = self.gh.dot(&self.dh);
                let mut old_f = previous_minimum;
                let mut stop = false;
                let mut good = false;
                let mut descente_iter = 0;
                let mut sort_bis;

                // Standard processing without boundary handling.
                if !sort {
                    while (f2 / previous_minimum > progres) && !stop {
                        if f2 < old_f && dy < 0.0 {
                            descente_iter += 1;
                            self.sol_save = self.sol.clone();
                            old_f = f2;
                            for i in 1..=ninc {
                                self.sol.set_value(
                                    i,
                                    self.sol.value(i) + ambda * self.dh.value(i),
                                );
                            }
                            for i in 1..=ninc {
                                self.is_divergent = self.is_divergent
                                    || self.sol.value(i) < inf_bound.value(i)
                                    || self.sol.value(i) > sup_bound.value(i);
                            }
                            if stop_on_divergent && self.is_divergent {
                                return;
                            }
                            let br = bounds(
                                inf_bound,
                                sup_bound,
                                &self.tol,
                                &mut self.sol,
                                &self.sol_save,
                                &mut constraints,
                                &mut self.delta,
                            );
                            stop = br.0;
                            is_new_sol = br.1;
                            ambda *= 1.7;
                        } else {
                            if f2 >= old_f || f2 >= previous_minimum {
                                good = false;
                                if descente_iter == 0 {
                                    // C'est le premier pas qui flanche.
                                    descente_iter += 1;
                                    good = minimize_direction_2pt(
                                        &self.sol_save,
                                        &mut self.delta,
                                        old_f,
                                        f2,
                                        &self.dhsave,
                                        &self.gh,
                                        &self.tol,
                                        &mut dir,
                                        f,
                                    );
                                } else if change_direction
                                    || descente_iter > 1
                                    || old_f > previous_minimum
                                {
                                    descente_iter += 1;
                                    good = minimize_direction_3pt(
                                        &self.previous_solution,
                                        &self.sol_save,
                                        &self.sol,
                                        old_f,
                                        &mut self.delta,
                                        &self.tol,
                                        &mut dir,
                                        f,
                                    );
                                }
                                if !good {
                                    self.sol = self.sol_save.clone();
                                    f2 = old_f;
                                } else {
                                    self.sol = self.sol_save.added(&self.delta);
                                    for i in 1..=ninc {
                                        self.is_divergent = self.is_divergent
                                            || self.sol.value(i) < inf_bound.value(i)
                                            || self.sol.value(i) > sup_bound.value(i);
                                    }
                                    if stop_on_divergent && self.is_divergent {
                                        return;
                                    }
                                    let br = bounds(
                                        inf_bound,
                                        sup_bound,
                                        &self.tol,
                                        &mut self.sol,
                                        &self.sol_save,
                                        &mut constraints,
                                        &mut self.delta,
                                    );
                                    sort = br.0;
                                    is_new_sol = br.1;
                                }
                                sort = false; // On a rejete le point sur la frontiere
                            }
                            stop = true;
                        }
                        self.dhsave = self.gh.clone();
                        if is_new_sol {
                            match dir.value_full(
                                f,
                                &self.sol,
                                &mut self.ff,
                                &mut self.df,
                                &mut self.gh,
                            ) {
                                Some(v) => {
                                    f2 = v.0;
                                    gnr1 = v.1;
                                }
                                None => {
                                    self.done = false;
                                    if !stop_on_divergent || !self.is_divergent {
                                        self.state = f.state_number();
                                    }
                                    return;
                                }
                            }
                        }
                        dy = self.gh.dot(&self.dh);
                        if dy.abs() <= eps {
                            if f2 > old_f {
                                self.sol = self.sol_save.clone();
                            }
                            self.done = false;
                            if !stop_on_divergent || !self.is_divergent {
                                self.done = true;
                                f.value(&self.sol, &mut self.ff);
                                self.state = f.state_number();
                            }
                            return;
                        }
                        if descente_iter >= 100 {
                            stop = true;
                        }
                    }
                }
                // ------------------------------------
                //  on passe au traitement des bords
                // ------------------------------------
                if sort {
                    stop = f2 > 1.001 * old_f; // Pour ne pas progresser sur le bord
                    sort_bis = sort;
                    descente_iter = 0;
                    while sort_bis && (f2 < old_f || descente_iter == 0) && !stop {
                        descente_iter += 1;
                        // On essaye de progresser sur le bord.
                        self.sol_save = self.sol.clone();
                        old_f = f2;
                        search_direction_constrained(
                            &self.df,
                            &self.gh,
                            &self.ff,
                            &constraints,
                            &self.sol,
                            change_direction,
                            &inv_length_max,
                            &mut self.dh,
                            &mut dy,
                        );
                        if dy < -eps {
                            // Pour eviter des calculs inutiles et des /0...
                            if change_direction {
                                ambda = ambda2 / (-dy).sqrt();
                                if ambda > 1.0 {
                                    ambda = 1.0;
                                }
                            } else {
                                ambda = 1.0;
                                ambda2 = 0.5 * ambda / self.dh.norm();
                            }
                            for i in 1..=ninc {
                                self.sol.set_value(
                                    i,
                                    self.sol.value(i) + ambda * self.dh.value(i),
                                );
                            }
                            for i in 1..=ninc {
                                self.is_divergent = self.is_divergent
                                    || self.sol.value(i) < inf_bound.value(i)
                                    || self.sol.value(i) > sup_bound.value(i);
                            }
                            if stop_on_divergent && self.is_divergent {
                                return;
                            }
                            let br = bounds(
                                inf_bound,
                                sup_bound,
                                &self.tol,
                                &mut self.sol,
                                &self.sol_save,
                                &mut constraints,
                                &mut self.delta,
                            );
                            sort_bis = br.0;
                            is_new_sol = br.1;

                            self.dhsave = self.gh.clone();
                            if is_new_sol {
                                match dir.value_full(
                                    f,
                                    &self.sol,
                                    &mut self.ff,
                                    &mut self.df,
                                    &mut self.gh,
                                ) {
                                    Some(v) => {
                                        f2 = v.0;
                                        gnr1 = v.1;
                                    }
                                    None => {
                                        self.done = false;
                                        if !stop_on_divergent || !self.is_divergent {
                                            self.state = f.state_number();
                                        }
                                        return;
                                    }
                                }
                            }
                            ambda2 = gnr1;
                        } else {
                            stop = true;
                        }

                        while (f2 / previous_minimum > progres) && (f2 < old_f) && !stop {
                            descente_iter += 1;
                            if f2 < old_f && dy < 0.0 {
                                // On essaye de progresser dans cette direction.
                                self.sol_save = self.sol.clone();
                                old_f = f2;
                                for i in 1..=ninc {
                                    self.sol.set_value(
                                        i,
                                        self.sol.value(i) + ambda * self.dh.value(i),
                                    );
                                }
                                for i in 1..=ninc {
                                    self.is_divergent = self.is_divergent
                                        || self.sol.value(i) < inf_bound.value(i)
                                        || self.sol.value(i) > sup_bound.value(i);
                                }
                                if stop_on_divergent && self.is_divergent {
                                    return;
                                }
                                let br = bounds(
                                    inf_bound,
                                    sup_bound,
                                    &self.tol,
                                    &mut self.sol,
                                    &self.sol_save,
                                    &mut constraints,
                                    &mut self.delta,
                                );
                                sort_bis = br.0;
                                is_new_sol = br.1;
                            }
                            self.dhsave = self.gh.clone();
                            if is_new_sol {
                                match dir.value_full(
                                    f,
                                    &self.sol,
                                    &mut self.ff,
                                    &mut self.df,
                                    &mut self.gh,
                                ) {
                                    Some(v) => {
                                        f2 = v.0;
                                        gnr1 = v.1;
                                    }
                                    None => {
                                        self.done = false;
                                        if !stop_on_divergent || !self.is_divergent {
                                            self.state = f.state_number();
                                        }
                                        return;
                                    }
                                }
                            }
                            ambda2 = gnr1;
                            dy = self.gh.dot(&self.dh);
                            stop = dy >= 0.0 || descente_iter >= 10 || sort_bis;
                        }
                        stop = dy >= 0.0 || descente_iter >= 10;
                    }
                    if ((f2 / previous_minimum > progres) && (f2 >= old_f))
                        || (f2 >= previous_minimum)
                    {
                        // On minimise par Brent.
                        descente_iter += 1;
                        good = minimize_direction_2pt(
                            &self.sol_save,
                            &mut self.delta,
                            old_f,
                            f2,
                            &self.dhsave,
                            &self.gh,
                            &self.tol,
                            &mut dir,
                            f,
                        );
                        if !good {
                            self.sol = self.sol_save.clone();
                            sort = false;
                        } else {
                            self.sol = self.sol_save.added(&self.delta);
                            for i in 1..=ninc {
                                self.is_divergent = self.is_divergent
                                    || self.sol.value(i) < inf_bound.value(i)
                                    || self.sol.value(i) > sup_bound.value(i);
                            }
                            if stop_on_divergent && self.is_divergent {
                                return;
                            }
                            let br = bounds(
                                inf_bound,
                                sup_bound,
                                &self.tol,
                                &mut self.sol,
                                &self.sol_save,
                                &mut constraints,
                                &mut self.delta,
                            );
                            sort = br.0;
                            is_new_sol = br.1;
                            if is_new_sol {
                                match dir.value_full(
                                    f,
                                    &self.sol,
                                    &mut self.ff,
                                    &mut self.df,
                                    &mut self.gh,
                                ) {
                                    Some(v) => {
                                        f2 = v.0;
                                        gnr1 = v.1;
                                    }
                                    None => {
                                        self.done = false;
                                        if !stop_on_divergent || !self.is_divergent {
                                            self.state = f.state_number();
                                        }
                                        return;
                                    }
                                }
                            }
                        }
                        dy = self.gh.dot(&self.dh);
                    }
                }

            }
            // ---------------------------------------------
            //  on passe aux tests d'ARRET
            // ---------------------------------------------
            self.save[self.kount as usize] = f2;
            let verif: bool;
            if change_direction {
                verif = true;
                // Gradient : Il faut eviter de boucler.
            } else {
                if self.kount > 1 {
                    // Pour accelerer les cas quasi-quadratique.
                    verif = self.save[(self.kount - 1) as usize]
                        < 1.0e-4 * self.save[(self.kount - 2) as usize];
                } else {
                    // Pour les cas dejas solutions.
                    verif = f2 < 1.0e-6 * self.save[0];
                }
            }
            if verif {
                for i in 1..=ninc {
                    self.delta.set_value(
                        i,
                        self.previous_solution.value(i) - self.sol.value(i),
                    );
                }
                if self.is_solution_reached() {
                    if previous_minimum < f2 {
                        self.sol = self.sol_save.clone();
                    }
                    self.done = false;
                    if !stop_on_divergent || !self.is_divergent {
                        self.done = true;
                        f.value(&self.sol, &mut self.ff);
                        self.state = f.state_number();
                    }
                    return;
                }
            }

            // Analyse de la progression...
            if (f2 - previous_minimum) <= tol_func {
                if self.kount > 5 {
                    // L'historique est il bon ?
                    if f2 >= 0.95 * self.save[(self.kount - 5) as usize] {
                        if !change_direction {
                            change_direction = true;
                        } else {
                            self.done = false;
                            if !stop_on_divergent || !self.is_divergent {
                                self.done = true;
                                self.state = f.state_number();
                            }
                            // si un gain inf a 5% on sort
                            return;
                        }
                    } else {
                        change_direction = false; // If yes we restart
                    }
                } else {
                    change_direction = false; // No history, we continue
                }
                // If the gradient does not decrease sufficiently with Newton, we
                // try the gradient method unless f decreases.
                if (gnr1 > 0.9 * oldgr) && (f2 > 0.5 * previous_minimum) {
                    change_direction = true;
                }

                if !change_direction && !verif {
                    for i in 1..=ninc {
                        self.delta.set_value(
                            i,
                            self.previous_solution.value(i) - self.sol.value(i),
                        );
                    }
                    if self.is_solution_reached() {
                        self.done = false;
                        if !stop_on_divergent || !self.is_divergent {
                            self.done = true;
                            f.value(&self.sol, &mut self.ff);
                            self.state = f.state_number();
                        }
                        return;
                    }
                }
            } else {
                // Cas de regression.
                if !change_direction {
                    // On passe au gradient.
                    change_direction = true;
                    self.sol = self.previous_solution.clone();
                    match dir.value_full(
                        f,
                        &self.sol,
                        &mut self.ff,
                        &mut self.df,
                        &mut self.gh,
                    ) {
                        Some(v) => {
                            f2 = v.0;
                            gnr1 = v.1;
                        }
                        None => {
                            self.done = false;
                            if !stop_on_divergent || !self.is_divergent {
                                self.state = f.state_number();
                            }
                            return;
                        }
                    }
                } else {
                    if !stop_on_divergent || !self.is_divergent {
                        self.state = f.state_number();
                    }
                    // y a plus d'issues
                    return;
                }
            }

            self.kount += 1;
        }
        if !stop_on_divergent || !self.is_divergent {
            self.state = f.state_number();
        }
    }
}
