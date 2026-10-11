//! `AdvApp2Var_Network` for `InitGrid(1)` and the cuts in `ComputePatches`
//! (`AdvApp2Var_Network.cxx:31-130`, `AdvApp2Var_ApproxAFunc2Var.cxx:239-255`).

struct NetPatch {
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    ord_u: i32,
    ord_v: i32,
    approximated: bool,
    equation: Vec<f64>,
    nb_coeff_u: i32,
    nb_coeff_v: i32,
    max_error: f64,
    moy_error: f64,
    iso_errors: [f64; 4],
}

pub(super) struct Network {
    patches: Vec<NetPatch>,
    u_params: Vec<f64>,
    v_params: Vec<f64>,
}

impl Network {
    pub(super) fn init(u0: f64, u1: f64, v0: f64, v1: f64, iu: i32, iv: i32) -> Self {
        Self {
            patches: vec![NetPatch {
                u0,
                u1,
                v0,
                v1,
                ord_u: iu,
                ord_v: iv,
                approximated: false,
                equation: Vec::new(),
                nb_coeff_u: 0,
                nb_coeff_v: 0,
                max_error: 0.0,
                moy_error: 0.0,
                iso_errors: [0.0; 4],
            }],
            u_params: vec![u0, u1],
            v_params: vec![v0, v1],
        }
    }

    pub(super) fn nb_patch_in_u(&self) -> i32 {
        self.u_params.len() as i32 - 1
    }

    pub(super) fn nb_patch_in_v(&self) -> i32 {
        self.v_params.len() as i32 - 1
    }

    /// `FirstNotApprox` (`AdvApp2Var_Network.cxx:43-58`). Index is 1-based.
    pub(super) fn first_not_approx(&self) -> Option<usize> {
        self.patches.iter().position(|patch| !patch.approximated)
    }

    pub(super) fn domain(&self, index: usize) -> Result<(f64, f64, f64, f64), String> {
        self.patches.get(index).map(|patch| (patch.u0, patch.u1, patch.v0, patch.v1)).ok_or_else(
            || "AdvApp2Var_Network::FirstNotApprox: AdvApp2Var_Network.cxx:50".to_string(),
        )
    }

    /// `UpdateInU` (`AdvApp2Var_Network.cxx:62-91`).
    pub(super) fn update_in_u(&mut self, cutting: f64) -> Result<(), String> {
        let mut i = 1i32;
        while (i as usize) <= self.u_params.len() && self.u_params[(i - 1) as usize] < cutting {
            i += 1;
        }
        if i < 1 || i as usize > self.u_params.len() {
            return Err("AdvApp2Var_Network::UpdateInU: AdvApp2Var_Network.cxx:66".into());
        }
        self.u_params.insert((i - 1) as usize, cutting);
        let nb_u = self.u_params.len() as i32 - 1;
        let mut j = 1i32;
        while j < self.v_params.len() as i32 {
            let patch_index = nb_u * (j - 1) + i - 1;
            let at = (patch_index - 1) as usize;
            let patch = self
                .patches
                .get(at)
                .ok_or_else(|| "AdvApp2Var_Network::UpdateInU: AdvApp2Var_Network.cxx:77".to_string())?;
            let (ord_u, ord_v) = (patch.ord_u, patch.ord_v);
            self.patches[at].u1 = cutting;
            self.patches[at].approximated = false;
            let u_next = *self.u_params.get(i as usize).ok_or_else(|| {
                "AdvApp2Var_Network::UpdateInU: AdvApp2Var_Network.cxx:83".to_string()
            })?;
            let v_j = *self.v_params.get((j - 1) as usize).ok_or_else(|| {
                "AdvApp2Var_Network::UpdateInU: AdvApp2Var_Network.cxx:84".to_string()
            })?;
            let v_next = *self.v_params.get(j as usize).ok_or_else(|| {
                "AdvApp2Var_Network::UpdateInU: AdvApp2Var_Network.cxx:85".to_string()
            })?;
            self.patches.insert(
                at + 1,
                NetPatch {
                    u0: cutting,
                    u1: u_next,
                    v0: v_j,
                    v1: v_next,
                    ord_u,
                    ord_v,
                    approximated: false,
                equation: Vec::new(),
                nb_coeff_u: 0,
                nb_coeff_v: 0,
                max_error: 0.0,
                moy_error: 0.0,
                iso_errors: [0.0; 4],
                },
            );
            j += 1;
        }
        Ok(())
    }

    /// `UpdateInV` (`AdvApp2Var_Network.cxx:95-130`).
    pub(super) fn update_in_v(&mut self, cutting: f64) -> Result<(), String> {
        let mut j = 1i32;
        while (j as usize) <= self.v_params.len() && self.v_params[(j - 1) as usize] < cutting {
            j += 1;
        }
        if j < 1 || j as usize > self.v_params.len() {
            return Err("AdvApp2Var_Network::UpdateInV: AdvApp2Var_Network.cxx:99".into());
        }
        self.v_params.insert((j - 1) as usize, cutting);
        let nb_u = self.u_params.len() as i32 - 1;
        for i in 1..=nb_u {
            let at = (nb_u * (j - 2) + i - 1) as usize;
            let patch = self.patches.get_mut(at).ok_or_else(|| {
                "AdvApp2Var_Network::UpdateInV: AdvApp2Var_Network.cxx:111".to_string()
            })?;
            patch.v1 = cutting;
            patch.approximated = false;
        }
        for i in 1..=nb_u {
            let source_at = (nb_u * (j - 2) + i - 1) as usize;
            let source = self.patches.get(source_at).ok_or_else(|| {
                "AdvApp2Var_Network::UpdateInV: AdvApp2Var_Network.cxx:120".to_string()
            })?;
            let (ord_u, ord_v) = (source.ord_u, source.ord_v);
            let u_i = *self.u_params.get((i - 1) as usize).ok_or_else(|| {
                "AdvApp2Var_Network::UpdateInV: AdvApp2Var_Network.cxx:121".to_string()
            })?;
            let u_next = *self.u_params.get(i as usize).ok_or_else(|| {
                "AdvApp2Var_Network::UpdateInV: AdvApp2Var_Network.cxx:122".to_string()
            })?;
            let v_next = *self.v_params.get(j as usize).ok_or_else(|| {
                "AdvApp2Var_Network::UpdateInV: AdvApp2Var_Network.cxx:124".to_string()
            })?;
            let insert_at = (nb_u * (j - 1) + i - 1) as usize;
            if insert_at > self.patches.len() {
                return Err("AdvApp2Var_Network::UpdateInV: AdvApp2Var_Network.cxx:128".into());
            }
            self.patches.insert(
                insert_at,
                NetPatch {
                    u0: u_i,
                    u1: u_next,
                    v0: cutting,
                    v1: v_next,
                    ord_u,
                    ord_v,
                    approximated: false,
                equation: Vec::new(),
                nb_coeff_u: 0,
                nb_coeff_v: 0,
                max_error: 0.0,
                moy_error: 0.0,
                iso_errors: [0.0; 4],
                },
            );
        }
        Ok(())
    }

    /// `myApprIsDone` after `MakeApprox` (`AdvApp2Var_Patch.cxx:887` and `:957`).
    pub(super) fn store_result(
        &mut self,
        index: usize,
        done: super::degree_error::FinishedPatch,
    ) -> Result<(), String> {
        let patch = self.patches.get_mut(index).ok_or_else(|| {
            "AdvApp2Var_Network::FirstNotApprox: AdvApp2Var_Network.cxx:50".to_string()
        })?;
        patch.approximated = true;
        patch.equation = done.equation;
        patch.nb_coeff_u = done.nb_coeff_u;
        patch.nb_coeff_v = done.nb_coeff_v;
        patch.max_error = done.max_error;
        patch.moy_error = done.moy_error;
        patch.iso_errors = done.iso_errors;
        Ok(())
    }

    /// `Compute3DErrors` (`AdvApp2Var_ApproxAFunc2Var.cxx:818-863`).
    /// `Perform` has already set `myDone` (`:341`). Returns that flag.
    pub(super) fn compute_3d_errors(
        &self,
        nb3: i32,
        tol3: &[f64],
        front3: &[[f64; 4]],
    ) -> Result<bool, String> {
        if nb3 < 1 {
            return Ok(true);
        }
        let npatch = self.patches.len();
        if npatch == 0 {
            return Err(
                "AdvApp2Var_ApproxAFunc2Var::Compute3DErrors: AdvApp2Var_ApproxAFunc2Var.cxx:855"
                    .into(),
            );
        }
        let mut done = true;
        for iesp in 1..=nb3 {
            if iesp != 1 {
                return Err(
                    "AdvApp2Var_ApproxAFunc2Var::Compute3DErrors: AdvApp2Var_ApproxAFunc2Var.cxx:845"
                        .into(),
                );
            }
            let mut error_max: f64 = 0.0;
            let mut error_moy: f64 = 0.0;
            let mut error_u0: f64 = 0.0;
            let mut error_v0: f64 = 0.0;
            let mut error_u1: f64 = 0.0;
            let mut error_v1: f64 = 0.0;
            let tol = *tol3.get((iesp - 1) as usize).ok_or_else(|| {
                "AdvApp2Var_ApproxAFunc2Var::Compute3DErrors: AdvApp2Var_ApproxAFunc2Var.cxx:838"
                    .to_string()
            })?;
            let front = front3.get((iesp - 1) as usize).ok_or_else(|| {
                "AdvApp2Var_ApproxAFunc2Var::Compute3DErrors: AdvApp2Var_ApproxAFunc2Var.cxx:839"
                    .to_string()
            })?;
            for patch in &self.patches {
                error_max = error_max.max(patch.max_error);
                error_u0 = error_u0.max(patch.iso_errors[2]);
                error_u1 = error_u1.max(patch.iso_errors[3]);
                error_v0 = error_v0.max(patch.iso_errors[0]);
                error_v1 = error_v1.max(patch.iso_errors[1]);
                error_moy += patch.moy_error;
            }
            let _u_front = error_u0.max(error_u1);
            let _v_front = error_v0.max(error_v1);
            error_moy /= npatch as f64;
            let _average = error_moy;
            if error_max > tol
                || error_u0 > front[2]
                || error_u1 > front[3]
                || error_v0 > front[0]
                || error_v1 > front[1]
            {
                done = false;
            }
        }
        Ok(done)
    }

    /// `SameDegree` (`AdvApp2Var_Network.cxx:135-156`). Counts only rise.
    pub(super) fn same_degree(&mut self, iu: i32, iv: i32) {
        let mut ncfu = 2 * iu + 2;
        let mut ncfv = 2 * iv + 2;
        for patch in &self.patches {
            ncfu = ncfu.max(patch.nb_coeff_u);
            ncfv = ncfv.max(patch.nb_coeff_v);
        }
        for patch in &mut self.patches {
            if patch.nb_coeff_u < ncfu {
                patch.nb_coeff_u = ncfu;
            }
            if patch.nb_coeff_v < ncfv {
                patch.nb_coeff_v = ncfv;
            }
        }
    }

    pub(super) fn u_parameters(&self) -> &[f64] {
        &self.u_params
    }

    pub(super) fn v_parameters(&self) -> &[f64] {
        &self.v_params
    }

    /// `Patch(i, j)` with 1-based indices (`AdvApp2Var_ApproxAFunc2Var.cxx:953`).
    fn patch_at(&self, i: i32, j: i32) -> Result<&NetPatch, String> {
        let nb_u = self.nb_patch_in_u();
        let at = (j - 1)
            .checked_mul(nb_u)
            .and_then(|row| row.checked_add(i - 1))
            .ok_or_else(patch_oob)?;
        let at = usize::try_from(at).map_err(|_| patch_oob())?;
        self.patches.get(at).ok_or_else(patch_oob)
    }

    pub(super) fn coeff_counts(&self, i: i32, j: i32) -> Result<(i32, i32), String> {
        let patch = self.patch_at(i, j)?;
        Ok((patch.nb_coeff_u, patch.nb_coeff_v))
    }

    pub(super) fn equation(&self, i: i32, j: i32) -> Result<&[f64], String> {
        Ok(&self.patch_at(i, j)?.equation)
    }
}

fn patch_oob() -> String {
    "AdvApp2Var_ApproxAFunc2Var::ConvertBS: AdvApp2Var_ApproxAFunc2Var.cxx:953".into()
}
