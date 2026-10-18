//! `math_PSO` (`math_PSO.cxx:22-267`, `math_PSO.hxx:60-96`) and
//! `math_PSOParticlesPool` (`math_PSOParticlesPool.cxx:22-66`,
//! `math_PSOParticlesPool.hxx:25-70`), including `PSO_Particle`
//! (`math_PSOParticlesPool.hxx:26-44`) and its `operator<`.
//!
//! Rust note: OCCT's `math_PSO` stores the objective as a
//! `math_MultipleVarFunction*`. Here the objective is passed to `perform*` as
//! `&dyn MathMultipleVarFunction`. The particle arithmetic, the
//! `math_BullardGenerator` draws and the iteration counts are unchanged.

use crate::{MathMultipleVarFunction, MathVector};
use occt_core::math_bullard::BullardGenerator;

/// `RealLast()`.
const REAL_LAST: f64 = f64::MAX;

/// `aBorderDivisor` (`math_PSO.cxx:22`).
const A_BORDER_DIVISOR: f64 = 1.0e+4;

/// `PSO_Particle` (`math_PSOParticlesPool.hxx:26-44`).
#[derive(Debug, Clone)]
pub struct PsoParticle {
    pub position: Vec<f64>,
    pub velocity: Vec<f64>,
    pub best_position: Vec<f64>,
    pub distance: f64,
    pub best_distance: f64,
}

impl PsoParticle {
    /// `PSO_Particle()` (`math_PSOParticlesPool.hxx:35-42`): `Distance` and
    /// `BestDistance` start at `RealLast()`; the three buffers start at 0
    /// (`math_PSOParticlesPool.cxx:32`, `myMemory.Init(0.)`).
    fn new(dimension_count: usize) -> Self {
        Self {
            position: vec![0.0; dimension_count],
            velocity: vec![0.0; dimension_count],
            best_position: vec![0.0; dimension_count],
            distance: REAL_LAST,
            best_distance: REAL_LAST,
        }
    }
}

/// `math_PSOParticlesPool` (`math_PSOParticlesPool.cxx:22-66`).
pub struct PsoParticlesPool {
    particles: Vec<PsoParticle>,
}

impl PsoParticlesPool {
    /// `math_PSOParticlesPool(theParticlesCount, theDimensionCount)`
    /// (`math_PSOParticlesPool.cxx:22-45`). Indexes are shifted to 1-based
    /// (`GetParticle` indexes `1..=myParticlesCount`).
    pub fn new(particles_count: usize, dimension_count: usize) -> Self {
        Self {
            particles: (0..particles_count)
                .map(|_| PsoParticle::new(dimension_count))
                .collect(),
        }
    }

    /// Number of particles held (`myParticlesCount`).
    pub fn len(&self) -> usize {
        self.particles.len()
    }

    /// Whether the pool is empty.
    pub fn is_empty(&self) -> bool {
        self.particles.is_empty()
    }

    /// `GetParticle(theIdx)` (`math_PSOParticlesPool.cxx:52-55`).
    pub fn particle(&self, idx: usize) -> &PsoParticle {
        &self.particles[idx - 1]
    }

    /// `GetParticle(theIdx)` (`math_PSOParticlesPool.cxx:52-55`), mutable.
    pub fn particle_mut(&mut self, idx: usize) -> &mut PsoParticle {
        &mut self.particles[idx - 1]
    }

    /// `GetBestParticle()` (`math_PSOParticlesPool.cxx:58-61`) —
    /// `std::min_element` under `PSO_Particle::operator<`
    /// (`math_PSOParticlesPool.hxx:38`). Ties keep the first particle, as
    /// `std::min_element` does. Returns the 1-based index.
    pub fn best_particle_index(&self) -> usize {
        let mut best = 0;
        for i in 1..self.particles.len() {
            if self.particles[i].distance < self.particles[best].distance {
                best = i;
            }
        }
        best + 1
    }

    /// `GetWorstParticle()` (`math_PSOParticlesPool.cxx:63-66`) —
    /// `std::max_element`. Returns the 1-based index.
    pub fn worst_particle_index(&self) -> usize {
        let mut worst = 0;
        for i in 1..self.particles.len() {
            if self.particles[worst].distance < self.particles[i].distance {
                worst = i;
            }
        }
        worst + 1
    }
}

/// `math_PSO` (`math_PSO.hxx:60-96`).
pub struct Pso {
    low_border: MathVector,
    upp_border: MathVector,
    steps: MathVector,
    n: usize,
    nb_particles: usize,
    nb_iter: usize,
}

impl Pso {
    /// `math_PSO(theFunc, theLowBorder, theUppBorder, theSteps, theNbParticles = 32,
    /// theNbIter = 100)` (`math_PSO.cxx:23-43`, `math_PSO.hxx:70-75`).
    ///
    /// `myN` is `theFunc->NbVariables()` in OCCT; here it is taken from the
    /// border vectors (`low_border.len()`), which the caller sizes the same way.
    pub fn new(
        low_border: MathVector,
        upp_border: MathVector,
        steps: MathVector,
        nb_particles: usize,
        nb_iter: usize,
    ) -> Self {
        let n = low_border.len();
        assert!(upp_border.len() == n && steps.len() == n, "math_PSO: dimension mismatch");
        Self {
            low_border,
            upp_border,
            steps,
            n,
            nb_particles,
            nb_iter,
        }
    }

    /// `aMinUV`/`aMaxUV` (`math_PSO.cxx:63-65`, `:131-132`): the borders shrunk
    /// by `(upp - low) / aBorderDivisor` on each side.
    fn shrunk_borders(&self) -> (MathVector, MathVector) {
        let mut min_uv = MathVector::new(1, self.n);
        let mut max_uv = MathVector::new(1, self.n);
        for i in 1..=self.n {
            let delta = self.upp_border.value(i) - self.low_border.value(i);
            min_uv.set_value(i, self.low_border.value(i) + delta / A_BORDER_DIVISOR);
            max_uv.set_value(i, self.upp_border.value(i) - delta / A_BORDER_DIVISOR);
        }
        (min_uv, max_uv)
    }

    /// `Perform(theParticles, theNbParticles, theValue, theOutPnt)`
    /// (`math_PSO.cxx:44-52`, `math_PSO.hxx:83-88`). The iteration count is
    /// `myNbIter` as set by the constructor, i.e. the `theNbIter = 100` default.
    pub fn perform_with_pool(
        &mut self,
        func: &dyn MathMultipleVarFunction,
        particles: &mut PsoParticlesPool,
        nb_particles: usize,
        value: &mut f64,
        out_pnt: &mut MathVector,
    ) {
        let nb_iter = self.nb_iter;
        self.perform_pso_with_given_particles(
            func,
            particles,
            nb_particles,
            value,
            out_pnt,
            nb_iter,
        );
    }

    /// `Perform(theSteps, theValue, theOutPnt, theNbIter)`
    /// (`math_PSO.cxx:57-121`, `math_PSO.hxx:78-81`): fills the particle pool
    /// from a regular grid over the shrunk borders, then runs the PSO loop.
    pub fn perform(
        &mut self,
        func: &dyn MathMultipleVarFunction,
        steps: &MathVector,
        value: &mut f64,
        out_pnt: &mut MathVector,
        nb_iter: usize,
    ) {
        let (min_uv, max_uv) = self.shrunk_borders();
        self.nb_iter = nb_iter;
        self.steps = steps.clone();

        let mut pool = PsoParticlesPool::new(self.nb_particles, self.n);

        let mut curr_point = MathVector::new(1, self.n);
        for i in 1..=self.n {
            curr_point.set_value(i, min_uv.value(i));
        }

        let mut particle_idx = pool.worst_particle_index();
        loop {
            let mut curr_value = 0.0;
            func.value(&curr_point, &mut curr_value);

            if curr_value < pool.particle(particle_idx).distance {
                let p = pool.particle_mut(particle_idx);
                for dim_idx in 0..self.n {
                    p.position[dim_idx] = curr_point.value(dim_idx + 1);
                    p.best_position[dim_idx] = curr_point.value(dim_idx + 1);
                }
                p.distance = curr_value;
                p.best_distance = curr_value;

                particle_idx = pool.worst_particle_index();
            }

            // Step (`math_PSO.cxx:98-118`).
            curr_point.set_value(1, curr_point.value(1) + self.steps.value(1).max(1.0e-15));
            for dim_idx in 1..self.n {
                if curr_point.value(dim_idx) > max_uv.value(dim_idx) {
                    curr_point.set_value(dim_idx, min_uv.value(dim_idx));
                    curr_point.set_value(
                        dim_idx + 1,
                        curr_point.value(dim_idx + 1) + self.steps.value(dim_idx + 1),
                    );
                } else {
                    break;
                }
            }

            if curr_point.value(self.n) > max_uv.value(self.n) {
                break;
            }
        }

        let nb_particles = self.nb_particles;
        self.perform_pso_with_given_particles(
            func,
            &mut pool,
            nb_particles,
            value,
            out_pnt,
            nb_iter,
        );
    }

    /// `performPSOWithGivenParticles` (`math_PSO.cxx:126-266`).
    fn perform_pso_with_given_particles(
        &mut self,
        func: &dyn MathMultipleVarFunction,
        particles: &mut PsoParticlesPool,
        nb_particles: usize,
        value: &mut f64,
        out_pnt: &mut MathVector,
        nb_iter: usize,
    ) {
        let (min_uv, max_uv) = self.shrunk_borders();
        self.nb_iter = nb_iter;
        self.nb_particles = nb_particles;

        // Current particle data.
        let mut curr_point = MathVector::new(1, self.n);
        let mut best_global_position = MathVector::new(1, self.n);

        // Generate initial particle velocities (`math_PSO.cxx:143-152`).
        let mut random = BullardGenerator::new();
        for part_idx in 1..=self.nb_particles {
            let p = particles.particle_mut(part_idx);
            for dim_idx in 1..=self.n {
                let ksi = random.next_real();
                p.velocity[dim_idx - 1] = self.steps.value(dim_idx) * (ksi - 0.5) * 2.0;
            }
        }

        let best_idx = particles.best_particle_index();
        for dim_idx in 0..self.n {
            best_global_position
                .set_value(dim_idx + 1, particles.particle(best_idx).position[dim_idx]);
        }
        let mut best_global_distance = particles.particle(best_idx).distance;

        // This velocity is used for detecting stagnation state.
        let mut termination_velocity = MathVector::new(1, self.n);
        for i in 1..=self.n {
            termination_velocity.set_value(i, self.steps.value(i) / 2048.0);
        }
        // This velocity shows minimal velocity on current step.
        let mut minimal_velocity = MathVector::new(1, self.n);

        // Run PSO iterations (`math_PSO.cxx:169-265`): `aStep` runs to `myNbIter - 1`.
        for step in 1..self.nb_iter {
            minimal_velocity.init(REAL_LAST);

            for part_idx in 1..=self.nb_particles {
                let ksi1 = random.next_real();
                let ksi2 = random.next_real();

                const A_RETENT_WEIGHT: f64 = 0.72900;
                const A_PERSON_WEIGHT: f64 = 1.49445;
                const A_SOCIAL_WEIGHT: f64 = 1.49445;

                let p = particles.particle_mut(part_idx);
                for dim_idx in 0..self.n {
                    p.velocity[dim_idx] = p.velocity[dim_idx] * A_RETENT_WEIGHT
                        + (p.best_position[dim_idx] - p.position[dim_idx])
                            * (A_PERSON_WEIGHT * ksi1)
                        + (best_global_position.value(dim_idx + 1) - p.position[dim_idx])
                            * (A_SOCIAL_WEIGHT * ksi2);

                    p.position[dim_idx] += p.velocity[dim_idx];
                    p.position[dim_idx] = p.position[dim_idx]
                        .max(min_uv.value(dim_idx + 1))
                        .min(max_uv.value(dim_idx + 1));
                    curr_point.set_value(dim_idx + 1, p.position[dim_idx]);

                    minimal_velocity.set_value(
                        dim_idx + 1,
                        minimal_velocity
                            .value(dim_idx + 1)
                            .min(p.velocity[dim_idx].abs()),
                    );
                }

                let mut a_distance = 0.0;
                func.value(&curr_point, &mut a_distance);
                particles.particle_mut(part_idx).distance = a_distance;

                let p = particles.particle_mut(part_idx);
                if p.distance < p.best_distance {
                    p.best_distance = p.distance;
                    for dim_idx in 0..self.n {
                        p.best_position[dim_idx] = p.position[dim_idx];
                    }

                    if p.distance < best_global_distance {
                        best_global_distance = p.distance;
                        for dim_idx in 0..self.n {
                            best_global_position.set_value(dim_idx + 1, p.position[dim_idx]);
                        }
                    }
                }
            }

            let mut is_terminal_velocity_reached = true;
            for dim_idx in 1..=self.n {
                if minimal_velocity.value(dim_idx) > termination_velocity.value(dim_idx) {
                    is_terminal_velocity_reached = false;
                    break;
                }
            }

            if is_terminal_velocity_reached {
                // Minimum number of steps.
                const A_MIN_STEPS: usize = 16;

                if step > A_MIN_STEPS {
                    break;
                }

                for part_idx in 1..=self.nb_particles {
                    let ksi = random.next_real();
                    let p = particles.particle_mut(part_idx);

                    for dim_idx in 0..self.n {
                        let at_min = p.position[dim_idx] == min_uv.value(dim_idx + 1);
                        let at_max = p.position[dim_idx] == max_uv.value(dim_idx + 1);
                        if at_min || at_max {
                            p.velocity[dim_idx] = if at_min {
                                self.steps.value(dim_idx + 1) * ksi
                            } else {
                                -self.steps.value(dim_idx + 1) * ksi
                            };
                        } else {
                            p.velocity[dim_idx] = self.steps.value(dim_idx + 1) * (ksi - 0.5) * 2.0;
                        }
                    }
                }
            }
        }

        *value = best_global_distance;
        for i in 1..=self.n {
            out_pnt.set_value(i, best_global_position.value(i));
        }
    }

    /// `NbParticles()`-equivalent state used by callers.
    pub fn nb_particles(&self) -> usize {
        self.nb_particles
    }

    /// `myN`.
    pub fn nb_variables(&self) -> usize {
        self.n
    }
}
