//! `GCPnts_UniformDeflection` and its `CPnts_UniformDeflection` iterator.
//!
//! Sources: OCCT 8.0.0 `GCPnts_UniformDeflection.hxx/.cxx` and
//! `CPnts_UniformDeflection.hxx/.cxx/.lxx`.

use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::{epsilon, CONFUSION};
use occt_geom2d::Curve2d;

use crate::curve::Curve;

#[derive(Clone, Copy, PartialEq, Eq)]
enum DeflectionType {
    Linear,
    Circular,
    Curved,
    Composite,
}

trait UniformCurve {
    fn d0(&self, u: f64) -> GpPnt;
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec);
    fn resolution(&self, tolerance: f64) -> f64;
    fn circle_radius(&self) -> Option<f64>;
    fn circle_parameter_span(&self, first: f64, last: f64) -> f64;
    fn deflection_type(&self) -> DeflectionType;
    fn nb_intervals(&self, continuity: u8) -> i32;
    fn intervals(&self, continuity: u8) -> Vec<f64>;
}

impl<T: Curve + ?Sized> UniformCurve for T {
    fn d0(&self, u: f64) -> GpPnt {
        Curve::d0(self, u)
    }

    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        Curve::d2(self, u)
    }

    fn resolution(&self, tolerance: f64) -> f64 {
        Curve::resolution(self, tolerance)
    }

    fn circle_radius(&self) -> Option<f64> {
        Curve::circle_radius(self)
    }

    fn circle_parameter_span(&self, first: f64, last: f64) -> f64 {
        if let Some((basis_first, basis_last)) = Curve::trimmed_basis_range(self) {
            let natural_span = Curve::last_parameter(self) - Curve::first_parameter(self);
            if natural_span.abs() > f64::EPSILON {
                return (basis_last - basis_first) * (last - first) / natural_span;
            }
        }
        last - first
    }

    fn deflection_type(&self) -> DeflectionType {
        if Curve::nb_intervals(self, 4) > 1 {
            return DeflectionType::Composite;
        }
        if Curve::gp_line(self).is_some() {
            return DeflectionType::Linear;
        }
        if Curve::gp_circ(self).is_some() {
            return DeflectionType::Circular;
        }
        if Curve::bspline_poles(self).is_some_and(|poles| poles.len() == 2)
            || Curve::bezier_poles(self).is_some_and(|poles| poles.len() == 2)
        {
            return DeflectionType::Linear;
        }
        DeflectionType::Curved
    }

    fn nb_intervals(&self, continuity: u8) -> i32 {
        Curve::nb_intervals(self, continuity)
    }

    fn intervals(&self, continuity: u8) -> Vec<f64> {
        Curve::parameter_intervals(self, continuity)
    }
}

struct Curve2dRange<'a> {
    curve: &'a dyn Curve2d,
    basis: &'a dyn Curve2d,
}

impl<'a> Curve2dRange<'a> {
    fn new(curve: &'a dyn Curve2d) -> Self {
        Self {
            curve,
            basis: curve.trimmed_basis().unwrap_or(curve),
        }
    }
}

impl UniformCurve for Curve2dRange<'_> {
    fn d0(&self, u: f64) -> GpPnt {
        let point = self.curve.d0(u);
        GpPnt::new(point.x(), point.y(), 0.0)
    }

    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        let (point, d1, d2) = self.curve.d2(u);
        (
            GpPnt::new(point.x(), point.y(), 0.0),
            GpVec::new(d1.x(), d1.y(), 0.0),
            GpVec::new(d2.x(), d2.y(), 0.0),
        )
    }

    fn resolution(&self, tolerance: f64) -> f64 {
        self.basis.resolution(tolerance)
    }

    fn circle_radius(&self) -> Option<f64> {
        self.basis.gp_circ2d().map(|circle| circle.radius())
    }

    fn circle_parameter_span(&self, first: f64, last: f64) -> f64 {
        last - first
    }

    fn deflection_type(&self) -> DeflectionType {
        if self.basis.nb_intervals(4) > 1 {
            return DeflectionType::Composite;
        }
        if self.basis.is_line() {
            return DeflectionType::Linear;
        }
        if self.basis.gp_circ2d().is_some() {
            return DeflectionType::Circular;
        }
        if self
            .basis
            .bspline2d()
            .is_some_and(|curve| curve.nb_poles() == 2)
            || self.basis.bezier_nb_poles() == Some(2)
        {
            return DeflectionType::Linear;
        }
        DeflectionType::Curved
    }

    fn nb_intervals(&self, continuity: u8) -> i32 {
        self.basis.nb_intervals(continuity)
    }

    fn intervals(&self, continuity: u8) -> Vec<f64> {
        self.basis.parameter_intervals(continuity)
    }
}

struct UniformDeflectionIterator<'a, C: UniformCurve + ?Sized> {
    curve: &'a C,
    done: bool,
    finish: bool,
    tolerance: f64,
    control: bool,
    point_index: i32,
    point_count: i32,
    params: [f64; 3],
    points: [GpPnt; 3],
    max_step: f64,
    deflection: f64,
    first: f64,
    last: f64,
    step: f64,
}

impl<'a, C: UniformCurve + ?Sized> UniformDeflectionIterator<'a, C> {
    fn new(
        curve: &'a C,
        deflection: f64,
        first: f64,
        last: f64,
        tolerance: f64,
        control: bool,
    ) -> Self {
        let first_param = first.min(last);
        let last_param = first.max(last);
        let max_step = last_param - first_param;
        let mut result = Self {
            curve,
            done: true,
            finish: false,
            tolerance,
            control,
            point_index: 0,
            point_count: 0,
            params: [0.0; 3],
            points: [GpPnt::zero(); 3],
            max_step,
            deflection,
            first: first_param,
            last: last_param,
            step: max_step / 2.0,
        };
        result.perform();
        result
    }

    fn perform(&mut self) {
        self.point_index = -1;
        self.point_count = -1;
        let parameter_epsilon = epsilon(self.first);

        while self.point_count < 2 && !self.finish {
            self.point_count += 1;
            let index = self.point_count as usize;
            self.params[index] = self.first;
            let (point, d1, d2) = self.curve.d2(self.first);
            self.points[index] = point;
            let norm_d1 = d1.magnitude();
            if norm_d1 < self.tolerance || d2.magnitude() < self.tolerance {
                self.step = self.max_step.min(1.5 * self.step);
            } else {
                let norm_d2 = d2.cross_magnitude(&d1);
                if norm_d2 / norm_d1 < self.deflection {
                    self.step = self.max_step.min(1.5 * self.step);
                } else {
                    self.step = (8.0 * self.deflection * norm_d1 / norm_d2).sqrt();
                    self.step = self.step.max(self.tolerance).min(self.max_step);
                }
            }

            if self.control {
                self.step = self.step.min(self.last - self.first);
                let end = self.curve.d0(self.first + self.step);
                let middle = self.curve.d0(self.first + self.step / 2.0);
                let chord = GpVec::from_pnts(&self.points[index], &end);
                let chord_length = chord.magnitude();
                if chord_length >= self.deflection {
                    let half_chord = GpVec::from_pnts(&self.points[index], &middle);
                    let mut deviation = half_chord.cross_magnitude(&chord) / chord_length;
                    if deviation > self.deflection / 5.0 {
                        deviation = deviation.max(1.1 * self.deflection);
                        self.step *= (self.deflection / deviation).sqrt();
                        self.step = self.step.max(self.tolerance).min(self.max_step);
                    }
                }
            }

            self.first += self.step;
            self.finish = self.last - self.first < self.tolerance
                || self.step.abs() < self.tolerance
                || self.step < parameter_epsilon;
        }

        if self.finish {
            self.finish_batch();
        }
    }

    fn finish_batch(&mut self) {
        if self.control && self.point_count == 1 {
            let previous = self.params[0];
            if self.last - previous < 0.33 * (self.last - self.first) {
                self.first = (self.last + previous) / 2.0;
                self.params[0] = self.first;
                self.params[1] = self.last;
                self.points[0] = self.curve.d0(self.params[0]);
                self.points[1] = self.curve.d0(self.params[1]);
                return;
            }

            let (last_point, _, _) = self.curve.d2(self.last);
            let previous_point = self.points[0];
            let last_to_previous = GpVec::from_pnts(&last_point, &previous_point);
            let distance = last_to_previous.magnitude();
            if distance < self.deflection {
                self.params[1] = self.last;
                self.points[1] = last_point;
                return;
            }

            self.first = (self.last * (self.params[1] - previous) + previous * self.step)
                / (self.first - previous);
            let adjusted = self.curve.d0(self.first);
            let adjusted_to_previous = GpVec::from_pnts(&adjusted, &previous_point);
            if last_to_previous.cross_magnitude(&adjusted_to_previous) / distance < self.deflection
                && previous >= self.last - self.max_step
            {
                self.params[1] = self.last;
                self.points[1] = last_point;
            } else {
                self.params[1] = self.first;
                self.points[1] = adjusted;
                self.params[2] = self.last;
                self.points[2] = last_point;
                self.point_count += 1;
            }
        } else {
            self.point_count += 1;
            self.point_count = self.point_count.min(2);
            let index = self.point_count as usize;
            self.params[index] = self.last;
            self.points[index] = self.curve.d0(self.last);
        }
    }

    fn more(&mut self) -> bool {
        if !self.done {
            false
        } else if self.point_index == self.point_count {
            if self.finish {
                false
            } else {
                self.perform();
                self.done
            }
        } else {
            self.point_index < self.point_count
        }
    }

    fn next(&mut self) {
        self.point_index += 1;
    }

    fn value(&self) -> f64 {
        self.params[(self.point_index + 1) as usize]
    }

    fn point(&self) -> GpPnt {
        self.points[(self.point_index + 1) as usize]
    }
}

fn hunt(intervals: &[f64], value: f64) -> i32 {
    if intervals.is_empty() {
        return 0;
    }
    if intervals[0] > value {
        return 0;
    }
    if intervals[intervals.len() - 1] < value {
        return intervals.len() as i32 + 1;
    }
    let mut low = 1i32;
    let mut high = intervals.len() as i32;
    while high - low != 1 {
        let middle = (high + low) / 2;
        if intervals[(middle - 1) as usize] < value {
            low = middle;
        } else {
            high = middle;
        }
    }
    low
}

fn perform_curve<C: UniformCurve + ?Sized>(
    curve: &C,
    params: &mut Vec<f64>,
    points: &mut Vec<GpPnt>,
    deflection: f64,
    first: f64,
    last: f64,
    tolerance: f64,
    control: bool,
) -> bool {
    let mut iterator =
        UniformDeflectionIterator::new(curve, deflection, first, last, tolerance, control);
    while iterator.more() {
        params.push(iterator.value());
        points.push(iterator.point());
        iterator.next();
    }
    iterator.done
}

fn perform_composite<C: UniformCurve + ?Sized>(
    curve: &C,
    params: &mut Vec<f64>,
    points: &mut Vec<GpPnt>,
    deflection: f64,
    first: f64,
    last: f64,
    tolerance: f64,
    control: bool,
) -> bool {
    let count = curve.nb_intervals(4);
    let intervals = curve.intervals(4);
    let mut index = hunt(&intervals, first);
    let mut lower = first;
    loop {
        let upper = if index + 1 <= intervals.len() as i32 {
            last.min(intervals[index as usize])
        } else {
            last
        };
        if !perform_curve(
            curve, params, points, deflection, lower, upper, tolerance, control,
        ) {
            return false;
        }
        index += 1;
        if index > count || last < intervals[(index - 1) as usize] {
            return true;
        }
        params.pop();
        points.pop();
        lower = upper;
    }
}

fn control_last<C: UniformCurve + ?Sized>(
    curve: &C,
    params: &mut [f64],
    points: &mut [GpPnt],
    last: f64,
) {
    let count = points.len();
    if count > 2 {
        let previous = params[count - 3];
        let current = params[count - 2];
        if last - current < 0.33 * (last - previous) {
            let adjusted = (last + previous) * 0.5;
            params[count - 2] = adjusted;
            points[count - 2] = curve.d0(adjusted);
        }
    }
}

/// OCCT `GCPnts_UniformDeflection` for 3D and 2D curve adaptors.
pub struct UniformDeflection {
    done: bool,
    deflection: f64,
    params: Vec<f64>,
    points: Vec<GpPnt>,
}

impl Default for UniformDeflection {
    fn default() -> Self {
        Self::new()
    }
}

impl UniformDeflection {
    pub fn new() -> Self {
        Self {
            done: false,
            deflection: 0.0,
            params: Vec::new(),
            points: Vec::new(),
        }
    }

    pub fn initialize(
        &mut self,
        curve: &dyn Curve,
        deflection: f64,
        first: f64,
        last: f64,
        with_control: bool,
    ) {
        self.initialize_curve(curve, deflection, first, last, with_control);
    }

    pub fn initialize_whole(&mut self, curve: &dyn Curve, deflection: f64, with_control: bool) {
        self.initialize(
            curve,
            deflection,
            curve.first_parameter(),
            curve.last_parameter(),
            with_control,
        );
    }

    pub fn initialize_2d(
        &mut self,
        curve: &dyn Curve2d,
        deflection: f64,
        first: f64,
        last: f64,
        with_control: bool,
    ) {
        self.initialize_curve(
            &Curve2dRange::new(curve),
            deflection,
            first,
            last,
            with_control,
        );
    }

    pub fn initialize_2d_whole(
        &mut self,
        curve: &dyn Curve2d,
        deflection: f64,
        with_control: bool,
    ) {
        self.initialize_2d(
            curve,
            deflection,
            curve.first_parameter(),
            curve.last_parameter(),
            with_control,
        );
    }

    fn initialize_curve<C: UniformCurve + ?Sized>(
        &mut self,
        curve: &C,
        deflection: f64,
        first: f64,
        last: f64,
        with_control: bool,
    ) {
        let tolerance = curve.resolution(CONFUSION);
        self.deflection = deflection;
        self.done = false;
        self.params.clear();
        self.points.clear();

        let lower = first.min(last);
        let upper = first.max(last);
        self.done = match curve.deflection_type() {
            DeflectionType::Linear => {
                self.params.extend([lower, upper]);
                self.points.extend([curve.d0(lower), curve.d0(upper)]);
                true
            }
            DeflectionType::Circular => {
                let radius = curve.circle_radius().unwrap();
                let mut angle = (1.0 - deflection / radius).max(0.0);
                angle = 2.0 * angle.acos();
                let angular_span = curve.circle_parameter_span(lower, upper).abs();
                let point_count = (angular_span / angle) as i32 + 2;
                let step = (upper - lower) / (point_count - 1) as f64;
                let mut parameter = lower;
                for _ in 1..=point_count {
                    self.params.push(parameter);
                    self.points.push(curve.d0(parameter));
                    parameter += step;
                }
                true
            }
            DeflectionType::Curved => perform_curve(
                curve,
                &mut self.params,
                &mut self.points,
                deflection,
                lower,
                upper,
                tolerance,
                with_control,
            ),
            DeflectionType::Composite => perform_composite(
                curve,
                &mut self.params,
                &mut self.points,
                deflection,
                lower,
                upper,
                tolerance,
                with_control,
            ),
        };
        control_last(curve, &mut self.params, &mut self.points, upper);
    }

    pub fn is_done(&self) -> bool {
        self.done
    }

    pub fn nb_points(&self) -> Option<i32> {
        self.done.then_some(self.params.len() as i32)
    }

    pub fn parameter(&self, index: i32) -> Option<f64> {
        if !self.done || index < 1 {
            return None;
        }
        self.params.get((index - 1) as usize).copied()
    }

    pub fn value(&self, index: i32) -> Option<GpPnt> {
        if !self.done || index < 1 {
            return None;
        }
        self.points.get((index - 1) as usize).copied()
    }

    pub fn deflection(&self) -> Option<f64> {
        self.done.then_some(self.deflection)
    }
}
