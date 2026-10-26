//! `IntPolyh_StartPoint` (`IntPolyh_StartPoint.hxx:23-126`, `.cxx:22-340`).
//!
//! One point of an intersection polyline: 3D position, the UV parameters on
//! both surfaces, the triangle/edge it was found on with the barycentric
//! lambda, the incidence angle and the chain link index.

/// `MyConfusionPrecision` of `IntPolyh_StartPoint.cxx:21`.
const MY_CONFUSION_PRECISION: f64 = 10.0e-12;

/// `IntPolyh_StartPoint` (`hxx:23`).
#[derive(Debug, Clone, Copy)]
pub struct IntPolyhStartPoint {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub u1: f64,
    pub v1: f64,
    pub u2: f64,
    pub v2: f64,
    pub lambda1: f64,
    pub lambda2: f64,
    pub angle: f64,
    pub t1: i32,
    pub e1: i32,
    pub t2: i32,
    pub e2: i32,
    pub chainlist: i32,
}

impl Default for IntPolyhStartPoint {
    /// `IntPolyh_StartPoint()` (`cxx:24-40`).
    fn default() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            u1: 0.0,
            v1: 0.0,
            u2: 0.0,
            v2: 0.0,
            lambda1: -1.0,
            lambda2: -1.0,
            angle: -2.0,
            t1: -1,
            e1: -2,
            t2: -1,
            e2: -2,
            chainlist: -1,
        }
    }
}

impl IntPolyhStartPoint {
    /// `IntPolyh_StartPoint(x, y, z, u1, v1, u2, v2, t1, e1, lambda1, t2,
    /// e2, lambda2, chainlist)` (`cxx:42-75`).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        x: f64,
        y: f64,
        z: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        t1: i32,
        e1: i32,
        lambda1: f64,
        t2: i32,
        e2: i32,
        lambda2: f64,
        chainlist: i32,
    ) -> Self {
        Self {
            x,
            y,
            z,
            u1,
            v1,
            u2,
            v2,
            lambda1,
            lambda2,
            angle: -2.0,
            t1,
            e1,
            t2,
            e2,
            chainlist,
        }
    }

    /// `GetEdgePoints(Triangle, FirstEdgePoint, SecondEdgePoint, LastPoint)`
    /// (`cxx:150-204`): returns the surface id (1 when the point lies on a
    /// triangle edge of the first surface, 2 for the second, 0 otherwise) and
    /// the three point indices of that triangle starting at the edge.
    pub fn get_edge_points(&self, triangle: &super::triangle::IntPolyhTriangle) -> (i32, i32, i32, i32) {
        let (first, second, third) = triangle.points3();
        let mut first_edge_point = 0;
        let mut second_edge_point = 0;
        let mut last_point = 0;
        let surf_id;

        if self.e1 != -1 {
            match self.e1 {
                1 => {
                    first_edge_point = first;
                    second_edge_point = second;
                    last_point = third;
                }
                2 => {
                    first_edge_point = second;
                    second_edge_point = third;
                    last_point = first;
                }
                3 => {
                    first_edge_point = third;
                    second_edge_point = first;
                    last_point = second;
                }
                _ => {}
            }
            surf_id = 1;
        } else if self.e2 != -1 {
            match self.e2 {
                1 => {
                    first_edge_point = first;
                    second_edge_point = second;
                    last_point = third;
                }
                2 => {
                    first_edge_point = second;
                    second_edge_point = third;
                    last_point = first;
                }
                3 => {
                    first_edge_point = third;
                    second_edge_point = first;
                    last_point = second;
                }
                _ => {}
            }
            surf_id = 2;
        } else {
            surf_id = 0;
        }

        (surf_id, first_edge_point, second_edge_point, last_point)
    }

    /// `CheckSameSP(SP)` (`cxx:261-299`): 1 when the two start points are the
    /// same for chaining purposes.
    pub fn check_same_sp(&self, sp: &Self) -> i32 {
        let mut test = 0;
        if (self.e1 >= -1 && self.e1 == sp.e1) || (self.e2 >= -1 && self.e2 == sp.e2) {
            if (self.lambda1 > -MY_CONFUSION_PRECISION
                && (self.lambda1 - sp.lambda1).abs() < MY_CONFUSION_PRECISION)
                || (self.lambda2 > -MY_CONFUSION_PRECISION
                    && (self.lambda2 - sp.lambda2).abs() < MY_CONFUSION_PRECISION)
            {
                test = 1;
            }
        }
        if test == 0 && (self.e1 == -1 || self.e2 == -1) {
            if (sp.u1 - self.u1).abs() < MY_CONFUSION_PRECISION
                && (sp.v1 - self.v1).abs() < MY_CONFUSION_PRECISION
            {
                test = 1;
            }
        } else if self.e1 == -2 && self.e2 == -2 {
            // `cxx:290-294`: the `.cxx` dumps and prints a diagnostic here;
            // the return value is unchanged.
        }
        test
    }

    /// `SetXYZ(XX, YY, ZZ)` (`cxx:206-212`).
    pub fn set_xyz(&mut self, xx: f64, yy: f64, zz: f64) {
        self.x = xx;
        self.y = yy;
        self.z = zz;
    }

    /// `SetUV1(UU1, VV1)` (`cxx:214-218`).
    pub fn set_uv1(&mut self, uu1: f64, vv1: f64) {
        self.u1 = uu1;
        self.v1 = vv1;
    }

    /// `SetUV2(UU2, VV2)` (`cxx:220-224`).
    pub fn set_uv2(&mut self, uu2: f64, vv2: f64) {
        self.u2 = uu2;
        self.v2 = vv2;
    }

    /// `SetCoupleValue(IT1, IT2)` (`cxx:247-251`).
    pub fn set_couple_value(&mut self, it1: i32, it2: i32) {
        self.t1 = it1;
        self.t2 = it2;
    }
}
