//! Integer vector with 1-based indexing.
//! Source: `math_IntegerVector.hxx` (`math_VectorBase<int>`)

/// Dynamic i32 vector, 1-based indexing. Source: `math_IntegerVector.hxx:48`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MathIntVector {
    data: Vec<i32>,
    lower: usize,
}

impl MathIntVector {
    pub fn new(lower: usize, upper: usize) -> Self {
        assert!(upper >= lower);
        let len = upper - lower + 1;
        Self { data: vec![0; len], lower }
    }

    pub fn with_init(lower: usize, upper: usize, init: i32) -> Self {
        let mut v = Self::new(lower, upper);
        v.init(init);
        v
    }

    pub fn from_slice(values: &[i32]) -> Self { Self { data: values.to_vec(), lower: 1 } }

    #[inline] pub fn len(&self) -> usize { self.data.len() }
    #[inline] pub fn lower(&self) -> usize { self.lower }
    #[inline] pub fn upper(&self) -> usize { self.lower + self.data.len() - 1 }

    pub fn init(&mut self, val: i32) { self.data.fill(val); }

    #[inline] pub fn value(&self, i: usize) -> i32 {
        assert!(i >= self.lower && i <= self.upper()); self.data[i - self.lower]
    }
    #[inline] pub fn set_value(&mut self, i: usize, val: i32) {
        assert!(i >= self.lower && i <= self.upper()); self.data[i - self.lower] = val;
    }

    pub fn max_index(&self) -> usize {
        let (i, _) = self.data.iter().enumerate().fold((0, i32::MIN), |(mi, mv), (i, &v)| if v > mv { (i, v) } else { (mi, mv) });
        self.lower + i
    }

    pub fn min_index(&self) -> usize {
        let (i, _) = self.data.iter().enumerate().fold((0, i32::MAX), |(mi, mv), (i, &v)| if v < mv { (i, v) } else { (mi, mv) });
        self.lower + i
    }
}
