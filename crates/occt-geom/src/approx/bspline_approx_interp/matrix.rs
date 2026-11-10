use occt_math::MathMatrix;

/// `math_Matrix::Multiplied(math_Matrix)`.
pub(super) fn matrix_multiply(a: &MathMatrix, b: &MathMatrix) -> MathMatrix {
    let rows = a.row_count();
    let cols = b.col_count();
    let inner = a.col_count();
    let mut out = MathMatrix::with_init(1, rows, 1, cols, 0.0);
    for i in 1..=rows {
        for j in 1..=cols {
            let mut sum = 0.0;
            for k in 1..=inner {
                sum += a.value(i, k) * b.value(b.lower_row() + k - 1, j);
            }
            out.set_value(i, j, sum);
        }
    }
    out
}
