use approx::{AbsDiffEq, RelativeEq, UlpsEq};
use crate::tuple::Tuple;
use auto_ops::*;

#[derive(Clone, Debug, PartialEq)]
pub struct Matrix {
    data: Vec<f64>,
    cols: usize,
    rows: usize,
}

impl Matrix {
    pub fn new(cols: usize, rows: usize, init: f64) -> Matrix {
        let data = vec![init; cols * rows];
        Matrix { data, cols, rows }
    }

    pub fn from_vec(cols: usize, rows: usize, data: Vec<f64>) -> Matrix {
        Matrix { data, cols, rows }
    }

    pub fn get(&self, x: usize, y: usize) -> f64 {
        assert!(x < self.cols && y < self.rows);
        self.data[y + self.rows * x]
    }

    pub fn set(&mut self, x: usize, y: usize, value: f64) {
        assert!(x < self.cols && y < self.rows);
        self.data[y + self.rows * x] = value;
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn width(&self) -> usize {
        self.cols
    }

    pub fn height(&self) -> usize {
        self.rows
    }

    pub fn get_col(&self, col: usize) -> Tuple {
        assert!(col < self.cols);
        let mut result = Tuple::new(self.rows, 0.);
        for row in 0..self.rows {
            result.set(row, self.get(col, row));
        }
        result
    }

    pub fn get_row(&self, row: usize) -> Tuple {
        assert!(row < self.rows);
        let mut result = Tuple::new(self.cols, 0.);
        for col in 0..self.cols {
            result.set(col, self.get(col, row));
        }
        result
    }

    pub fn identity(dimension: usize) -> Matrix {
        let mut result = Matrix::new(dimension, dimension, 0.);
        for i in 0..dimension {
            result.set(i, i, 1.);
        }
        result
    }

    pub fn translation(x: f64, y: f64, z: f64) -> Matrix {
        let mut result = Matrix::identity(4);
        result.set(3, 0, x);
        result.set(3, 1, y);
        result.set(3, 2, z);
        result
    }

    pub fn scaling(x: f64, y: f64, z: f64) -> Matrix {
        let mut result = Matrix::identity(4);
        result.set(0, 0, x);
        result.set(1, 1, y);
        result.set(2, 2, z);
        result
    }

    pub fn shearing(xy: f64, xz: f64, yx: f64, yz: f64, zx: f64, zy: f64) -> Matrix {
        let mut result = Matrix::identity(4);
        result.set(1, 0, xy);
        result.set(2, 0, xz);
        result.set(0, 1, yx);
        result.set(2, 1, yz);
        result.set(0, 2, zx);
        result.set(1, 2, zy);
        result
    }

    pub fn transpose(&self) -> Matrix {
        let mut result = Matrix::new(self.rows, self.cols, 0.);
        for col in 0..self.cols {
            for row in 0..self.rows {
                result.set(row, col, self.get(col, row));
            }
        }
        result
    }

    pub fn submatrix(&self, x: usize, y: usize) -> Matrix {
        assert_eq!(self.cols, self.rows);
        let mut result = Matrix::new(self.cols - 1, self.rows - 1, 0.);
        for col in 0..self.cols {
            for row in 0..self.rows {
                if row != y && col != x {
                    let new_x = if col > x { col - 1 } else { col };
                    let new_y = if row > y { row - 1 } else { row };
                    result.set(new_x, new_y, self.get(col, row));
                }
            }
        }
        result
    }

    pub fn translate(&self, x: f64, y: f64, z: f64) -> Matrix {
        Matrix::translation(x, y, z) * self
    }

    pub fn scale(&self, x: f64, y: f64, z: f64) -> Matrix {
        Matrix::scaling(x, y, z) * self
    }

    pub fn shear(&self, xy: f64, xz: f64, yx: f64, yz: f64, zx: f64, zy: f64) -> Matrix {
        Matrix::shearing(xy, xz, yx, yz, zx, zy) * self
    }

    // function determinant(M)
    //   det ← 0
    //
    //   if M.size = 2
    //     det ← M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
    //
    //   else
    //     for column ← 0 to M.size - 1
    //       det ← det + M[0, column] * cofactor(M, 0, column)
    //     end for
    //   end if
    //
    //   return det
    // end function
    pub fn determinant(&self) -> f64 {
        assert_eq!(self.cols, self.rows);
        if self.cols == 2 {
            self.get(0, 0) * self.get(1, 1) - self.get(1, 0) * self.get(0, 1)
        } else {
            let mut det = 0.;
            for col in 0..self.cols {
                det += self.get(col, 0) * self.cofactor(col, 0);
            }
            det
        }
    }

    pub fn is_invertible(&self) -> bool {
        self.determinant() != 0.
    }

    pub fn minor(&self, x: usize, y: usize) -> f64 {
        self.submatrix(x, y).determinant()
    }

    pub fn cofactor(&self, x: usize, y: usize) -> f64 {
        let minor = self.minor(x, y);
        if (y + x) % 2 == 0 {
            minor
        } else {
            -minor
        }
    }

    // function inverse(M)
    //   fail if M is not invertible
    //
    //   M2 ← new matrix of same size as M
    //
    //   for row ← 0 to M.size - 1
    //     for col ← 0 to M.size - 1
    //       c ← cofactor(M, row, col)
    //
    //       # note that "col, row" here, instead of "row, col",
    //       # accomplishes the transpose operation!
    //       M2[col, row] ← c / determinant(M)
    //     end for
    //   end for
    //
    //   return M2
    // end function
    pub fn inverse(&self) -> Matrix {
        let mut result = Matrix::new(self.cols, self.rows, 0.);
        let det = self.determinant();
        for col in 0..self.cols {
            for row in 0..self.rows {
                result.set(row, col, self.cofactor(col, row) / det);
            }
        }
        result
    }

    pub fn rotation_x(radians: f64) -> Matrix {
        let mut result = Matrix::identity(4);
        result.set(1, 1, f64::cos(radians));
        result.set(1, 2, f64::sin(radians));
        result.set(2, 1, -f64::sin(radians));
        result.set(2, 2, f64::cos(radians));
        result
    }

    pub fn rotation_y(radians: f64) -> Matrix {
        let mut result = Matrix::identity(4);
        result.set(0, 0, f64::cos(radians));
        result.set(0, 2, -f64::sin(radians));
        result.set(2, 0, f64::sin(radians));
        result.set(2, 2, f64::cos(radians));
        result
    }

    pub fn rotation_z(radians: f64) -> Matrix {
        let mut result = Matrix::identity(4);
        result.set(0, 0, f64::cos(radians));
        result.set(0, 1, f64::sin(radians));
        result.set(1, 0, -f64::sin(radians));
        result.set(1, 1, f64::cos(radians));
        result
    }

    pub fn rotate_x(&self, radians: f64) -> Matrix {
        Matrix::rotation_x(radians) * self
    }

    pub fn rotate_y(&self, radians: f64) -> Matrix {
        Matrix::rotation_y(radians) * self
    }

    pub fn rotate_z(&self, radians: f64) -> Matrix {
        Matrix::rotation_z(radians) * self
    }
}

// function matrix_multiply(A, B)
//   M ← matrix()
//
//   for row ← 0 to 3
//     for col ← 0 to 3
//       M[row, col] ← A[row, 0] * B[0, col] +
//                      A[row, 1] * B[1, col] +
//                      A[row, 2] * B[2, col] +
//                      A[row, 3] * B[3, col]
//     end for
//   end for
//
//   return M
// end function
impl_op_ex!(*|a: &Matrix, b: &Matrix| -> Matrix {
    assert_eq!(a.cols, b.rows);
    let mut result = Matrix::new(a.cols, b.rows, 0.);
    for col in 0..b.cols {
        for row in 0..a.rows {
            let tuple_row = a.get_row(row);
            let tuple_col = b.get_col(col);
            result.set(col, row, tuple_col.dot(&tuple_row));
        }
    }
    result
});

impl_op_ex!(*|a: &Matrix, b: &Tuple| -> Tuple {
    assert_eq!(a.rows, b.dimension());
    let mut result = Tuple::new(a.cols(), 0.);
    for row in 0..a.rows {
        let tuple_row = a.get_row(row);
        result.set(row, b.dot(&tuple_row));
    }
    result
});

// constant EPSILON ← 0.00001
//
// function equal(a, b)
//   if abs(a - b) < EPSILON
//     return true
//   else
//     return false
//   end if
// end function
// Feature: Matrices
impl AbsDiffEq for Matrix {
    type Epsilon = f64;

    fn default_epsilon() -> f64 {
        1e-5_f64
    }

    fn abs_diff_eq(&self, other: &Matrix, epsilon: f64) -> bool {
        if self.cols != other.cols || self.rows != other.rows {
            return false;
        }
        for x in 0..self.cols {
            for y in 0..self.rows {
                if !self.get(y, x).abs_diff_eq(&other.get(y, x), epsilon) {
                    return false;
                }
            }
        }
        true
    }
}

impl RelativeEq for Matrix {
    fn default_max_relative() -> f64 {
        f64::default_max_relative()
    }

    fn relative_eq(&self, other: &Matrix, epsilon: f64, max_relative: f64) -> bool {
        if self.cols != other.cols || self.rows != other.rows {
            return false;
        }
        for x in 0..self.cols {
            for y in 0..self.rows {
                if !self
                    .get(y, x)
                    .relative_eq(&other.get(y, x), epsilon, max_relative)
                {
                    return false;
                }
            }
        }
        true
    }
}

impl UlpsEq for Matrix {
    fn default_max_ulps() -> u32 {
        f64::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Matrix, epsilon: f64, max_ulps: u32) -> bool {
        if self.cols != other.cols || self.rows != other.rows {
            return false;
        }
        for x in 0..self.cols {
            for y in 0..self.rows {
                if !self.get(y, x).ulps_eq(&other.get(y, x), epsilon, max_ulps) {
                    return false;
                }
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use std::f64::consts::PI;
    use std::f64::consts::SQRT_2;

    use super::*;

    // Scenario: Constructing and inspecting a 4x4 matrix
    //   Given the following 4x4 matrix M:
    //     |  1   |  2   |  3   |  4   |
    //     |  5.5 |  6.5 |  7.5 |  8.5 |
    //     |  9   | 10   | 11   | 12   |
    //     | 13.5 | 14.5 | 15.5 | 16.5 |
    //   Then M[0,0] = 1
    //     And M[0,3] = 4
    //     And M[1,0] = 5.5
    //     And M[1,2] = 7.5
    //     And M[2,2] = 11
    //     And M[3,0] = 13.5
    //     And M[3,2] = 15.5
    #[test]
    fn constructing_4x4_matrix() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                1., 5.5, 9., 13.5, 2., 6.5, 10., 14.5, 3., 7.5, 11., 15.5, 4., 8.5, 12., 16.5,
            ],
        );
        assert_eq!(matrix.get(0, 0), 1.);
        assert_eq!(matrix.get(3, 0), 4.);
        assert_eq!(matrix.get(0, 1), 5.5);
        assert_eq!(matrix.get(2, 1), 7.5);
        assert_eq!(matrix.get(2, 2), 11.);
        assert_eq!(matrix.get(0, 3), 13.5);
        assert_eq!(matrix.get(2, 3), 15.5);
    }

    // Scenario: A 2x2 matrix ought to be representable
    //   Given the following 2x2 matrix M:
    //     | -3 |  5 |
    //     |  1 | -2 |
    //   Then M[0,0] = -3
    //     And M[0,1] = 5
    //     And M[1,0] = 1
    //     And M[1,1] = -2
    #[test]
    fn constructing_2x2_matrix() {
        let matrix = Matrix::from_vec(2, 2, vec![-3., 1., 5., -2.]);
        assert_eq!(matrix.get(0, 0), -3.);
        assert_eq!(matrix.get(1, 0), 5.);
        assert_eq!(matrix.get(0, 1), 1.);
        assert_eq!(matrix.get(1, 1), -2.);
    }

    // Scenario: A 3x3 matrix ought to be representable
    //   Given the following 3x3 matrix M:
    //     | -3 |  5 |  0 |
    //     |  1 | -2 | -7 |
    //     |  0 |  1 |  1 |
    //   Then M[0,0] = -3
    //     And M[1,1] = -2
    //     And M[2,2] = 1
    #[test]
    fn construction_3x3_matrix() {
        let matrix = Matrix::from_vec(3, 3, vec![-3., 1., 0., 5., -2., 1., 0., -7., 1.]);
        assert_eq!(matrix.get(0, 0), -3.);
        assert_eq!(matrix.get(1, 1), -2.);
        assert_eq!(matrix.get(2, 2), 1.);
    }

    // Scenario: Matrix equality with identical matrices
    //   Given the following matrix A:
    //       | 1 | 2 | 3 | 4 |
    //       | 5 | 6 | 7 | 8 |
    //       | 9 | 8 | 7 | 6 |
    //       | 5 | 4 | 3 | 2 |
    //     And the following matrix B:
    //       | 1 | 2 | 3 | 4 |
    //       | 5 | 6 | 7 | 8 |
    //       | 9 | 8 | 7 | 6 |
    //       | 5 | 4 | 3 | 2 |
    //   Then A = B
    #[test]
    fn matrix_equality_with_identical_matrices() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                1., 2., 3., 4., 5., 6., 7., 8., 9., 8., 7., 6., 5., 4., 3., 2.,
            ],
        );
        let expected = Matrix::from_vec(
            4,
            4,
            vec![
                1., 2., 3., 4., 5., 6., 7., 8., 9., 8., 7., 6., 5., 4., 3., 2.,
            ],
        );
        assert_eq!(matrix, expected);
    }

    // Scenario: Matrix equality with different matrices
    //   Given the following matrix A:
    //       | 1 | 2 | 3 | 4 |
    //       | 5 | 6 | 7 | 8 |
    //       | 9 | 8 | 7 | 6 |
    //       | 5 | 4 | 3 | 2 |
    //     And the following matrix B:
    //       | 2 | 3 | 4 | 5 |
    //       | 6 | 7 | 8 | 9 |
    //       | 8 | 7 | 6 | 5 |
    //       | 4 | 3 | 2 | 1 |
    //   Then A != B
    #[test]
    fn matrix_equality_with_different_matrices() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                1., 2., 3., 4., 5., 6., 7., 8., 9., 8., 7., 6., 5., 4., 3., 2.,
            ],
        );
        let expected = Matrix::from_vec(
            4,
            4,
            vec![
                2., 3., 4., 5., 6., 7., 8., 9., 8., 7., 6., 5., 4., 3., 2., 1.,
            ],
        );
        assert_ne!(matrix, expected);
    }

    // Scenario: Multiplying two matrices
    //   Given the following matrix A:
    //       | 1 | 2 | 3 | 4 |
    //       | 5 | 6 | 7 | 8 |
    //       | 9 | 8 | 7 | 6 |
    //       | 5 | 4 | 3 | 2 |
    //     And the following matrix B:
    //       | -2 | 1 | 2 |  3 |
    //       |  3 | 2 | 1 | -1 |
    //       |  4 | 3 | 6 |  5 |
    //       |  1 | 2 | 7 |  8 |
    //   Then A * B is the following 4x4 matrix:
    //       | 20|  22 |  50 |  48 |
    //       | 44|  54 | 114 | 108 |
    //       | 40|  58 | 110 | 102 |
    //       | 16|  26 |  46 |  42 |
    #[test]
    fn multiplying_two_matrices() {
        let matrix_a = Matrix::from_vec(
            4,
            4,
            vec![
                1., 5., 9., 5., 2., 6., 8., 4., 3., 7., 7., 3., 4., 8., 6., 2.,
            ],
        );
        let matrix_b = Matrix::from_vec(
            4,
            4,
            vec![
                -2., 3., 4., 1., 1., 2., 3., 2., 2., 1., 6., 7., 3., -1., 5., 8.,
            ],
        );
        let expected = Matrix::from_vec(
            4,
            4,
            vec![
                20., 44., 40., 16., 22., 54., 58., 26., 50., 114., 110., 46., 48., 108., 102., 42.,
            ],
        );
        assert_eq!(matrix_a * matrix_b, expected);
    }

    // Scenario: A matrix multiplied by a tuple
    //   Given the following matrix A:
    //       | 1 | 2 | 3 | 4 |
    //       | 2 | 4 | 4 | 2 |
    //       | 8 | 6 | 4 | 1 |
    //       | 0 | 0 | 0 | 1 |
    //     And b ← tuple(1, 2, 3, 1)
    //   Then A * b = tuple(18, 24, 33, 1)
    #[test]
    fn multiplying_matrix_by_tuple() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                1., 2., 8., 0., 2., 4., 6., 0., 3., 4., 4., 0., 4., 2., 1., 1.,
            ],
        );
        let point = Tuple::point(1., 2., 3.);
        let expected = Tuple::point(18., 24., 33.);
        assert_eq!(matrix * point, expected);
    }

    #[test]
    fn get_row_returns_correct_row() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                1., 5., 9., 3., 2., 6., 0., 4., 3., 7., 1., 5., 4., 8., 2., 6.,
            ],
        );
        let expected = Tuple::from_vec(vec![1., 2., 3., 4.]);
        assert_eq!(matrix.get_row(0), expected);
    }

    #[test]
    fn get_col_returns_correct_col() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                1., 5., 9., 3., 2., 6., 0., 4., 3., 7., 1., 5., 4., 8., 2., 6.,
            ],
        );
        let expected = Tuple::from_vec(vec![1., 5., 9., 3.]);
        assert_eq!(matrix.get_col(0), expected);
    }

    // Scenario: Multiplying a matrix by the identity matrix
    //   Given the following matrix A:
    //     | 0 | 1 |  2 |  4 |
    //     | 1 | 2 |  4 |  8 |
    //     | 2 | 4 |  8 | 16 |
    //     | 4 | 8 | 16 | 32 |
    //   Then A * identity_matrix = A
    #[test]
    fn multiplying_matrix_by_identity_matrix() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                0., 1., 2., 4., 1., 2., 4., 8., 2., 4., 8., 16., 4., 8., 16., 32.,
            ],
        );
        assert_eq!(&matrix * Matrix::identity(4), matrix);
    }

    // Scenario: Multiplying the identity matrix by a tuple
    //   Given a ← tuple(1, 2, 3, 4)
    //   Then identity_matrix * a = a
    #[test]
    fn multiplying_the_identity_matrix_by_a_tuple() {
        let identity = Matrix::identity(4);
        let a = Tuple::from_vec(vec![1., 2., 3., 4.]);
        let result = &identity * &a;
        assert_eq!(result, a);
    }

    // Scenario: Transposing a matrix
    //   Given the following matrix A:
    //     | 0 | 9 | 3 | 0 |
    //     | 9 | 8 | 0 | 8 |
    //     | 1 | 8 | 5 | 3 |
    //     | 0 | 0 | 5 | 8 |
    //   Then transpose(A) is the following matrix:
    //     | 0 | 9 | 1 | 0 |
    //     | 9 | 8 | 8 | 0 |
    //     | 3 | 0 | 5 | 5 |
    //     | 0 | 8 | 3 | 8 |
    #[test]
    fn transposing_matrix() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                0., 9., 3., 0., 9., 8., 0., 8., 1., 8., 5., 3., 0., 0., 5., 8.,
            ],
        );
        let expected = Matrix::from_vec(
            4,
            4,
            vec![
                0., 9., 1., 0., 9., 8., 8., 0., 3., 0., 5., 5., 0., 8., 3., 8.,
            ],
        );
        assert_eq!(matrix.transpose(), expected);
    }

    // Scenario: Transposing the identity matrix
    //   Given A ← transpose(identity_matrix)
    //   Then A = identity_matrix
    #[test]
    fn transposing_identity_matrix() {
        let matrix = Matrix::identity(4);
        assert_eq!(matrix.transpose(), matrix);
    }

    // Scenario: Calculating the determinant of a 2x2 matrix
    //   Given the following 2x2 matrix A:
    //     |  1 | 5 |
    //     | -3 | 2 |
    //   Then determinant(A) = 17
    #[test]
    fn calculating_determinant_of_2x2_matrix() {
        let matrix = Matrix::from_vec(2, 2, vec![1., -3., 5., 2.]);
        assert_eq!(matrix.determinant(), 17.);
    }

    // Scenario: A submatrix of a 3x3 matrix is a 2x2 matrix
    //   Given the following 3x3 matrix A:
    //     |  1 | 5 |  0 |
    //     | -3 | 2 |  7 |
    //     |  0 | 6 | -3 |
    //   Then submatrix(A, 0, 2) is the following 2x2 matrix:
    //     | -3 | 2 |
    //     |  0 | 6 |
    #[test]
    fn submatrix_of_3x3_is_2x2() {
        let matrix = Matrix::from_vec(3, 3, vec![1., -3., 0., 5., 2., 6., 0., 7., -3.]);
        let expected = Matrix::from_vec(2, 2, vec![-3., 0., 2., 6.]);
        assert_eq!(matrix.submatrix(2, 0), expected);
    }

    // Scenario: A submatrix of a 4x4 matrix is a 3x3 matrix
    //   Given the following 4x4 matrix A:
    //     | -6 |  1 |  1 |  6 |
    //     | -8 |  5 |  8 |  6 |
    //     | -1 |  0 |  8 |  2 |
    //     | -7 |  1 | -1 |  1 |
    //   Then submatrix(A, 2, 1) is the following 3x3 matrix:
    //     | -6 |  1 | 6 |
    //     | -8 |  8 | 6 |
    //     | -7 | -1 | 1 |
    #[test]
    fn submatrix_of_4x4_is_3x3() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                -6., -8., -1., -7., 1., 5., 0., 1., 1., 8., 8., -1., 6., 6., 2., 1.,
            ],
        );
        let expected = Matrix::from_vec(3, 3, vec![-6., -8., -7., 1., 8., -1., 6., 6., 1.]);
        assert_eq!(matrix.submatrix(1, 2), expected);
    }

    // Scenario: Calculating a minor of a 3x3 matrix
    //   Given the following 3x3 matrix A:
    //       |  3 |  5 |  0 |
    //       |  2 | -1 | -7 |
    //       |  6 | -1 |  5 |
    //     And B ← submatrix(A, 1, 0)
    //   Then determinant(B) = 25
    //     And minor(A, 1, 0) = 25
    #[test]
    fn calculating_minor_of_3x3_matrix() {
        let matrix = Matrix::from_vec(3, 3, vec![1., 2., 6., 5., -1., -1., 0., -7., 5.]);
        let submatrix = matrix.submatrix(0, 1);
        assert_eq!(submatrix.determinant(), 25.);
        assert_eq!(matrix.minor(0, 1), 25.);
    }

    // Scenario: Calculating a cofactor of a 3x3 matrix
    //   Given the following 3x3 matrix A:
    //       |  3 |  5 |  0 |
    //       |  2 | -1 | -7 |
    //       |  6 | -1 |  5 |
    //   Then minor(A, 0, 0) = -12
    //     And cofactor(A, 0, 0) = -12
    //     And minor(A, 1, 0) = 25
    //     And cofactor(A, 1, 0) = -25
    #[test]
    fn calculating_cofactor_of_3x3_matrix() {
        let matrix = Matrix::from_vec(3, 3, vec![3., 2., 6., 5., -1., -1., 0., -7., 5.]);
        assert_eq!(matrix.minor(0, 0), -12.);
        assert_eq!(matrix.cofactor(0, 0), -12.);
        assert_eq!(matrix.minor(0, 1), 25.);
        assert_eq!(matrix.cofactor(0, 1), -25.);
    }

    // Scenario: Calculating the determinant of a 3x3 matrix
    //   Given the following 3x3 matrix A:
    //     |  1 |  2 |  6 |
    //     | -5 |  8 | -4 |
    //     |  2 |  6 |  4 |
    //   Then cofactor(A, 0, 0) = 56
    //     And cofactor(A, 0, 1) = 12
    //     And cofactor(A, 0, 2) = -46
    //     And determinant(A) = -196
    #[test]
    fn calculating_determinant_of_3x3_matrix() {
        let matrix = Matrix::from_vec(3, 3, vec![1., -5., 2., 2., 8., 6., 6., -4., 4.]);
        assert_eq!(matrix.cofactor(0, 0), 56.);
        assert_eq!(matrix.cofactor(1, 0), 12.);
        assert_eq!(matrix.cofactor(2, 0), -46.);
        assert_eq!(matrix.determinant(), -196.);
    }

    // Scenario: Calculating the determinant of a 4x4 matrix
    //   Given the following 4x4 matrix A:
    //     | -2 | -8 |  3 |  5 |
    //     | -3 |  1 |  7 |  3 |
    //     |  1 |  2 | -9 |  6 |
    //     | -6 |  7 |  7 | -9 |
    //   Then cofactor(A, 0, 0) = 690
    //     And cofactor(A, 0, 1) = 447
    //     And cofactor(A, 0, 2) = 210
    //     And cofactor(A, 0, 3) = 51
    //     And determinant(A) = -4071
    #[test]
    fn calculating_determinant_of_4x4_matrix() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                -2., -3., 1., -6., -8., 1., 2., 7., 3., 7., -9., 7., 5., 3., 6., -9.,
            ],
        );
        assert_eq!(matrix.cofactor(0, 0), 690.);
        assert_eq!(matrix.cofactor(1, 0), 447.);
        assert_eq!(matrix.cofactor(2, 0), 210.);
        assert_eq!(matrix.cofactor(3, 0), 51.);
        assert_eq!(matrix.determinant(), -4071.);
    }

    // Scenario: Testing an invertible matrix for invertibility
    //   Given the following 4x4 matrix A:
    //     |  6 |  4 |  4 |  4 |
    //     |  5 |  5 |  7 |  6 |
    //     |  4 | -9 |  3 | -7 |
    //     |  9 |  1 |  7 | -6 |
    //   Then determinant(A) = -2120
    //     And A is invertible
    #[test]
    fn testing_invertible_matrix_for_invertibility() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                6., 5., 4., 9., 4., 5., -9., 1., 4., 7., 3., 7., 4., 6., -7., -6.,
            ],
        );
        assert_eq!(matrix.determinant(), -2120.);
        assert!(matrix.is_invertible());
    }

    // Scenario: Testing a noninvertible matrix for invertibility
    //   Given the following 4x4 matrix A:
    //     | -4 |  2 | -2 | -3 |
    //     |  9 |  6 |  2 |  6 |
    //     |  0 | -5 |  1 | -5 |
    //     |  0 |  0 |  0 |  0 |
    //   Then determinant(A) = 0
    //     And A is not invertible
    #[test]
    fn testing_non_invertible_matrix_for_invertibility() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                -4., 9., 0., 0., 2., 6., -5., 0., -2., 2., 1., 0., -3., 6., -5., 0.,
            ],
        );
        assert_eq!(matrix.determinant(), 0.);
        assert!(!matrix.is_invertible());
    }

    // Scenario: Calculating the inverse of a matrix
    //   Given the following 4x4 matrix A:
    //       | -5 |  2 |  6 | -8 |
    //       |  1 | -5 |  1 |  8 |
    //       |  7 |  7 | -6 | -7 |
    //       |  1 | -3 |  7 |  4 |
    //     And B ← inverse(A)
    //   Then determinant(A) = 532
    //     And cofactor(A, 2, 3) = -160
    //     And B[3,2] = -160/532
    //     And cofactor(A, 3, 2) = 105
    //     And B[2,3] = 105/532
    //     And B is the following 4x4 matrix:
    //       |  0.21805 |  0.45113 |  0.24060 | -0.04511 |
    //       | -0.80827 | -1.45677 | -0.44361 |  0.52068 |
    //       | -0.07895 | -0.22368 | -0.05263 |  0.19737 |
    //       | -0.52256 | -0.81391 | -0.30075 |  0.30639 |
    #[test]
    fn calculating_inverse_of_matrix() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                -5., 1., 7., 1., 2., -5., 7., -3., 6., 1., -6., 7., -8., 8., -7., 4.,
            ],
        );
        let inv = matrix.inverse();
        let expected = Matrix::from_vec(
            4,
            4,
            vec![
                116. / 532.,
                -430. / 532.,
                -42. / 532.,
                -278. / 532.,
                240. / 532.,
                -775. / 532.,
                -119. / 532.,
                -433. / 532.,
                128. / 532.,
                -236. / 532.,
                -28. / 532.,
                -160. / 532.,
                -24. / 532.,
                277. / 532.,
                105. / 532.,
                163. / 532.,
            ],
        );
        assert_eq!(matrix.determinant(), 532.);
        assert_eq!(matrix.cofactor(3, 2), -160.);
        assert_eq!(inv.get(2, 3), -160. / 532.);
        assert_eq!(matrix.cofactor(2, 3), 105.);
        assert_eq!(inv.get(3, 2), 105. / 532.);
        assert_eq!(inv, expected);
    }

    // Scenario: Calculating the inverse of another matrix
    //   Given the following 4x4 matrix A:
    //     |  8 | -5 |  9 |  2 |
    //     |  7 |  5 |  6 |  1 |
    //     | -6 |  0 |  9 |  6 |
    //     | -3 |  0 | -9 | -4 |
    //   Then inverse(A) is the following 4x4 matrix:
    //     | -0.15385 | -0.15385 | -0.28205 | -0.53846 |
    //     | -0.07692 |  0.12308 |  0.02564 |  0.03077 |
    //     |  0.35897 |  0.35897 |  0.43590 |  0.92308 |
    //     | -0.69231 | -0.69231 | -0.76923 | -1.92308 |
    #[test]
    fn calculating_inverse_of_matrix_2() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                8., 7., -6., -3., -5., 5., 0., 0., 9., 6., 9., -9., 2., 1., 6., -4.,
            ],
        );
        let expected = Matrix::from_vec(
            4,
            4,
            vec![
                -90. / 585.,
                -45. / 585.,
                210. / 585.,
                -405. / 585.,
                -90. / 585.,
                72. / 585.,
                210. / 585.,
                -405. / 585.,
                -165. / 585.,
                15. / 585.,
                255. / 585.,
                -450. / 585.,
                -315. / 585.,
                18. / 585.,
                540. / 585.,
                -1125. / 585.,
            ],
        );
        assert_eq!(matrix.inverse(), expected);
    }

    // Scenario: Calculating the inverse of a third matrix
    //   Given the following 4x4 matrix A:
    //     |  9 |  3 |  0 |  9 |
    //     | -5 | -2 | -6 | -3 |
    //     | -4 |  9 |  6 |  4 |
    //     | -7 |  6 |  6 |  2 |
    //   Then inverse(A) is the following 4x4 matrix:
    //     | -0.04074 | -0.07778 |  0.14444 | -0.22222 |
    //     | -0.07778 |  0.03333 |  0.36667 | -0.33333 |
    //     | -0.02901 | -0.14630 | -0.10926 |  0.12963 |
    //     |  0.17778 |  0.06667 | -0.26667 |  0.33333 |
    #[test]
    fn calculating_inverse_of_matrix_3() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                9., -5., -4., -7., 3., -2., 9., 6., 0., -6., 6., 6., 9., -3., 4., 2.,
            ],
        );
        let expected = Matrix::from_vec(
            4,
            4,
            vec![
                -66. / 1620.,
                -126. / 1620.,
                -47. / 1620.,
                288. / 1620.,
                -126. / 1620.,
                54. / 1620.,
                -237. / 1620.,
                108. / 1620.,
                234. / 1620.,
                594. / 1620.,
                -177. / 1620.,
                -432. / 1620.,
                -360. / 1620.,
                -540. / 1620.,
                210. / 1620.,
                540. / 1620.,
            ],
        );
        assert_eq!(matrix.inverse(), expected);
    }

    // Scenario: Multiplying a product by its inverse
    //   Given the following 4x4 matrix A:
    //       |  3 | -9 |  7 |  3 |
    //       |  3 | -8 |  2 | -9 |
    //       | -4 |  4 |  4 |  1 |
    //       | -6 |  5 | -1 |  1 |
    //     And the following 4x4 matrix B:
    //       |  8 |  2 |  2 |  2 |
    //       |  3 | -1 |  7 |  0 |
    //       |  7 |  0 |  5 |  4 |
    //       |  6 | -2 |  0 |  5 |
    //     And C ← A * B
    //   Then C * inverse(B) = A
    #[test]
    fn multiplying_product_by_inverse() {
        let matrix_a = Matrix::from_vec(
            4,
            4,
            vec![
                3., 3., -4., -6., -9., -8., 4., 5., 7., 2., 4., -1., 3., -9., 1., 1.,
            ],
        );
        let matrix_b = Matrix::from_vec(
            4,
            4,
            vec![
                8., 3., 7., 6., 2., -1., 0., -2., 2., 7., 5., 0., 2., 0., 4., 5.,
            ],
        );
        let product = &matrix_a * &matrix_b;
        let result = product * matrix_b.inverse();
        assert_relative_eq!(matrix_a, result);
    }

    // Scenario: Multiplying by a translation matrix
    //   Given transform ← translation(5, -3, 2)
    //     And p ← point(-3, 4, 5)
    //    Then transform * p = point(2, 1, 7)
    #[test]
    fn multiplying_by_translation_matrix() {
        let transform = Matrix::translation(5., -3., 2.);
        let point = Tuple::point(-3., 4., 5.);
        let expected = Tuple::point(2., 1., 7.);
        assert_eq!(transform * point, expected);
    }

    // Scenario: Multiplying by the inverse of a translation matrix
    //   Given transform ← translation(5, -3, 2)
    //     And inv ← inverse(transform)
    //     And p ← point(-3, 4, 5)
    //    Then inv * p = point(-8, 7, 3)
    #[test]
    fn multiplying_by_inverse_of_translation_matrix() {
        let transform = Matrix::translation(5., -3., 2.);
        let inv = transform.inverse();
        let point = Tuple::point(-3., 4., 5.);
        let expected = Tuple::point(-8., 7., 3.);
        assert_eq!(inv * point, expected);
    }

    // Scenario: Translation does not affect vectors
    //   Given transform ← translation(5, -3, 2)
    //     And v ← vector(-3, 4, 5)
    //    Then transform * v = v
    #[test]
    fn translation_does_not_affect_vectors() {
        let transform = Matrix::translation(5., -3., 2.);
        let vector = Tuple::vector(-3., 4., 5.);
        assert_eq!(transform * &vector, vector);
    }

    // Scenario: A scaling matrix applied to a point
    //   Given transform ← scaling(2, 3, 4)
    //     And p ← point(-4, 6, 8)
    //    Then transform * p = point(-8, 18, 32)
    #[test]
    fn scaling_matrix_applied_to_point() {
        let transform = Matrix::scaling(2., 3., 4.);
        let point = Tuple::point(-4., 6., 8.);
        let expected = Tuple::point(-8., 18., 32.);
        assert_eq!(transform * point, expected);
    }

    // Scenario: A scaling matrix applied to a vector
    //   Given transform ← scaling(2, 3, 4)
    //     And v ← vector(-4, 6, 8)
    //    Then transform * v = vector(-8, 18, 32)
    #[test]
    fn scaling_matrix_applied_to_vector() {
        let transform = Matrix::scaling(2., 3., 4.);
        let vector = Tuple::vector(-4., 6., 8.);
        let expected = Tuple::vector(-8., 18., 32.);
        assert_eq!(transform * vector, expected);
    }

    // Scenario: Multiplying by the inverse of a scaling matrix
    //   Given transform ← scaling(2, 3, 4)
    //     And inv ← inverse(transform)
    //     And v ← vector(-4, 6, 8)
    //    Then inv * v = vector(-2, 2, 2)
    #[test]
    fn multiplying_by_inverse_of_scaling_matrix() {
        let transform = Matrix::scaling(2., 3., 4.);
        let inv = transform.inverse();
        let vector = Tuple::vector(-4., 6., 8.);
        let expected = Tuple::vector(-2., 2., 2.);
        assert_eq!(inv * vector, expected);
    }

    // Scenario: Reflection is scaling by a negative value
    //   Given transform ← scaling(-1, 1, 1)
    //     And p ← point(2, 3, 4)
    //    Then transform * p = point(-2, 3, 4)
    #[test]
    fn reflection_is_scaling_by_negative_value() {
        let transform = Matrix::scaling(-1., 1., 1.);
        let point = Tuple::point(2., 3., 4.);
        let expected = Tuple::point(-2., 3., 4.);
        assert_eq!(transform * point, expected);
    }

    // Scenario: Rotating a point around the x axis
    //   Given p ← point(0, 1, 0)
    //     And half_quarter ← rotation_x(π / 4)
    //     And full_quarter ← rotation_x(π / 2)
    //   Then half_quarter * p = point(0, √2/2, √2/2)
    //     And full_quarter * p = point(0, 0, 1)
    #[test]
    fn rotating_point_around_x_axis() {
        let point = Tuple::point(0., 1., 0.);
        let half_quarter = Matrix::rotation_x(PI / 4.);
        let full_quarter = Matrix::rotation_x(PI / 2.);
        let expected_half = Tuple::point(0., SQRT_2 / 2., SQRT_2 / 2.);
        let expected_full = Tuple::point(0., 0., 1.);
        let result_half = half_quarter * &point;
        let result_full = full_quarter * point;
        assert_relative_eq!(expected_half, result_half, epsilon = 1e-5f64);
        assert_relative_eq!(expected_full, result_full, epsilon = 1e-5f64);
    }

    // Scenario: The inverse of an x-rotation rotates in the opposite direction
    //   Given p ← point(0, 1, 0)
    //     And half_quarter ← rotation_x(π / 4)
    //     And inv ← inverse(half_quarter)
    //   Then inv * p = point(0, √2/2, -√2/2)
    #[test]
    fn inverse_of_x_rotation_rotates_in_opposite_direction() {
        let point = Tuple::point(0., 1., 0.);
        let half_quarter = Matrix::rotation_x(PI / 4.);
        let inv = half_quarter.inverse();
        let expected = Tuple::point(0., SQRT_2 / 2., -SQRT_2 / 2.);
        let result = inv * point;
        assert_relative_eq!(expected, result, epsilon = 1e-5f64);
    }

    // Scenario: Rotating a point around the y axis
    //   Given p ← point(0, 0, 1)
    //     And half_quarter ← rotation_y(π / 4)
    //     And full_quarter ← rotation_y(π / 2)
    //   Then half_quarter * p = point(√2/2, 0, √2/2)
    //     And full_quarter * p = point(1, 0, 0)
    #[test]
    fn rotating_point_around_y_axis() {
        let point = Tuple::point(0., 0., 1.);
        let half_quarter = Matrix::rotation_y(PI / 4.);
        let full_quarter = Matrix::rotation_y(PI / 2.);
        let expected_half = Tuple::point(SQRT_2 / 2., 0., SQRT_2 / 2.);
        let expected_full = Tuple::point(1., 0., 0.);
        let result_half = half_quarter * &point;
        let result_full = full_quarter * point;
        assert_relative_eq!(expected_half, result_half, epsilon = 1e-5f64);
        assert_relative_eq!(expected_full, result_full, epsilon = 1e-5f64);
    }

    // Scenario: Rotating a point around the z axis
    //   Given p ← point(0, 1, 0)
    //     And half_quarter ← rotation_z(π / 4)
    //     And full_quarter ← rotation_z(π / 2)
    //   Then half_quarter * p = point(-√2/2, √2/2, 0)
    //     And full_quarter * p = point(-1, 0, 0)
    #[test]
    fn rotating_point_around_z_axis() {
        let point = Tuple::point(0., 1., 0.);
        let half_quarter = Matrix::rotation_z(PI / 4.);
        let full_quarter = Matrix::rotation_z(PI / 2.);
        let expected_half = Tuple::point(-SQRT_2 / 2., SQRT_2 / 2., 0.);
        let expected_full = Tuple::point(-1., 0., 0.);
        let result_half = half_quarter * &point;
        let result_full = full_quarter * point;
        assert_relative_eq!(expected_half, result_half, epsilon = 1e-5f64);
        assert_relative_eq!(expected_full, result_full, epsilon = 1e-5f64);
    }

    // Scenario: A shearing transformation moves x in proportion to y
    //   Given transform ← shearing(1, 0, 0, 0, 0, 0)
    //     And p ← point(2, 3, 4)
    //   Then transform * p = point(5, 3, 4)
    #[test]
    fn shearing_transformation_moves_x_in_proportion_to_y() {
        let transform = Matrix::shearing(1., 0., 0., 0., 0., 0.);
        let point = Tuple::point(2., 3., 4.);
        let expected = Tuple::point(5., 3., 4.);
        assert_eq!(transform * point, expected);
    }

    // Scenario: A shearing transformation moves x in proportion to z
    //   Given transform ← shearing(0, 1, 0, 0, 0, 0)
    //     And p ← point(2, 3, 4)
    //   Then transform * p = point(6, 3, 4)
    #[test]
    fn shearing_transformation_moves_x_in_proportion_to_z() {
        let transform = Matrix::shearing(0., 1., 0., 0., 0., 0.);
        let point = Tuple::point(2., 3., 4.);
        let expected = Tuple::point(6., 3., 4.);
        assert_eq!(transform * point, expected);
    }

    // Scenario: A shearing transformation moves y in proportion to x
    //   Given transform ← shearing(0, 0, 1, 0, 0, 0)
    //     And p ← point(2, 3, 4)
    //   Then transform * p = point(2, 5, 4)
    #[test]
    fn shearing_transformation_moves_y_in_proportion_to_x() {
        let transform = Matrix::shearing(0., 0., 1., 0., 0., 0.);
        let point = Tuple::point(2., 3., 4.);
        let expected = Tuple::point(2., 5., 4.);
        assert_eq!(transform * point, expected);
    }

    // Scenario: A shearing transformation moves y in proportion to z
    //   Given transform ← shearing(0, 0, 0, 1, 0, 0)
    //     And p ← point(2, 3, 4)
    //   Then transform * p = point(2, 7, 4)
    #[test]
    fn shearing_transformation_moves_y_in_proportion_to_z() {
        let transform = Matrix::shearing(0., 0., 0., 1., 0., 0.);
        let point = Tuple::point(2., 3., 4.);
        let expected = Tuple::point(2., 7., 4.);
        assert_eq!(transform * point, expected);
    }

    // Scenario: A shearing transformation moves z in proportion to x
    //   Given transform ← shearing(0, 0, 0, 0, 1, 0)
    //     And p ← point(2, 3, 4)
    //   Then transform * p = point(2, 3, 6)
    #[test]
    fn shearing_transformation_moves_z_in_proportion_to_x() {
        let transform = Matrix::shearing(0., 0., 0., 0., 1., 0.);
        let point = Tuple::point(2., 3., 4.);
        let expected = Tuple::point(2., 3., 6.);
        assert_eq!(transform * point, expected);
    }

    // Scenario: A shearing transformation moves z in proportion to y
    //   Given transform ← shearing(0, 0, 0, 0, 0, 1)
    //     And p ← point(2, 3, 4)
    //   Then transform * p = point(2, 3, 7)
    #[test]
    fn shearing_transformation_moves_z_in_proportion_to_y() {
        let transform = Matrix::shearing(0., 0., 0., 0., 0., 1.);
        let point = Tuple::point(2., 3., 4.);
        let expected = Tuple::point(2., 3., 7.);
        assert_eq!(transform * point, expected);
    }

    // Scenario: Individual transformations are applied in sequence
    //   Given p ← point(1, 0, 1)
    //     And A ← rotation_x(π / 2)
    //     And B ← scaling(5, 5, 5)
    //     And C ← translation(10, 5, 7)
    //   # apply rotation first
    //   When p2 ← A * p
    //   Then p2 = point(1, -1, 0)
    //   # then apply scaling
    //   When p3 ← B * p2
    //   Then p3 = point(5, -5, 0)
    //   # then apply translation
    //   When p4 ← C * p3
    //   Then p4 = point(15, 0, 7)
    #[test]
    fn individual_transformations_are_applied_in_sequence() {
        let point = Tuple::point(1., 0., 1.);
        let rotate = Matrix::rotation_x(PI / 2.);
        let scale = Matrix::scaling(5., 5., 5.);
        let translate = Matrix::translation(10., 5., 7.);
        let rotated_point = rotate * point;
        assert_relative_eq!(Tuple::point(1., -1., 0.), rotated_point, epsilon = 1e-5f64);
        let scaled_point = scale * rotated_point;
        assert_relative_eq!(Tuple::point(5., -5., 0.), scaled_point, epsilon = 1e-5f64);
        let translated_point = translate * scaled_point;
        assert_eq!(translated_point, Tuple::point(15., 0., 7.));
    }

    #[test]
    fn fluent_translate() {
        let transform = Matrix::identity(4).translate(10., 5., 7.);
        let result = Matrix::translation(10., 5., 7.);
        assert_eq!(transform, result);
    }

    #[test]
    fn fluent_scale() {
        let transform = Matrix::identity(4).scale(5., 5., 5.);
        let result = Matrix::scaling(5., 5., 5.);
        assert_eq!(transform, result);
    }

    #[test]
    fn fluent_shear() {
        let transform = Matrix::identity(4).shear(1., 2., 3., 4., 5., 6.);
        let result = Matrix::shearing(1., 2., 3., 4., 5., 6.);
        assert_eq!(transform, result);
    }

    #[test]
    fn fluent_rotate_x() {
        let transform = Matrix::identity(4).rotate_x(PI / 2.);
        let result = Matrix::rotation_x(PI / 2.);
        assert_eq!(transform, result);
    }

    #[test]
    fn fluent_rotate_y() {
        let transform = Matrix::identity(4).rotate_y(PI / 2.);
        let result = Matrix::rotation_y(PI / 2.);
        assert_eq!(transform, result);
    }

    #[test]
    fn fluent_rotate_z() {
        let transform = Matrix::identity(4).rotate_z(PI / 2.);
        let result = Matrix::rotation_z(PI / 2.);
        assert_eq!(transform, result);
    }

    // Scenario: Chained transformations must be applied in reverse order
    //   Given p ← point(1, 0, 1)
    //     And A ← rotation_x(π / 2)
    //     And B ← scaling(5, 5, 5)
    //     And C ← translation(10, 5, 7)
    //   When T ← C * B * A
    //   Then T * p = point(15, 0, 7)
    #[test]
    fn chained_transformations_must_be_applied_in_reverse_order() {
        let point = Tuple::point(1., 0., 1.);
        let rotate = Matrix::rotation_x(PI / 2.);
        let scale = Matrix::scaling(5., 5., 5.);
        let translate = Matrix::translation(10., 5., 7.);
        let transformation = translate * scale * rotate;
        assert_eq!(transformation * point, Tuple::point(15., 0., 7.));
    }

    #[test]
    fn portrait_matrix() {
        let mut matrix = Matrix::new(2, 5, 0.);
        matrix.set(0, 0, 0.);
        matrix.set(0, 1, 1.);
        matrix.set(0, 2, 2.);
        matrix.set(0, 3, 3.);
        matrix.set(0, 4, 4.);
        matrix.set(1, 0, 5.);
        matrix.set(1, 1, 6.);
        matrix.set(1, 2, 7.);
        matrix.set(1, 3, 8.);
        matrix.set(1, 4, 9.);
        let expected = Matrix::from_vec(2, 5, vec![0., 1., 2., 3., 4., 5., 6., 7., 8., 9.]);
        assert_eq!(matrix, expected);
    }

    #[test]
    fn landscape_matrix() {
        let mut matrix = Matrix::new(5, 2, 0.);
        matrix.set(0, 0, 0.);
        matrix.set(0, 1, 1.);
        matrix.set(1, 0, 2.);
        matrix.set(1, 1, 3.);
        matrix.set(2, 0, 4.);
        matrix.set(2, 1, 5.);
        matrix.set(3, 0, 6.);
        matrix.set(3, 1, 7.);
        matrix.set(4, 0, 8.);
        matrix.set(4, 1, 9.);
        let expected = Matrix::from_vec(5, 2, vec![0., 1., 2., 3., 4., 5., 6., 7., 8., 9.]);
        assert_eq!(matrix, expected);
    }
}
