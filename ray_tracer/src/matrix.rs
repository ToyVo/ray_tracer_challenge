use std::ops::Mul;
use approx::{AbsDiffEq, RelativeEq, UlpsEq};

use crate::Tuple;

#[derive(Clone, Debug, PartialEq)]
pub struct Matrix {
    data: Vec<f64>,
    rows: usize,
    cols: usize,
}

impl Matrix {
    pub fn new(rows: usize, cols: usize, init: f64) -> Matrix {
        let data = vec![init; rows * cols];
        Matrix { data, rows, cols }
    }

    pub fn from_vec(rows: usize, cols: usize, data: Vec<f64>) -> Matrix {
        Matrix { data, rows, cols }
    }

    pub fn get(&self, row: usize, col: usize) -> f64 {
        assert!(row < self.rows && col < self.cols);
        self.data[row * self.cols + col]
    }

    pub fn set(&mut self, row: usize, col: usize, value: f64) {
        assert!(row < self.rows && col < self.cols);
        self.data[row * self.cols + col] = value;
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn width(&self) -> usize {
        self.cols
    }

    pub fn height(&self) -> usize {
        self.rows
    }

    pub fn get_row(&self, row: usize) -> Tuple {
        assert!(row < self.rows);
        let mut result = Tuple::new(self.cols, 0.);
        for col in 0..self.cols {
            result.set(col, self.get(row, col));
        }
        result
    }

    pub fn get_col(&self, col: usize) -> Tuple {
        assert!(col < self.cols);
        let mut result = Tuple::new(self.rows, 0.);
        for row in 0..self.rows {
            result.set(row, self.get(row, col));
        }
        result
    }

    pub fn identity(size: usize) -> Matrix {
        let mut result = Matrix::new(size, size, 0.);
        for i in 0..size {
            result.set(i, i, 1.);
        }
        result
    }

    pub fn translation(x: f64, y: f64, z: f64) -> Matrix {
        let mut result = Matrix::identity(4);
        result.set(0, 3, x);
        result.set(1, 3, y);
        result.set(2, 3, z);
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
        result.set(0, 1, xy);
        result.set(0, 2, xz);
        result.set(1, 0, yx);
        result.set(1, 2, yz);
        result.set(2, 0, zx);
        result.set(2, 1, zy);
        result
    }

    pub fn transpose(&self) -> Matrix {
        let mut result = Matrix::new(self.cols, self.rows, 0.);
        for row in 0..self.rows {
            for col in 0..self.cols {
                result.set(col, row, self.get(row, col));
            }
        }
        result
    }

    pub fn submatrix(&self, row: usize, col: usize) -> Matrix {
        assert_eq!(self.rows, self.cols);
        let mut result = Matrix::new(self.rows - 1, self.cols - 1, 0.);
        for i in 0..self.rows {
            for j in 0..self.cols {
                if i != row && j != col {
                    let new_i = if i > row { i - 1 } else { i };
                    let new_j = if j > col { j - 1 } else { j };
                    result.set(new_i, new_j, self.get(i, j));
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

    pub fn determinant(&self) -> f64 {
        assert_eq!(self.rows, self.cols);
        if self.cols == 2 {
            self.get(0, 0) * self.get(1, 1) - self.get(0, 1) * self.get(1, 0)
        } else {
            let mut det = 0.;
            for col in 0..self.cols {
                det = det + self.get(0, col) * self.cofactor(0, col);
            }
            det
        }
    }

    pub fn is_invertible(&self) -> bool {
        self.determinant() != 0.
    }

    pub fn minor(&self, row: usize, col: usize) -> f64 {
        self.submatrix(row, col).determinant()
    }

    pub fn cofactor(&self, row: usize, col: usize) -> f64 {
        let minor = self.minor(row, col);
        if (row + col) % 2 == 0 {
            minor
        } else {
            -minor
        }
    }

    pub fn inverse(&self) -> Matrix {
        let mut result = Matrix::new(self.rows, self.cols, 0.);
        let det = self.determinant();
        for row in 0..self.rows {
            for col in 0..self.cols {
                result.set(col, row, self.cofactor(row, col) / det);
            }
        }
        result
    }

    pub fn rotation_x(rad: f64) -> Matrix {
        let mut result = Matrix::identity(4);
        result.set(1, 1, f64::cos(rad));
        result.set(1, 2, -f64::sin(rad));
        result.set(2, 1, f64::sin(rad));
        result.set(2, 2, f64::cos(rad));
        result
    }

    pub fn rotation_y(rad: f64) -> Matrix {
        let mut result = Matrix::identity(4);
        result.set(0, 0, f64::cos(rad));
        result.set(0, 2, f64::sin(rad));
        result.set(2, 0, -f64::sin(rad));
        result.set(2, 2, f64::cos(rad));
        result
    }

    pub fn rotation_z(rad: f64) -> Matrix {
        let mut result = Matrix::identity(4);
        result.set(0, 0, f64::cos(rad));
        result.set(0, 1, -f64::sin(rad));
        result.set(1, 0, f64::sin(rad));
        result.set(1, 1, f64::cos(rad));
        result
    }

    pub fn rotate_x(&self, rad: f64) -> Matrix {
        Matrix::rotation_x(rad) * self
    }

    pub fn rotate_y(&self, rad: f64) -> Matrix {
        Matrix::rotation_y(rad) * self
    }

    pub fn rotate_z(&self, rad: f64) -> Matrix {
        Matrix::rotation_z(rad) * self
    }
}

impl Mul<&Matrix> for &Matrix {
    type Output = Matrix;

    fn mul(self, other: &Matrix) -> Matrix {
        assert_eq!(self.rows, other.cols);
        let mut result = Matrix::new(other.rows, self.cols, 0.);
        for row in 0..other.rows {
            for col in 0..self.cols {
                let tuple_row = self.get_row(row);
                let tuple_col = other.get_col(col);
                result.set(row, col, tuple_col.dot(&tuple_row));
            }
        }
        result
    }
}

impl Mul<Matrix> for Matrix {
    type Output = Matrix;

    fn mul(self, other: Matrix) -> Matrix {
        &self * &other
    }
}

impl Mul<&Matrix> for Matrix {
    type Output = Matrix;

    fn mul(self, other: &Matrix) -> Matrix {
        &self * other
    }
}

impl Mul<Matrix> for &Matrix {
    type Output = Matrix;

    fn mul(self, other: Matrix) -> Matrix {
        self * &other
    }
}

impl Mul<&Tuple> for &Matrix {
    type Output = Tuple;

    fn mul(self, other: &Tuple) -> Tuple {
        assert_eq!(self.rows, other.dimension());
        let mut result = Tuple::new(self.cols(), 0.);
        for row in 0..self.rows {
            let tuple_row = self.get_row(row);
            result.set(row, other.dot(&tuple_row));
        }
        result
    }
}

impl Mul<Tuple> for Matrix {
    type Output = Tuple;

    fn mul(self, other: Tuple) -> Tuple {
        &self * &other
    }
}

impl Mul<&Tuple> for Matrix {
    type Output = Tuple;

    fn mul(self, other: &Tuple) -> Tuple {
        &self * other
    }
}

impl Mul<Tuple> for &Matrix {
    type Output = Tuple;

    fn mul(self, other: Tuple) -> Tuple {
        self * &other
    }
}

impl AbsDiffEq for Matrix {
    type Epsilon = f64;

    fn default_epsilon() -> f64 {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Matrix, epsilon: f64) -> bool {
        if self.cols != other.cols || self.rows != other.rows {
            return false;
        }
        for x in 0..self.cols {
            for y in 0..self.rows {
                if !self.get(x, y).abs_diff_eq(&other.get(x, y), epsilon) {
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

    fn relative_eq(
        &self,
        other: &Matrix,
        epsilon: f64,
        max_relative: f64,
    ) -> bool {
        if self.cols != other.cols || self.rows != other.rows {
            return false;
        }
        for x in 0..self.cols {
            for y in 0..self.rows {
                if !self.get(x, y).relative_eq(&other.get(x, y), epsilon, max_relative) {
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
                if !self.get(x, y).ulps_eq(&other.get(x, y), epsilon, max_ulps) {
                    return false;
                }
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;
    use std::f64::consts::SQRT_2;
    use approx::assert_relative_eq;

    use super::*;

    #[test]
    fn constructing_4x4_matrix() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                1., 2., 3., 4., 5.5, 6.5, 7.5, 8.5, 9., 10., 11., 12., 13.5, 14.5, 15.5, 16.5,
            ],
        );
        assert_eq!(matrix.get(0, 0), 1.);
        assert_eq!(matrix.get(0, 3), 4.);
        assert_eq!(matrix.get(1, 0), 5.5);
        assert_eq!(matrix.get(1, 2), 7.5);
        assert_eq!(matrix.get(2, 2), 11.);
        assert_eq!(matrix.get(3, 0), 13.5);
        assert_eq!(matrix.get(3, 2), 15.5);
    }

    #[test]
    fn constructing_2x2_matrix() {
        let matrix = Matrix::from_vec(2, 2, vec![-3., 5., 1., -2.]);
        assert_eq!(matrix.get(0, 0), -3.);
        assert_eq!(matrix.get(0, 1), 5.);
        assert_eq!(matrix.get(1, 0), 1.);
        assert_eq!(matrix.get(1, 1), -2.);
    }

    #[test]
    fn construction_3x3_matrix() {
        let matrix = Matrix::from_vec(3, 3, vec![-3., 5., 0., 1., -2., -7., 0., 1., 1.]);
        assert_eq!(matrix.get(0, 0), -3.);
        assert_eq!(matrix.get(1, 1), -2.);
        assert_eq!(matrix.get(2, 2), 1.);
    }

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

    #[test]
    fn multiplying_two_matrices() {
        let matrix_a = Matrix::from_vec(
            4,
            4,
            vec![
                1., 2., 3., 4., 5., 6., 7., 8., 9., 8., 7., 6., 5., 4., 3., 2.,
            ],
        );
        let matrix_b = Matrix::from_vec(
            4,
            4,
            vec![
                -2., 1., 2., 3., 3., 2., 1., -1., 4., 3., 6., 5., 1., 2., 7., 8.,
            ],
        );
        let expected = Matrix::from_vec(
            4,
            4,
            vec![
                20., 22., 50., 48., 44., 54., 114., 108., 40., 58., 110., 102., 16., 26., 46., 42.,
            ],
        );
        assert_eq!(matrix_a * matrix_b, expected);
    }

    #[test]
    fn multiplying_matrix_by_tuple() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                1., 2., 3., 4., 2., 4., 4., 2., 8., 6., 4., 1., 0., 0., 0., 1.,
            ],
        );
        let point = Tuple::from_vec(vec![1., 2., 3., 1.]);
        let expected = Tuple::from_vec(vec![18., 24., 33., 1.]);
        assert_eq!(matrix * point, expected);
    }

    #[test]
    fn get_row_returns_correct_row() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                1., 2., 3., 4., 5., 6., 7., 8., 9., 0., 1., 2., 3., 4., 5., 6.,
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
                1., 2., 3., 4., 5., 6., 7., 8., 9., 0., 1., 2., 3., 4., 5., 6.,
            ],
        );
        let expected = Tuple::from_vec(vec![1., 5., 9., 3.]);
        assert_eq!(matrix.get_col(0), expected);
    }

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

    #[test]
    fn transposing_identity_matrix() {
        let matrix = Matrix::identity(4);
        assert_eq!(matrix.transpose(), matrix);
    }

    #[test]
    fn calculating_determinant_of_2x2_matrix() {
        let matrix = Matrix::from_vec(2, 2, vec![1., 5., -3., 2.]);
        assert_eq!(matrix.determinant(), 17.);
    }

    #[test]
    fn submatrix_of_3x3_is_2x2() {
        let matrix = Matrix::from_vec(3, 3, vec![1., 5., 0., -3., 2., 7., 0., 6., -3.]);
        let expected = Matrix::from_vec(2, 2, vec![-3., 2., 0., 6.]);
        assert_eq!(matrix.submatrix(0, 2), expected);
    }

    #[test]
    fn submatrix_of_4x4_is_3x3() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                -6., 1., 1., 6., -8., 5., 8., 6., -1., 0., 8., 2., -7., 1., -1., 1.,
            ],
        );
        let expected = Matrix::from_vec(3, 3, vec![-6., 1., 6., -8., 8., 6., -7., -1., 1.]);
        assert_eq!(matrix.submatrix(2, 1), expected);
    }

    #[test]
    fn calculating_minor_of_3x3_matrix() {
        let matrix = Matrix::from_vec(3, 3, vec![1., 5., 0., 2., -1., -7., 6., -1., 5.]);
        let submatrix = matrix.submatrix(1, 0);
        assert_eq!(submatrix.determinant(), 25.);
        assert_eq!(matrix.minor(1, 0), 25.);
    }

    #[test]
    fn calculating_cofactor_of_3x3_matrix() {
        let matrix = Matrix::from_vec(3, 3, vec![3., 5., 0., 2., -1., -7., 6., -1., 5.]);
        assert_eq!(matrix.minor(0, 0), -12.);
        assert_eq!(matrix.cofactor(0, 0), -12.);
        assert_eq!(matrix.minor(1, 0), 25.);
        assert_eq!(matrix.cofactor(1, 0), -25.);
    }

    #[test]
    fn calculating_determinant_of_3x3_matrix() {
        let matrix = Matrix::from_vec(3, 3, vec![1., 2., 6., -5., 8., -4., 2., 6., 4.]);
        assert_eq!(matrix.cofactor(0, 0), 56.);
        assert_eq!(matrix.cofactor(0, 1), 12.);
        assert_eq!(matrix.cofactor(0, 2), -46.);
        assert_eq!(matrix.determinant(), -196.);
    }

    #[test]
    fn calculating_determinant_of_4x4_matrix() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                -2., -8., 3., 5., -3., 1., 7., 3., 1., 2., -9., 6., -6., 7., 7., -9.,
            ],
        );
        assert_eq!(matrix.cofactor(0, 0), 690.);
        assert_eq!(matrix.cofactor(0, 1), 447.);
        assert_eq!(matrix.cofactor(0, 2), 210.);
        assert_eq!(matrix.cofactor(0, 3), 51.);
        assert_eq!(matrix.determinant(), -4071.);
    }

    #[test]
    fn testing_invertible_matrix_for_invertibility() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                6., 4., 4., 4., 5., 5., 7., 6., 4., -9., 3., -7., 9., 1., 7., -6.,
            ],
        );
        assert_eq!(matrix.determinant(), -2120.);
        assert!(matrix.is_invertible());
    }

    #[test]
    fn testing_non_invertible_matrix_for_invertibility() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                -4., 2., -2., -3., 9., 6., 2., 6., 0., -5., 1., -5., 0., 0., 0., 0.,
            ],
        );
        assert_eq!(matrix.determinant(), 0.);
        assert!(!matrix.is_invertible());
    }

    #[test]
    fn calculating_inverse_of_matrix() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                -5., 2., 6., -8., 1., -5., 1., 8., 7., 7., -6., -7., 1., -3., 7., 4.,
            ],
        );
        let inv = matrix.inverse();
        let expected = Matrix::from_vec(
            4,
            4,
            vec![
                116. / 532.,
                240. / 532.,
                128. / 532.,
                -24. / 532.,
                -430. / 532.,
                -775. / 532.,
                -236. / 532.,
                277. / 532.,
                -42. / 532.,
                -119. / 532.,
                -28. / 532.,
                105. / 532.,
                -278. / 532.,
                -433. / 532.,
                -160. / 532.,
                163. / 532.,
            ],
        );
        assert_eq!(matrix.determinant(), 532.);
        assert_eq!(matrix.cofactor(2, 3), -160.);
        assert_eq!(inv.get(3, 2), -160. / 532.);
        assert_eq!(matrix.cofactor(3, 2), 105.);
        assert_eq!(inv.get(2, 3), 105. / 532.);
        assert_eq!(inv, expected);
    }

    #[test]
    fn calculating_inverse_of_matrix_2() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                8., -5., 9., 2., 7., 5., 6., 1., -6., 0., 9., 6., -3., 0., -9., -4.,
            ],
        );
        let expected = Matrix::from_vec(
            4,
            4,
            vec![
                -90. / 585.,
                -90. / 585.,
                -165. / 585.,
                -315. / 585.,
                -45. / 585.,
                72. / 585.,
                15. / 585.,
                18. / 585.,
                210. / 585.,
                210. / 585.,
                255. / 585.,
                540. / 585.,
                -405. / 585.,
                -405. / 585.,
                -450. / 585.,
                -1125. / 585.,
            ],
        );
        assert_eq!(matrix.inverse(), expected);
    }

    #[test]
    fn calculating_inverse_of_matrix_3() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                9., 3., 0., 9., -5., -2., -6., -3., -4., 9., 6., 4., -7., 6., 6., 2.,
            ],
        );
        let expected = Matrix::from_vec(
            4,
            4,
            vec![
                -66. / 1620.,
                -126. / 1620.,
                234. / 1620.,
                -360. / 1620.,
                -126. / 1620.,
                54. / 1620.,
                594. / 1620.,
                -540. / 1620.,
                -47. / 1620.,
                -237. / 1620.,
                -177. / 1620.,
                210. / 1620.,
                288. / 1620.,
                108. / 1620.,
                -432. / 1620.,
                540. / 1620.,
            ],
        );
        assert_eq!(matrix.inverse(), expected);
    }

    #[test]
    fn multiplying_product_by_inverse() {
        let matrix_a = Matrix::from_vec(
            4,
            4,
            vec![
                3., -9., 7., 3., 3., -8., 2., -9., -4., 4., 4., 1., -6., 5., -1., 1.,
            ],
        );
        let matrix_b = Matrix::from_vec(
            4,
            4,
            vec![
                8., 2., 2., 2., 3., -1., 7., 0., 7., 0., 5., 4., 6., -2., 0., 5.,
            ],
        );
        let product = &matrix_a * &matrix_b;
        let result = product * matrix_b.inverse();
        assert_relative_eq!(matrix_a, result, epsilon = 1e-5f64);
    }

    #[test]
    fn multiplying_by_translation_matrix() {
        let transform = Matrix::translation(5., -3., 2.);
        let point = Tuple::point(-3., 4., 5.);
        let expected = Tuple::point(2., 1., 7.);
        assert_eq!(transform * point, expected);
    }

    #[test]
    fn multiplying_by_inverse_of_translation_matrix() {
        let transform = Matrix::translation(5., -3., 2.);
        let inv = transform.inverse();
        let point = Tuple::point(-3., 4., 5.);
        let expected = Tuple::point(-8., 7., 3.);
        assert_eq!(inv * point, expected);
    }

    #[test]
    fn translation_does_not_affect_vectors() {
        let transform = Matrix::translation(5., -3., 2.);
        let vector = Tuple::vector(-3., 4., 5.);
        assert_eq!(transform * &vector, vector);
    }

    #[test]
    fn scaling_matrix_applied_to_point() {
        let transform = Matrix::scaling(2., 3., 4.);
        let point = Tuple::point(-4., 6., 8.);
        let expected = Tuple::point(-8., 18., 32.);
        assert_eq!(transform * point, expected);
    }

    #[test]
    fn scaling_matrix_applied_to_vector() {
        let transform = Matrix::scaling(2., 3., 4.);
        let vector = Tuple::vector(-4., 6., 8.);
        let expected = Tuple::vector(-8., 18., 32.);
        assert_eq!(transform * vector, expected);
    }

    #[test]
    fn multiplying_by_inverse_of_scaling_matrix() {
        let transform = Matrix::scaling(2., 3., 4.);
        let inv = transform.inverse();
        let vector = Tuple::vector(-4., 6., 8.);
        let expected = Tuple::vector(-2., 2., 2.);
        assert_eq!(inv * vector, expected);
    }

    #[test]
    fn reflection_is_scaling_by_negative_value() {
        let transform = Matrix::scaling(-1., 1., 1.);
        let point = Tuple::point(2., 3., 4.);
        let expected = Tuple::point(-2., 3., 4.);
        assert_eq!(transform * point, expected);
    }

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

    #[test]
    fn inverse_of_x_rotation_rotates_in_opposite_direction() {
        let point = Tuple::point(0., 1., 0.);
        let half_quarter = Matrix::rotation_x(PI / 4.);
        let inv = half_quarter.inverse();
        let expected = Tuple::point(0., SQRT_2 / 2., -SQRT_2 / 2.);
        let result = inv * point;
        assert_relative_eq!(expected, result, epsilon = 1e-5f64);
    }

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

    #[test]
    fn shearing_transformation_moves_x_in_proportion_to_y() {
        let transform = Matrix::shearing(1., 0., 0., 0., 0., 0.);
        let point = Tuple::point(2., 3., 4.);
        let expected = Tuple::point(5., 3., 4.);
        assert_eq!(transform * point, expected);
    }

    #[test]
    fn shearing_transformation_moves_x_in_proportion_to_z() {
        let transform = Matrix::shearing(0., 1., 0., 0., 0., 0.);
        let point = Tuple::point(2., 3., 4.);
        let expected = Tuple::point(6., 3., 4.);
        assert_eq!(transform * point, expected);
    }

    #[test]
    fn shearing_transformation_moves_y_in_proportion_to_x() {
        let transform = Matrix::shearing(0., 0., 1., 0., 0., 0.);
        let point = Tuple::point(2., 3., 4.);
        let expected = Tuple::point(2., 5., 4.);
        assert_eq!(transform * point, expected);
    }

    #[test]
    fn shearing_transformation_moves_y_in_proportion_to_z() {
        let transform = Matrix::shearing(0., 0., 0., 1., 0., 0.);
        let point = Tuple::point(2., 3., 4.);
        let expected = Tuple::point(2., 7., 4.);
        assert_eq!(transform * point, expected);
    }

    #[test]
    fn shearing_transformation_moves_z_in_proportion_to_x() {
        let transform = Matrix::shearing(0., 0., 0., 0., 1., 0.);
        let point = Tuple::point(2., 3., 4.);
        let expected = Tuple::point(2., 3., 6.);
        assert_eq!(transform * point, expected);
    }

    #[test]
    fn shearing_transformation_moves_z_in_proportion_to_y() {
        let transform = Matrix::shearing(0., 0., 0., 0., 0., 1.);
        let point = Tuple::point(2., 3., 4.);
        let expected = Tuple::point(2., 3., 7.);
        assert_eq!(transform * point, expected);
    }

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

    #[test]
    fn landscape_matrix() {
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
}
