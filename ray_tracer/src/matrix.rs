use approx::{AbsDiffEq, RelativeEq, UlpsEq};
use std::ops::Mul;

use crate::Tuple;

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

impl Mul<&Matrix> for &Matrix {
    type Output = Matrix;

    fn mul(self, other: &Matrix) -> Matrix {
        assert_eq!(self.cols, other.rows);
        let mut result = Matrix::new(self.cols, other.rows, 0.);
        for col in 0..self.cols {
            for row in 0..other.rows {
                let tuple_row = self.get_row(row);
                let tuple_col = other.get_col(col);
                result.set(col, row, tuple_col.dot(&tuple_row));
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

    #[test]
    fn constructing_2x2_matrix() {
        let matrix = Matrix::from_vec(2, 2, vec![-3., 1., 5., -2.]);
        assert_eq!(matrix.get(0, 0), -3.);
        assert_eq!(matrix.get(1, 0), 5.);
        assert_eq!(matrix.get(0, 1), 1.);
        assert_eq!(matrix.get(1, 1), -2.);
    }

    #[test]
    fn construction_3x3_matrix() {
        let matrix = Matrix::from_vec(3, 3, vec![-3., 1., 0., 5., -2., 1., 0., -7., 1.]);
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

    #[test]
    fn multiplying_matrix_by_tuple() {
        let matrix = Matrix::from_vec(
            4,
            4,
            vec![
                1., 2., 8., 0., 2., 4., 6., 0., 3., 4., 4., 0., 4., 2., 1., 1.,
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
        let matrix = Matrix::from_vec(2, 2, vec![1., -3., 5., 2.]);
        assert_eq!(matrix.determinant(), 17.);
    }

    #[test]
    fn submatrix_of_3x3_is_2x2() {
        let matrix = Matrix::from_vec(3, 3, vec![1., -3., 0., 5., 2., 6., 0., 7., -3.]);
        let expected = Matrix::from_vec(2, 2, vec![-3., 0., 2., 6.]);
        assert_eq!(matrix.submatrix(2, 0), expected);
    }

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

    #[test]
    fn calculating_minor_of_3x3_matrix() {
        let matrix = Matrix::from_vec(3, 3, vec![1., 2., 6., 5., -1., -1., 0., -7., 5.]);
        let submatrix = matrix.submatrix(0, 1);
        assert_eq!(submatrix.determinant(), 25.);
        assert_eq!(matrix.minor(0, 1), 25.);
    }

    #[test]
    fn calculating_cofactor_of_3x3_matrix() {
        let matrix = Matrix::from_vec(3, 3, vec![3., 2., 6., 5., -1., -1., 0., -7., 5.]);
        assert_eq!(matrix.minor(0, 0), -12.);
        assert_eq!(matrix.cofactor(0, 0), -12.);
        assert_eq!(matrix.minor(0, 1), 25.);
        assert_eq!(matrix.cofactor(0, 1), -25.);
    }

    #[test]
    fn calculating_determinant_of_3x3_matrix() {
        let matrix = Matrix::from_vec(3, 3, vec![1., -5., 2., 2., 8., 6., 6., -4., 4.]);
        assert_eq!(matrix.cofactor(0, 0), 56.);
        assert_eq!(matrix.cofactor(1, 0), 12.);
        assert_eq!(matrix.cofactor(2, 0), -46.);
        assert_eq!(matrix.determinant(), -196.);
    }

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
