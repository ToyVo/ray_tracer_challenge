use crate::Tuple;
use num::Num;
use std::iter::Sum;
use std::ops::{Mul, Neg};

#[derive(Clone, Debug, PartialEq)]
pub struct Matrix<T> {
    data: Vec<T>,
    rows: usize,
    cols: usize,
}

impl<T> Matrix<T>
where
    T: Clone,
{
    pub fn new(rows: usize, cols: usize, init: T) -> Matrix<T> {
        let data = vec![init; rows * cols];
        Matrix { data, rows, cols }
    }

    pub fn from_vec(rows: usize, cols: usize, data: Vec<T>) -> Matrix<T> {
        Matrix { data, rows, cols }
    }

    pub fn get(&self, row: usize, col: usize) -> T {
        assert!(row < self.rows && col < self.cols);
        self.data[row * self.cols + col].clone()
    }

    pub fn set(&mut self, row: usize, col: usize, value: T) {
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
}

impl<T> Matrix<T>
where
    T: Num + Copy + Clone,
{
    pub fn get_row(&self, row: usize) -> Tuple<T> {
        assert!(row < self.rows);
        let mut result = Tuple::new(self.cols, T::zero());
        for col in 0..self.cols {
            result.set(col, self.get(row, col));
        }
        result
    }

    pub fn get_col(&self, col: usize) -> Tuple<T> {
        assert!(col < self.cols);
        let mut result = Tuple::new(self.rows, T::zero());
        for row in 0..self.rows {
            result.set(row, self.get(row, col));
        }
        result
    }

    pub fn identity(size: usize) -> Matrix<T> {
        let mut result = Matrix::new(size, size, T::zero());
        for i in 0..size {
            result.set(i, i, T::one());
        }
        result
    }

    pub fn transpose(&self) -> Matrix<T> {
        let mut result = Matrix::new(self.cols, self.rows, T::zero());
        for row in 0..self.rows {
            for col in 0..self.cols {
                result.set(col, row, self.get(row, col));
            }
        }
        result
    }

    pub fn submatrix(&self, row: usize, col: usize) -> Matrix<T> {
        assert_eq!(self.rows, self.cols);
        let mut result = Matrix::new(self.rows - 1, self.cols - 1, T::zero());
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
}

impl<T> Matrix<T>
where
    T: Num + Copy + Clone + Neg<Output = T>,
{
    pub fn determinant(&self) -> T {
        assert_eq!(self.rows, self.cols);
        if self.cols == 2 {
            self.get(0, 0) * self.get(1, 1) - self.get(0, 1) * self.get(1, 0)
        } else {
            let mut det = T::zero();
            for col in 0..self.cols {
                det = det + self.get(0, col) * self.cofactor(0, col);
            }
            det
        }
    }

    pub fn is_invertible(&self) -> bool {
        self.determinant() != T::zero()
    }

    pub fn minor(&self, row: usize, col: usize) -> T {
        self.submatrix(row, col).determinant()
    }

    pub fn cofactor(&self, row: usize, col: usize) -> T {
        let minor = self.minor(row, col);
        if (row + col) % 2 == 0 {
            minor
        } else {
            -minor
        }
    }
}

impl Matrix<f64> {
    pub fn inverse(&self) -> Self {
        let mut result = Matrix::new(self.rows, self.cols, 0.);
        let det = self.determinant();
        for row in 0..self.rows {
            for col in 0..self.cols {
                result.set(col, row, self.cofactor(row, col) / det);
            }
        }
        result
    }

    pub fn equals(&self, other: &Self, epsilon: f64) -> bool {
        for row in 0..self.rows {
            for col in 0..self.cols {
                if (self.get(row, col) - other.get(row, col)).abs() > epsilon {
                    return false;
                }
            }
        }
        true
    }
}

impl<T> Mul<Matrix<T>> for Matrix<T>
where
    T: Num + Copy + Clone + Sum,
{
    type Output = Matrix<T>;

    fn mul(self, other: Matrix<T>) -> Matrix<T> {
        assert_eq!(self.rows, other.cols);
        let mut result = Matrix::new(other.rows, self.cols, T::zero());
        for row in 0..other.rows {
            for col in 0..self.cols {
                let tuple_row = self.get_row(row);
                let tuple_col = other.get_col(col);
                result.set(row, col, tuple_col.dot(tuple_row));
            }
        }
        result
    }
}

impl<T> Mul<Tuple<T>> for Matrix<T>
where
    T: Num + Copy + Clone + Sum,
{
    type Output = Tuple<T>;

    fn mul(self, other: Tuple<T>) -> Tuple<T> {
        assert_eq!(self.rows, other.dimension());
        let mut result = Tuple::new(self.cols(), T::zero());
        for row in 0..self.rows {
            let tuple_row = self.get_row(row);
            result.set(row, other.dot(tuple_row));
        }
        result
    }
}

impl std::fmt::Display for Matrix<Tuple<f64>> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut ppm = String::with_capacity(self.cols() * self.rows() * 12);
        ppm.push_str("P3\n");
        ppm.push_str(&format!("{} {}\n", self.cols(), self.rows()));
        ppm.push_str("255\n");

        for row in 0..self.rows() {
            let mut row_colors: Vec<u8> = Vec::with_capacity(self.cols() * 3);
            for col in 0..self.cols() {
                let color = self.get(row, col);
                row_colors.push((color.r() * 255.0).round() as u8);
                row_colors.push((color.g() * 255.0).round() as u8);
                row_colors.push((color.b() * 255.0).round() as u8);
            }
            let mut row = String::new();
            let mut line_length = 0;
            for (_i, color) in row_colors.iter().enumerate() {
                let color_str = format!("{} ", color);
                let color_len = color_str.len();
                if line_length + color_len > 70 {
                    row = row.trim_end().to_string();
                    row.push('\n');
                    line_length = 0;
                }
                row.push_str(&color_str);
                line_length += color_len;
            }
            ppm.push_str(row.trim_end());
            ppm.push('\n');
        }
        write!(f, "{}", ppm)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constructing_4x4_matrix() {
        let m = Matrix::from_vec(
            4,
            4,
            vec![
                1., 2., 3., 4., 5.5, 6.5, 7.5, 8.5, 9., 10., 11., 12., 13.5, 14.5, 15.5, 16.5,
            ],
        );
        assert_eq!(m.get(0, 0), 1.);
        assert_eq!(m.get(0, 3), 4.);
        assert_eq!(m.get(1, 0), 5.5);
        assert_eq!(m.get(1, 2), 7.5);
        assert_eq!(m.get(2, 2), 11.);
        assert_eq!(m.get(3, 0), 13.5);
        assert_eq!(m.get(3, 2), 15.5);
    }

    #[test]
    fn constructing_2x2_matrix() {
        let m = Matrix::from_vec(2, 2, vec![-3, 5, 1, -2]);
        assert_eq!(m.get(0, 0), -3);
        assert_eq!(m.get(0, 1), 5);
        assert_eq!(m.get(1, 0), 1);
        assert_eq!(m.get(1, 1), -2);
    }

    #[test]
    fn construction_3x3_matrix() {
        let m = Matrix::from_vec(3, 3, vec![-3, 5, 0, 1, -2, -7, 0, 1, 1]);
        assert_eq!(m.get(0, 0), -3);
        assert_eq!(m.get(1, 1), -2);
        assert_eq!(m.get(2, 2), 1);
    }

    #[test]
    fn matrix_equality_with_identical_matrices() {
        let a = Matrix::from_vec(4, 4, vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 8, 7, 6, 5, 4, 3, 2]);
        let b = Matrix::from_vec(4, 4, vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 8, 7, 6, 5, 4, 3, 2]);
        assert_eq!(a, b);
    }

    #[test]
    fn matrix_equality_with_different_matrices() {
        let a = Matrix::from_vec(4, 4, vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 8, 7, 6, 5, 4, 3, 2]);
        let b = Matrix::from_vec(4, 4, vec![2, 3, 4, 5, 6, 7, 8, 9, 8, 7, 6, 5, 4, 3, 2, 1]);
        assert_ne!(a, b);
    }

    #[test]
    fn multiplying_two_matrices() {
        let a = Matrix::from_vec(4, 4, vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 8, 7, 6, 5, 4, 3, 2]);
        let b = Matrix::from_vec(4, 4, vec![-2, 1, 2, 3, 3, 2, 1, -1, 4, 3, 6, 5, 1, 2, 7, 8]);
        let expected = Matrix::from_vec(
            4,
            4,
            vec![
                20, 22, 50, 48, 44, 54, 114, 108, 40, 58, 110, 102, 16, 26, 46, 42,
            ],
        );
        assert_eq!(a * b, expected);
    }

    #[test]
    fn multiplying_matrix_by_tuple() {
        let a = Matrix::from_vec(
            4,
            4,
            vec![
                1., 2., 3., 4., 2., 4., 4., 2., 8., 6., 4., 1., 0., 0., 0., 1.,
            ],
        );
        let b = Tuple::from_vec(vec![1., 2., 3., 1.]);
        let expected = Tuple::from_vec(vec![18., 24., 33., 1.]);
        assert_eq!(a * b, expected);
    }

    #[test]
    fn get_row_returns_correct_row() {
        let a = Matrix::from_vec(4, 4, vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6]);
        let expected = Tuple::from_vec(vec![1, 2, 3, 4]);
        assert_eq!(a.get_row(0), expected);
    }

    #[test]
    fn get_col_returns_correct_col() {
        let a = Matrix::from_vec(4, 4, vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6]);
        let expected = Tuple::from_vec(vec![1, 5, 9, 3]);
        assert_eq!(a.get_col(0), expected);
    }

    #[test]
    fn multiplying_matrix_by_identity_matrix() {
        let a = Matrix::from_vec(
            4,
            4,
            vec![0, 1, 2, 4, 1, 2, 4, 8, 2, 4, 8, 16, 4, 8, 16, 32],
        );
        assert_eq!(a.clone() * Matrix::identity(4), a);
    }

    #[test]
    fn transposing_matrix() {
        let a = Matrix::from_vec(4, 4, vec![0, 9, 3, 0, 9, 8, 0, 8, 1, 8, 5, 3, 0, 0, 5, 8]);
        let expected = Matrix::from_vec(4, 4, vec![0, 9, 1, 0, 9, 8, 8, 0, 3, 0, 5, 5, 0, 8, 3, 8]);
        assert_eq!(a.transpose(), expected);
    }

    #[test]
    fn transposing_identity_matrix() {
        let a: Matrix<u8> = Matrix::identity(4);
        assert_eq!(a.transpose(), a);
    }

    #[test]
    fn calculating_determinant_of_2x2_matrix() {
        let a = Matrix::from_vec(2, 2, vec![1, 5, -3, 2]);
        assert_eq!(a.determinant(), 17);
    }

    #[test]
    fn submatrix_of_3x3_is_2x2() {
        let a = Matrix::from_vec(3, 3, vec![1, 5, 0, -3, 2, 7, 0, 6, -3]);
        let expected = Matrix::from_vec(2, 2, vec![-3, 2, 0, 6]);
        assert_eq!(a.submatrix(0, 2), expected);
    }

    #[test]
    fn submatrix_of_4x4_is_3x3() {
        let a = Matrix::from_vec(
            4,
            4,
            vec![-6, 1, 1, 6, -8, 5, 8, 6, -1, 0, 8, 2, -7, 1, -1, 1],
        );
        let expected = Matrix::from_vec(3, 3, vec![-6, 1, 6, -8, 8, 6, -7, -1, 1]);
        assert_eq!(a.submatrix(2, 1), expected);
    }

    #[test]
    fn calculating_minor_of_3x3_matrix() {
        let a = Matrix::from_vec(3, 3, vec![1, 5, 0, 2, -1, -7, 6, -1, 5]);
        let b = a.submatrix(1, 0);
        assert_eq!(b.determinant(), 25);
        assert_eq!(a.minor(1, 0), 25);
    }

    #[test]
    fn calculating_cofactor_of_3x3_matrix() {
        let a = Matrix::from_vec(3, 3, vec![3, 5, 0, 2, -1, -7, 6, -1, 5]);
        assert_eq!(a.minor(0, 0), -12);
        assert_eq!(a.cofactor(0, 0), -12);
        assert_eq!(a.minor(1, 0), 25);
        assert_eq!(a.cofactor(1, 0), -25);
    }

    #[test]
    fn calculating_determinant_of_3x3_matrix() {
        let a = Matrix::from_vec(3, 3, vec![1, 2, 6, -5, 8, -4, 2, 6, 4]);
        assert_eq!(a.cofactor(0, 0), 56);
        assert_eq!(a.cofactor(0, 1), 12);
        assert_eq!(a.cofactor(0, 2), -46);
        assert_eq!(a.determinant(), -196);
    }

    #[test]
    fn calculating_determinant_of_4x4_matrix() {
        let a = Matrix::from_vec(
            4,
            4,
            vec![-2, -8, 3, 5, -3, 1, 7, 3, 1, 2, -9, 6, -6, 7, 7, -9],
        );
        assert_eq!(a.cofactor(0, 0), 690);
        assert_eq!(a.cofactor(0, 1), 447);
        assert_eq!(a.cofactor(0, 2), 210);
        assert_eq!(a.cofactor(0, 3), 51);
        assert_eq!(a.determinant(), -4071);
    }

    #[test]
    fn testing_invertible_matrix_for_invertibility() {
        let a = Matrix::from_vec(
            4,
            4,
            vec![6, 4, 4, 4, 5, 5, 7, 6, 4, -9, 3, -7, 9, 1, 7, -6],
        );
        assert_eq!(a.determinant(), -2120);
        assert!(a.is_invertible());
    }

    #[test]
    fn testing_non_invertible_matrix_for_invertibility() {
        let a = Matrix::from_vec(
            4,
            4,
            vec![-4, 2, -2, -3, 9, 6, 2, 6, 0, -5, 1, -5, 0, 0, 0, 0],
        );
        assert_eq!(a.determinant(), 0);
        assert!(!a.is_invertible());
    }

    #[test]
    fn calculating_inverse_of_matrix() {
        let a = Matrix::from_vec(
            4,
            4,
            vec![
                -5., 2., 6., -8., 1., -5., 1., 8., 7., 7., -6., -7., 1., -3., 7., 4.,
            ],
        );
        let b = a.inverse();
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
        assert_eq!(a.determinant(), 532.);
        assert_eq!(a.cofactor(2, 3), -160.);
        assert_eq!(b.get(3, 2), -160. / 532.);
        assert_eq!(a.cofactor(3, 2), 105.);
        assert_eq!(b.get(2, 3), 105. / 532.);
        assert_eq!(b, expected);
    }

    #[test]
    fn calculating_inverse_of_matrix_2() {
        let a = Matrix::from_vec(
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
        assert_eq!(a.inverse(), expected);
    }

    #[test]
    fn calculating_inverse_of_matrix_3() {
        let a = Matrix::from_vec(
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
        assert_eq!(a.inverse(), expected);
    }

    #[test]
    fn multiplying_product_by_inverse() {
        let a = Matrix::from_vec(
            4,
            4,
            vec![
                3., -9., 7., 3., 3., -8., 2., -9., -4., 4., 4., 1., -6., 5., -1., 1.,
            ],
        );
        let b = Matrix::from_vec(
            4,
            4,
            vec![
                8., 2., 2., 2., 3., -1., 7., 0., 7., 0., 5., 4., 6., -2., 0., 5.,
            ],
        );
        let c = a.clone() * b.clone();
        let result = c * b.inverse();
        let epsilon = 0.00000000000001;
        assert!(a.equals(&result, epsilon));
    }

    #[test]
    fn portrait_matrix() {
        let mut a = Matrix::new(5, 2, 0);
        a.set(0, 0, 0);
        a.set(0, 1, 1);
        a.set(1, 0, 2);
        a.set(1, 1, 3);
        a.set(2, 0, 4);
        a.set(2, 1, 5);
        a.set(3, 0, 6);
        a.set(3, 1, 7);
        a.set(4, 0, 8);
        a.set(4, 1, 9);
        let expected = Matrix::from_vec(5, 2, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
        assert_eq!(a, expected);
    }

    #[test]
    fn landscape_matrix() {
        let mut a = Matrix::new(2, 5, 0);
        a.set(0, 0, 0);
        a.set(0, 1, 1);
        a.set(0, 2, 2);
        a.set(0, 3, 3);
        a.set(0, 4, 4);
        a.set(1, 0, 5);
        a.set(1, 1, 6);
        a.set(1, 2, 7);
        a.set(1, 3, 8);
        a.set(1, 4, 9);
        let expected = Matrix::from_vec(2, 5, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
        assert_eq!(a, expected);
    }

    #[test]
    fn creating_canvas() {
        let c = Matrix::new(20, 10, Tuple::color(0.0, 0.0, 0.0));
        assert_eq!(c.width(), 10);
        assert_eq!(c.height(), 20);
    }

    #[test]
    fn writing_pixels_to_canvas() {
        let mut c = Matrix::new(20, 10, Tuple::color(0.0, 0.0, 0.0));
        let red = Tuple::color(1.0, 0.0, 0.0);
        c.set(2, 3, red.clone());
        assert_eq!(c.get(2, 3), red);
    }

    #[test]
    fn constructing_ppm_header() {
        let c = Matrix::new(3, 5, Tuple::color(0.0, 0.0, 0.0));
        let ppm = format!("{}", c);
        let lines: Vec<&str> = ppm.split('\n').collect();
        assert_eq!(lines[0], "P3");
        assert_eq!(lines[1], "5 3");
        assert_eq!(lines[2], "255");
    }

    #[test]
    fn constructing_ppm_pixel_data() {
        let mut c = Matrix::new(3, 5, Tuple::color(0.0, 0.0, 0.0));
        let c1 = Tuple::color(1.5, 0.0, 0.0);
        let c2 = Tuple::color(0.0, 0.5, 0.0);
        let c3 = Tuple::color(-0.5, 0.0, 1.0);
        c.set(0, 0, c1);
        c.set(1, 2, c2);
        c.set(2, 4, c3);
        let ppm = format!("{}", c);
        let lines: Vec<&str> = ppm.split('\n').collect();
        assert_eq!(lines[3], "255 0 0 0 0 0 0 0 0 0 0 0 0 0 0");
        assert_eq!(lines[4], "0 0 0 0 0 0 0 128 0 0 0 0 0 0 0");
        assert_eq!(lines[5], "0 0 0 0 0 0 0 0 0 0 0 0 0 0 255");
    }

    #[test]
    fn splitting_long_lines_in_ppm_files() {
        let c = Matrix::new(2, 10, Tuple::color(1.0, 0.8, 0.6));
        let ppm = format!("{}", c);
        let lines: Vec<&str> = ppm.split('\n').collect();
        assert_eq!(
            lines[3],
            "255 204 153 255 204 153 255 204 153 255 204 153 255 204 153 255 204"
        );
        assert_eq!(
            lines[4],
            "153 255 204 153 255 204 153 255 204 153 255 204 153"
        );
        assert_eq!(
            lines[5],
            "255 204 153 255 204 153 255 204 153 255 204 153 255 204 153 255 204"
        );
        assert_eq!(
            lines[6],
            "153 255 204 153 255 204 153 255 204 153 255 204 153"
        );
    }

    #[test]
    fn ppm_files_are_terminated_by_a_newline_character() {
        let c = Matrix::new(3, 5, Tuple::color(0.0, 0.0, 0.0));
        let ppm = format!("{}", c);
        let lines: Vec<&str> = ppm.split('\n').collect();
        assert_eq!(lines[lines.len() - 1], "");
    }
}
