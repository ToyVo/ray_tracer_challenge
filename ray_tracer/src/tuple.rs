use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::Matrix;

#[derive(Clone, Debug, PartialEq)]
pub struct Tuple {
    data: Vec<f64>,
}

impl Tuple {
    pub fn new(dimension: usize, init: f64) -> Tuple {
        Tuple {
            data: vec![init; dimension],
        }
    }

    pub fn from_vec(data: Vec<f64>) -> Tuple {
        Tuple { data }
    }

    pub fn color(r: f64, g: f64, b: f64) -> Tuple {
        Tuple::from_vec(vec![r, g, b])
    }

    pub fn point(x: f64, y: f64, z: f64) -> Tuple {
        Tuple::from_vec(vec![x, y, z, 1.])
    }

    pub fn vector(x: f64, y: f64, z: f64) -> Tuple {
        Tuple::from_vec(vec![x, y, z, 0.])
    }

    pub fn translate(&self, x: f64, y: f64, z: f64) -> Tuple {
        Matrix::translation(x, y, z) * self
    }

    pub fn scale(&self, x: f64, y: f64, z: f64) -> Tuple {
        Matrix::scaling(x, y, z) * self
    }

    pub fn shear(&self, xy: f64, xz: f64, yx: f64, yz: f64, zx: f64, zy: f64) -> Tuple {
        Matrix::shearing(xy, xz, yx, yz, zx, zy) * self
    }

    pub fn rotate_x(&self, rad: f64) -> Tuple {
        Matrix::rotation_x(rad) * self
    }

    pub fn rotate_y(&self, rad: f64) -> Tuple {
        Matrix::rotation_y(rad) * self
    }

    pub fn rotate_z(&self, rad: f64) -> Tuple {
        Matrix::rotation_z(rad) * self
    }

    pub fn to_vector(&self) -> Tuple {
        Tuple::from_vec(self.data.iter().enumerate().map(|(i, &v)| if i == self.data.len() - 1 { 0. } else { v }).collect())
    }

    pub fn dimension(&self) -> usize {
        self.data.len()
    }

    pub fn set(&mut self, index: usize, value: f64) {
        assert!(index < self.dimension());
        self.data[index] = value;
    }
    pub fn get(&self, index: usize) -> f64 {
        assert!(index < self.dimension());
        self.data[index]
    }

    pub fn r(&self) -> f64 {
        self.data[0]
    }

    pub fn g(&self) -> f64 {
        self.data[1]
    }

    pub fn b(&self) -> f64 {
        self.data[2]
    }

    pub fn a(&self) -> f64 {
        self.data[3]
    }

    pub fn x(&self) -> f64 {
        self.data[0]
    }

    pub fn y(&self) -> f64 {
        self.data[1]
    }

    pub fn z(&self) -> f64 {
        self.data[2]
    }

    pub fn w(&self) -> f64 {
        self.data[3]
    }
    pub fn dot(&self, other: &Tuple) -> f64 {
        self.data
            .iter()
            .zip(other.data.iter())
            .map(|(x, y)| *x * *y)
            .sum()
    }

    pub fn is_point(&self) -> bool {
        self.w() == 1.
    }

    pub fn is_vector(&self) -> bool {
        self.w() == 0.
    }

    pub fn cross(&self, other: &Tuple) -> Tuple {
        Tuple::vector(
            self.y() * other.z() - self.z() * other.y(),
            self.z() * other.x() - self.x() * other.z(),
            self.x() * other.y() - self.y() * other.x(),
        )
    }

    pub fn magnitude(&self) -> f64 {
        (self.data.iter().map(|x| x * x).sum::<f64>()).sqrt()
    }

    pub fn normalize(&self) -> Tuple {
        let magnitude = self.magnitude();
        Tuple::from_vec(self.data.iter().map(|x| x / magnitude).collect())
    }

    pub fn nearly_equals(&self, other: &Tuple, delta: f64) -> bool {
        for (x, y) in self.data.iter().zip(other.data.iter()) {
            if (x - y).abs() > delta {
                return false;
            }
        }
        true
    }

    pub fn reflect(&self, normal: &Tuple) -> Tuple {
        self - &(normal * 2. * self.dot(normal))
    }
}

impl Add<&Tuple> for &Tuple {
    type Output = Tuple;

    fn add(self, other: &Tuple) -> Tuple {
        Tuple::from_vec(
            self.data
                .iter()
                .zip(other.data.iter())
                .map(|(x, y)| x + y)
                .collect(),
        )
    }
}

impl Add<Tuple> for Tuple {
    type Output = Tuple;

    fn add(self, other: Tuple) -> Tuple {
        &self + &other
    }
}

impl Add<&Tuple> for Tuple {
    type Output = Tuple;

    fn add(self, other: &Tuple) -> Tuple {
        &self + other
    }
}

impl Add<Tuple> for &Tuple {
    type Output = Tuple;

    fn add(self, other: Tuple) -> Tuple {
        self + &other
    }
}

impl Sub<&Tuple> for &Tuple {
    type Output = Tuple;

    fn sub(self, other: &Tuple) -> Tuple {
        Tuple::from_vec(
            self.data
                .iter()
                .zip(other.data.iter())
                .map(|(x, y)| *x - *y)
                .collect(),
        )
    }
}

impl Sub<Tuple> for Tuple {
    type Output = Tuple;

    fn sub(self, other: Tuple) -> Tuple {
        &self - &other
    }
}

impl Sub<&Tuple> for Tuple {
    type Output = Tuple;

    fn sub(self, other: &Tuple) -> Tuple {
        &self - other
    }
}

impl Sub<Tuple> for &Tuple {
    type Output = Tuple;

    fn sub(self, other: Tuple) -> Tuple {
        self - &other
    }
}

impl Neg for &Tuple {
    type Output = Tuple;

    fn neg(self) -> Tuple {
        Tuple::from_vec(self.data.iter().map(|x| -*x).collect())
    }
}

impl Neg for Tuple {
    type Output = Tuple;

    fn neg(self) -> Tuple {
        -&self
    }
}

impl Mul<f64> for &Tuple {
    type Output = Tuple;

    fn mul(self, scalar: f64) -> Tuple {
        Tuple::from_vec(self.data.iter().map(|x| *x * scalar).collect())
    }
}

impl Mul<f64> for Tuple {
    type Output = Tuple;

    fn mul(self, scalar: f64) -> Tuple {
        &self * scalar
    }
}

impl Mul<&Tuple> for &Tuple {
    type Output = Tuple;

    fn mul(self, other: &Tuple) -> Tuple {
        Tuple::from_vec(
            self.data
                .iter()
                .zip(other.data.iter())
                .map(|(x, y)| *x * *y)
                .collect(),
        )
    }
}

impl Mul<Tuple> for Tuple {
    type Output = Tuple;

    fn mul(self, other: Tuple) -> Tuple {
        &self * &other
    }
}

impl Mul<&Tuple> for Tuple {
    type Output = Tuple;

    fn mul(self, other: &Tuple) -> Tuple {
        &self * other
    }
}

impl Mul<Tuple> for &Tuple {
    type Output = Tuple;

    fn mul(self, other: Tuple) -> Tuple {
        self * &other
    }
}

impl Div<f64> for &Tuple {
    type Output = Tuple;

    fn div(self, scalar: f64) -> Tuple {
        Tuple::from_vec(self.data.iter().map(|x| *x / scalar).collect())
    }
}

impl Div<f64> for Tuple {
    type Output = Tuple;

    fn div(self, scalar: f64) -> Tuple {
        &self / scalar
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tuple_with_w_1_is_point() {
        let point = Tuple::from_vec(vec![4.3, -4.2, 3.1, 1.0]);
        assert!(point.is_point());
    }

    #[test]
    fn tuple_with_w_0_is_vector() {
        let vector = Tuple::from_vec(vec![4.3, -4.2, 3.1, 0.0]);
        assert!(vector.is_vector());
    }

    #[test]
    fn point_creates_tuples_with_w_1() {
        let point = Tuple::point(4.0, -4.0, 3.0);
        assert!(point.is_point());
    }

    #[test]
    fn vector_creates_tuples_with_w_0() {
        let vector = Tuple::vector(4.0, -4.0, 3.0);
        assert!(vector.is_vector());
    }

    #[test]
    fn equality_with_identical_tuples() {
        let tuple_a = Tuple::from_vec(vec![1.0, 2.0, 3.0, 4.0]);
        let tuple_b = Tuple::from_vec(vec![1.0, 2.0, 3.0, 4.0]);
        assert_eq!(tuple_a, tuple_b);
    }

    #[test]
    fn equality_with_different_tuples() {
        let tuple_a = Tuple::from_vec(vec![1.0, 2.0, 3.0, 4.0]);
        let tuple_b = Tuple::from_vec(vec![2.0, 3.0, 4.0, 5.0]);
        assert_ne!(tuple_a, tuple_b);
    }

    #[test]
    fn adding_two_tuples() {
        let tuple_a = Tuple::from_vec(vec![3.0, -2.0, 5.0, 1.0]);
        let tuple_b = Tuple::from_vec(vec![-2.0, 3.0, 1.0, 0.0]);
        let expected = Tuple::from_vec(vec![1.0, 1.0, 6.0, 1.0]);
        assert_eq!(tuple_a + tuple_b, expected);
    }

    #[test]
    fn subtracting_two_points() {
        let point_a = Tuple::point(3.0, 2.0, 1.0);
        let point_b = Tuple::point(5.0, 6.0, 7.0);
        let expected = Tuple::vector(-2.0, -4.0, -6.0);
        assert_eq!(point_a - point_b, expected);
    }

    #[test]
    fn subtracting_vector_from_point() {
        let point = Tuple::point(3.0, 2.0, 1.0);
        let vector = Tuple::vector(5.0, 6.0, 7.0);
        let expected = Tuple::point(-2.0, -4.0, -6.0);
        assert_eq!(point - vector, expected);
    }

    #[test]
    fn subtracting_two_vectors() {
        let vector_a = Tuple::vector(3.0, 2.0, 1.0);
        let vector_b = Tuple::vector(5.0, 6.0, 7.0);
        let expected = Tuple::vector(-2.0, -4.0, -6.0);
        assert_eq!(vector_a - vector_b, expected);
    }

    #[test]
    fn subtracting_vector_from_zero_vector() {
        let zero = Tuple::vector(0.0, 0.0, 0.0);
        let vector = Tuple::vector(1.0, -2.0, 3.0);
        let expected = Tuple::vector(-1.0, 2.0, -3.0);
        assert_eq!(zero - vector, expected);
    }

    #[test]
    fn negating_a_tuple() {
        let vector = Tuple::from_vec(vec![1.0, -2.0, 3.0, -4.0]);
        let expected = Tuple::from_vec(vec![-1.0, 2.0, -3.0, 4.0]);
        assert_eq!(-vector, expected);
    }

    #[test]
    fn multiplying_tuple_by_scalar() {
        let vector = Tuple::from_vec(vec![1.0, -2.0, 3.0, -4.0]);
        let expected = Tuple::from_vec(vec![3.5, -7.0, 10.5, -14.0]);
        assert_eq!(vector * 3.5, expected);
    }

    #[test]
    fn multiplying_tuple_by_fraction() {
        let vector = Tuple::from_vec(vec![1.0, -2.0, 3.0, -4.0]);
        let expected = Tuple::from_vec(vec![0.5, -1.0, 1.5, -2.0]);
        assert_eq!(vector * 0.5, expected);
    }

    #[test]
    fn dividing_tuple_by_scalar() {
        let vector = Tuple::from_vec(vec![1.0, -2.0, 3.0, -4.0]);
        let expected = Tuple::from_vec(vec![0.5, -1.0, 1.5, -2.0]);
        assert_eq!(vector / 2.0, expected);
    }

    #[test]
    fn computing_magnitude_of_vector_1_0_0() {
        let vector = Tuple::vector(1.0, 0.0, 0.0);
        assert_eq!(vector.magnitude(), 1.0);
    }

    #[test]
    fn computing_magnitude_of_vector_0_1_0() {
        let vector = Tuple::vector(0.0, 1.0, 0.0);
        assert_eq!(vector.magnitude(), 1.0);
    }

    #[test]
    fn computing_magnitude_of_vector_0_0_1() {
        let vector = Tuple::vector(0.0, 0.0, 1.0);
        assert_eq!(vector.magnitude(), 1.0);
    }

    #[test]
    fn computing_magnitude_of_vector_1_2_3() {
        let vector = Tuple::vector(1.0, 2.0, 3.0);
        assert_eq!(vector.magnitude(), 14.0_f64.sqrt());
    }

    #[test]
    fn computing_magnitude_of_vector_neg1_neg2_neg3() {
        let vector = Tuple::vector(-1.0, -2.0, -3.0);
        assert_eq!(vector.magnitude(), 14.0_f64.sqrt());
    }

    #[test]
    fn normalizing_vector_4_0_0_gives_1_0_0() {
        let vector = Tuple::vector(4.0, 0.0, 0.0);
        let expected = Tuple::vector(1.0, 0.0, 0.0);
        assert_eq!(vector.normalize(), expected);
    }

    #[test]
    fn normalizing_vector_1_2_3() {
        let vector = Tuple::vector(1.0, 2.0, 3.0);
        let expected = Tuple::vector(
            1.0 / 14.0_f64.sqrt(),
            2.0 / 14.0_f64.sqrt(),
            3.0 / 14.0_f64.sqrt(),
        );
        assert_eq!(vector.normalize(), expected);
    }

    #[test]
    fn magnitude_of_normalized_vector() {
        let vector = Tuple::vector(1.0, 2.0, 3.0);
        let normalized = vector.normalize();
        assert_eq!(normalized.magnitude(), 1.0);
    }

    #[test]
    fn dot_product_of_two_tuples() {
        let vector_a = Tuple::vector(1.0, 2.0, 3.0);
        let vector_b = Tuple::vector(2.0, 3.0, 4.0);
        assert_eq!(vector_a.dot(&vector_b), 20.0);
    }

    #[test]
    fn cross_product_of_two_vectors() {
        let vector_a = Tuple::vector(1.0, 2.0, 3.0);
        let vector_b = Tuple::vector(2.0, 3.0, 4.0);
        let expected_a = Tuple::vector(-1.0, 2.0, -1.0);
        let expected_b = Tuple::vector(1.0, -2.0, 1.0);
        assert_eq!(vector_a.cross(&vector_b), expected_a);
        assert_eq!(vector_b.cross(&vector_a), expected_b);
    }

    #[test]
    fn creating_color() {
        let color = Tuple::color(-0.5, 0.4, 1.7);
        assert_eq!(color.r(), -0.5);
        assert_eq!(color.g(), 0.4);
        assert_eq!(color.b(), 1.7);
    }

    #[test]
    fn adding_colors() {
        let color_a = Tuple::color(0.9, 0.6, 0.75);
        let color_b = Tuple::color(0.7, 0.1, 0.25);
        assert_eq!(color_a + color_b, Tuple::color(1.6, 0.7, 1.0));
    }

    #[test]
    fn subtracting_colors() {
        let color_a = Tuple::color(0.9, 0.6, 0.75);
        let color_b = Tuple::color(0.7, 0.1, 0.25);
        let expected = Tuple::color(0.2, 0.5, 0.5);
        let result = color_a - color_b;
        assert!(result.nearly_equals(&expected, 1e-3f64))
    }

    #[test]
    fn multiplying_color_by_scalar() {
        let color = Tuple::color(0.2, 0.3, 0.4);
        assert_eq!(color * 2.0, Tuple::color(0.4, 0.6, 0.8));
    }

    #[test]
    fn multiplying_colors() {
        let color_a = Tuple::color(1.0, 0.2, 0.4);
        let color_b = Tuple::color(0.9, 1.0, 0.1);
        let expected = Tuple::color(0.9, 0.2, 0.04);
        let result = color_a * color_b;
        assert!(result.nearly_equals(&expected, 1e-3f64))
    }

    #[test]
    fn to_vector_makes_the_last_element_0() {
        let point = Tuple::point(1., 2., 3.);
        let expected = Tuple::vector(1., 2., 3.);
        assert_eq!(point.to_vector(), expected);
    }

    #[test]
    fn reflect_vector_approaching_at_45_degrees() {
        let vector = Tuple::vector(1., -1., 0.);
        let normal = Tuple::vector(0., 1., 0.);
        let ray = vector.reflect(&normal);
        let expected = Tuple::vector(1., 1., 0.);
        assert_eq!(ray, expected);
    }

    #[test]
    fn reflect_vector_off_slanted_surface() {
        let vector = Tuple::vector(0., -1., 0.);
        let normal = Tuple::vector(2.0_f64.sqrt() / 2., 2.0_f64.sqrt() / 2., 0.);
        let ray = vector.reflect(&normal);
        let expected = Tuple::vector(1., 0., 0.);
        assert!(ray.nearly_equals(&expected, 1e-3f64));
    }
}
