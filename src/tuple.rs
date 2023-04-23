use num::Num;
use std::iter::Sum;
use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Clone, Debug, PartialEq)]
pub struct Tuple<T> {
    data: Vec<T>,
}

impl<T> Tuple<T>
where
    T: Clone,
{
    pub fn new(dimension: usize, init: T) -> Self {
        Self {
            data: vec![init; dimension],
        }
    }

    pub fn from_vec(data: Vec<T>) -> Self {
        Self { data }
    }

    pub fn dimension(&self) -> usize {
        self.data.len()
    }

    pub fn set(&mut self, index: usize, value: T) {
        assert!(index < self.dimension());
        self.data[index] = value;
    }
}

impl<T> Tuple<T>
where
    T: Copy + Clone,
{
    pub fn get(&self, index: usize) -> T {
        assert!(index < self.dimension());
        self.data[index]
    }

    pub fn r(&self) -> T {
        self.data[0]
    }

    pub fn g(&self) -> T {
        self.data[1]
    }

    pub fn b(&self) -> T {
        self.data[2]
    }

    pub fn a(&self) -> T {
        self.data[3]
    }

    pub fn x(&self) -> T {
        self.data[0]
    }

    pub fn y(&self) -> T {
        self.data[1]
    }

    pub fn z(&self) -> T {
        self.data[2]
    }

    pub fn w(&self) -> T {
        self.data[3]
    }
}

impl<T> Tuple<T>
where
    T: Num + Copy + Clone,
{
    pub fn dot(&self, other: Self) -> T
    where
        T: Sum,
    {
        self.data
            .iter()
            .zip(other.data.iter())
            .map(|(x, y)| *x * *y)
            .sum()
    }

    pub fn color(r: T, g: T, b: T) -> Self {
        Self::from_vec(vec![r, g, b])
    }

    pub fn point(x: T, y: T, z: T) -> Self {
        Self::from_vec(vec![x, y, z, T::one()])
    }

    pub fn vector(x: T, y: T, z: T) -> Self {
        Self::from_vec(vec![x, y, z, T::zero()])
    }

    pub fn is_point(&self) -> bool {
        self.data[3] == T::one()
    }

    pub fn is_vector(&self) -> bool {
        self.data[3] == T::zero()
    }

    pub fn cross(&self, other: &Self) -> Self {
        Self::vector(
            self.y() * other.z() - self.z() * other.y(),
            self.z() * other.x() - self.x() * other.z(),
            self.x() * other.y() - self.y() * other.x(),
        )
    }
}

impl Tuple<f64> {
    pub fn magnitude(&self) -> f64 {
        (self.data.iter().map(|x| x * x).sum::<f64>()).sqrt()
    }

    pub fn normalize(&self) -> Self {
        let magnitude = self.magnitude();
        Self::from_vec(self.data.iter().map(|x| x / magnitude).collect())
    }
}

impl<T> Add for Tuple<T>
where
    T: Num + Copy + Clone,
{
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self::from_vec(
            self.data
                .iter()
                .zip(other.data.iter())
                .map(|(x, y)| *x + *y)
                .collect(),
        )
    }
}

impl<T> Sub for Tuple<T>
where
    T: Num + Copy + Clone,
{
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self::from_vec(
            self.data
                .iter()
                .zip(other.data.iter())
                .map(|(x, y)| *x - *y)
                .collect(),
        )
    }
}

impl<T> Neg for Tuple<T>
where
    T: Num + Copy + Clone + Neg<Output = T>,
{
    type Output = Self;

    fn neg(self) -> Self {
        Self::from_vec(self.data.iter().map(|x| -*x).collect())
    }
}

impl<T> Mul<T> for Tuple<T>
where
    T: Num + Copy + Clone,
{
    type Output = Self;

    fn mul(self, scalar: T) -> Self {
        Self::from_vec(self.data.iter().map(|x| *x * scalar).collect())
    }
}

impl<T> Mul<Tuple<T>> for Tuple<T>
where
    T: Num + Copy + Clone,
{
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Self::from_vec(
            self.data
                .iter()
                .zip(other.data.iter())
                .map(|(x, y)| *x * *y)
                .collect(),
        )
    }
}

impl<T> Div<T> for Tuple<T>
where
    T: Num + Copy + Clone,
{
    type Output = Self;

    fn div(self, scalar: T) -> Self {
        Self::from_vec(self.data.iter().map(|x| *x / scalar).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tuple_with_w_1_is_point() {
        let a = Tuple::from_vec(vec![4.3, -4.2, 3.1, 1.0]);
        assert!(a.is_point());
    }

    #[test]
    fn tuple_with_w_0_is_vector() {
        let a = Tuple::from_vec(vec![4.3, -4.2, 3.1, 0.0]);
        assert!(a.is_vector());
    }

    #[test]
    fn point_creates_tuples_with_w_1() {
        let a = Tuple::point(4.0, -4.0, 3.0);
        assert!(a.is_point());
    }

    #[test]
    fn vector_creates_tuples_with_w_0() {
        let a = Tuple::vector(4.0, -4.0, 3.0);
        assert!(a.is_vector());
    }

    #[test]
    fn equality_with_identical_tuples() {
        let a = Tuple::from_vec(vec![1.0, 2.0, 3.0, 4.0]);
        let b = Tuple::from_vec(vec![1.0, 2.0, 3.0, 4.0]);
        assert_eq!(a, b);
    }

    #[test]
    fn equality_with_different_tuples() {
        let a = Tuple::from_vec(vec![1.0, 2.0, 3.0, 4.0]);
        let b = Tuple::from_vec(vec![2.0, 3.0, 4.0, 5.0]);
        assert_ne!(a, b);
    }

    #[test]
    fn adding_two_tuples() {
        let a1 = Tuple::from_vec(vec![3.0, -2.0, 5.0, 1.0]);
        let a2 = Tuple::from_vec(vec![-2.0, 3.0, 1.0, 0.0]);
        let expected = Tuple::from_vec(vec![1.0, 1.0, 6.0, 1.0]);
        assert_eq!(a1 + a2, expected);
    }

    #[test]
    fn subtracting_two_points() {
        let p1 = Tuple::point(3.0, 2.0, 1.0);
        let p2 = Tuple::point(5.0, 6.0, 7.0);
        let expected = Tuple::vector(-2.0, -4.0, -6.0);
        assert_eq!(p1 - p2, expected);
    }

    #[test]
    fn subtracting_vector_from_point() {
        let p = Tuple::point(3.0, 2.0, 1.0);
        let v = Tuple::vector(5.0, 6.0, 7.0);
        let expected = Tuple::point(-2.0, -4.0, -6.0);
        assert_eq!(p - v, expected);
    }

    #[test]
    fn subtracting_two_vectors() {
        let v1 = Tuple::vector(3.0, 2.0, 1.0);
        let v2 = Tuple::vector(5.0, 6.0, 7.0);
        let expected = Tuple::vector(-2.0, -4.0, -6.0);
        assert_eq!(v1 - v2, expected);
    }

    #[test]
    fn subtracting_vector_from_zero_vector() {
        let zero = Tuple::vector(0.0, 0.0, 0.0);
        let v = Tuple::vector(1.0, -2.0, 3.0);
        let expected = Tuple::vector(-1.0, 2.0, -3.0);
        assert_eq!(zero - v, expected);
    }

    #[test]
    fn negating_a_tuple() {
        let a = Tuple::from_vec(vec![1.0, -2.0, 3.0, -4.0]);
        let expected = Tuple::from_vec(vec![-1.0, 2.0, -3.0, 4.0]);
        assert_eq!(-a, expected);
    }

    #[test]
    fn multiplying_tuple_by_scalar() {
        let a = Tuple::from_vec(vec![1.0, -2.0, 3.0, -4.0]);
        let expected = Tuple::from_vec(vec![3.5, -7.0, 10.5, -14.0]);
        assert_eq!(a * 3.5, expected);
    }

    #[test]
    fn multiplying_tuple_by_fraction() {
        let a = Tuple::from_vec(vec![1.0, -2.0, 3.0, -4.0]);
        let expected = Tuple::from_vec(vec![0.5, -1.0, 1.5, -2.0]);
        assert_eq!(a * 0.5, expected);
    }

    #[test]
    fn dividing_tuple_by_scalar() {
        let a = Tuple::from_vec(vec![1.0, -2.0, 3.0, -4.0]);
        let expected = Tuple::from_vec(vec![0.5, -1.0, 1.5, -2.0]);
        assert_eq!(a / 2.0, expected);
    }

    #[test]
    fn computing_magnitude_of_vector_1_0_0() {
        let v = Tuple::vector(1.0, 0.0, 0.0);
        assert_eq!(v.magnitude(), 1.0);
    }

    #[test]
    fn computing_magnitude_of_vector_0_1_0() {
        let v = Tuple::vector(0.0, 1.0, 0.0);
        assert_eq!(v.magnitude(), 1.0);
    }

    #[test]
    fn computing_magnitude_of_vector_0_0_1() {
        let v = Tuple::vector(0.0, 0.0, 1.0);
        assert_eq!(v.magnitude(), 1.0);
    }

    #[test]
    fn computing_magnitude_of_vector_1_2_3() {
        let v = Tuple::vector(1.0, 2.0, 3.0);
        assert_eq!(v.magnitude(), 14.0_f64.sqrt());
    }

    #[test]
    fn computing_magnitude_of_vector_neg1_neg2_neg3() {
        let v = Tuple::vector(-1.0, -2.0, -3.0);
        assert_eq!(v.magnitude(), 14.0_f64.sqrt());
    }

    #[test]
    fn normalizing_vector_4_0_0_gives_1_0_0() {
        let v = Tuple::vector(4.0, 0.0, 0.0);
        let expected = Tuple::vector(1.0, 0.0, 0.0);
        assert_eq!(v.normalize(), expected);
    }

    #[test]
    fn normalizing_vector_1_2_3() {
        let v = Tuple::vector(1.0, 2.0, 3.0);
        let expected = Tuple::vector(
            1.0 / 14.0_f64.sqrt(),
            2.0 / 14.0_f64.sqrt(),
            3.0 / 14.0_f64.sqrt(),
        );
        assert_eq!(v.normalize(), expected);
    }

    #[test]
    fn magnitude_of_normalized_vector() {
        let v = Tuple::vector(1.0, 2.0, 3.0);
        let normalized = v.normalize();
        assert_eq!(normalized.magnitude(), 1.0);
    }

    #[test]
    fn dot_product_of_two_tuples() {
        let a = Tuple::vector(1.0, 2.0, 3.0);
        let b = Tuple::vector(2.0, 3.0, 4.0);
        assert_eq!(a.dot(b), 20.0);
    }

    #[test]
    fn cross_product_of_two_vectors() {
        let a = Tuple::vector(1.0, 2.0, 3.0);
        let b = Tuple::vector(2.0, 3.0, 4.0);
        let expected_a = Tuple::vector(-1.0, 2.0, -1.0);
        let expected_b = Tuple::vector(1.0, -2.0, 1.0);
        assert_eq!(a.cross(&b), expected_a);
        assert_eq!(b.cross(&a), expected_b);
    }

    #[test]
    fn creating_color() {
        let c = Tuple::color(-0.5, 0.4, 1.7);
        assert_eq!(c.r(), -0.5);
        assert_eq!(c.g(), 0.4);
        assert_eq!(c.b(), 1.7);
    }

    #[test]
    fn adding_colors() {
        let c1 = Tuple::color(0.9, 0.6, 0.75);
        let c2 = Tuple::color(0.7, 0.1, 0.25);
        assert_eq!(c1 + c2, Tuple::color(1.6, 0.7, 1.0));
    }

    #[test]
    fn subtracting_colors() {
        let c1 = Tuple::color(0.9, 0.6, 0.75);
        let c2 = Tuple::color(0.7, 0.1, 0.25);
        // check each component with a small epsilon
        let result = c1 - c2;
        assert!(result.r() - 0.2 < 0.00001);
        assert!(result.g() - 0.5 < 0.00001);
        assert!(result.b() - 0.5 < 0.00001);
    }

    #[test]
    fn multiplying_color_by_scalar() {
        let c = Tuple::color(0.2, 0.3, 0.4);
        assert_eq!(c * 2.0, Tuple::color(0.4, 0.6, 0.8));
    }

    #[test]
    fn multiplying_colors() {
        let c1 = Tuple::color(1.0, 0.2, 0.4);
        let c2 = Tuple::color(0.9, 1.0, 0.1);
        // check each component with a small epsilon
        let result = c1 * c2;
        assert!(result.r() - 0.9 < 0.00001);
        assert!(result.g() - 0.2 < 0.00001);
        assert!(result.b() - 0.04 < 0.00001);
    }
}
