use std::ops::{Add, Mul, Sub, Div, Neg};

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Tuple {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub w: f64,
}

impl Tuple {
    pub fn new(x: f64, y: f64, z: f64, w: f64) -> Tuple {
        Tuple { x, y, z, w }
    }

    pub fn point(x: f64, y: f64, z: f64) -> Tuple {
        Tuple::new(x, y, z, 1.0)
    }

    pub fn vector(x: f64, y: f64, z: f64) -> Tuple {
        Tuple::new(x, y, z, 0.0)
    }

    pub fn is_point(&self) -> bool {
        self.w == 1.0
    }

    pub fn is_vector(&self) -> bool {
        self.w == 0.0
    }

    pub fn magnitude(&self) -> f64 {
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2) + self.w.powi(2)).sqrt()
    }

    pub fn normalize(&self) -> Tuple {
        let magnitude = self.magnitude();
        Tuple::new(self.x / magnitude, self.y / magnitude, self.z / magnitude, self.w / magnitude)
    }

    pub fn dot(&self, other: Tuple) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z + self.w * other.w
    }

    pub fn cross(&self, other: &Tuple) -> Tuple {
        Tuple::vector(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }
}

impl Add for Tuple {
    type Output = Tuple;

    fn add(self, other: Tuple) -> Tuple {
        Tuple::new(
            self.x + other.x,
            self.y + other.y,
            self.z + other.z,
            self.w + other.w,
        )
    }
}

impl Sub for Tuple {
    type Output = Tuple;

    fn sub(self, other: Tuple) -> Tuple {
        Tuple::new(
            self.x - other.x,
            self.y - other.y,
            self.z - other.z,
            self.w - other.w,
        )
    }
}

impl Neg for Tuple {
    type Output = Tuple;

    fn neg(self) -> Tuple {
        Tuple::new(-self.x, -self.y, -self.z, -self.w)
    }
}

impl Mul<f64> for Tuple {
    type Output = Tuple;

    fn mul(self, scalar: f64) -> Tuple {
        Tuple::new(self.x * scalar, self.y * scalar, self.z * scalar, self.w * scalar)
    }
}

impl Div<f64> for Tuple {
    type Output = Tuple;

    fn div(self, scalar: f64) -> Tuple {
        Tuple::new(self.x / scalar, self.y / scalar, self.z / scalar, self.w / scalar)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tuple_with_w_1_is_point() {
        let a = Tuple::new(4.3, -4.2, 3.1, 1.0);
        assert!(a.is_point());
    }

    #[test]
    fn tuple_with_w_0_is_vector() {
        let a = Tuple::new(4.3, -4.2, 3.1, 0.0);
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
        let a = Tuple::new(1.0, 2.0, 3.0, 4.0);
        let b = Tuple::new(1.0, 2.0, 3.0, 4.0);
        assert!(a == b);
    }

    #[test]
    fn equality_with_different_tuples() {
        let a = Tuple::new(1.0, 2.0, 3.0, 4.0);
        let b = Tuple::new(2.0, 3.0, 4.0, 5.0);
        assert!(a != b);
    }

    #[test]
    fn adding_two_tuples() {
        let a1 = Tuple::new(3.0, -2.0, 5.0, 1.0);
        let a2 = Tuple::new(-2.0, 3.0, 1.0, 0.0);
        let expected = Tuple::new(1.0, 1.0, 6.0, 1.0);
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
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);
        let expected = Tuple::new(-1.0, 2.0, -3.0, 4.0);
        assert_eq!(-a, expected);
    }

    #[test]
    fn multiplying_tuple_by_scalar() {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);
        let expected = Tuple::new(3.5, -7.0, 10.5, -14.0);
        assert_eq!(a * 3.5, expected);
    }

    #[test]
    fn multiplying_tuple_by_fraction() {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);
        let expected = Tuple::new(0.5, -1.0, 1.5, -2.0);
        assert_eq!(a * 0.5, expected);
    }

    #[test]
    fn dividing_tuple_by_scalar() {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);
        let expected = Tuple::new(0.5, -1.0, 1.5, -2.0);
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
        assert_eq!(v.magnitude(), (14.0 as f64).sqrt());
    }

    #[test]
    fn computing_magnitude_of_vector_neg1_neg2_neg3() {
        let v = Tuple::vector(-1.0, -2.0, -3.0);
        assert_eq!(v.magnitude(), (14.0 as f64).sqrt());
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
        let expected = Tuple::vector(1.0 / (14.0 as f64).sqrt(), 2.0 / (14.0 as f64).sqrt(), 3.0 / (14.0 as f64).sqrt());
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

}
