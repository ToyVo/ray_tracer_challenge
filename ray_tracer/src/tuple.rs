use crate::Matrix;
use approx::{AbsDiffEq, RelativeEq, UlpsEq};
use auto_ops::*;
use std::ops::Neg;

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
        Tuple::from_vec(
            self.data
                .iter()
                .enumerate()
                .map(|(i, &v)| if i == self.data.len() - 1 { 0. } else { v })
                .collect(),
        )
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

    // function dot(a, b)
    //   return a.x * b.x +
    //          a.y * b.y +
    //          a.z * b.z +
    //          a.w * b.w
    // end function
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

    // function cross(a, b)
    //   return vector(a.y * b.z - a.z * b.y,
    //                 a.z * b.x - a.x * b.z,
    //                 a.x * b.y - a.y * b.x)
    // end function
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

    // function normalize(v)
    //   return tuple(v.x / magnitude(v),
    //                v.y / magnitude(v),
    //                v.z / magnitude(v),
    //                v.w / magnitude(v))
    // end function
    pub fn normalize(&self) -> Tuple {
        let magnitude = self.magnitude();
        Tuple::from_vec(self.data.iter().map(|x| x / magnitude).collect())
    }

    // function reflect(in, normal)
    //   return in - normal * 2 * dot(in, normal)
    // end function
    pub fn reflect(&self, normal: &Tuple) -> Tuple {
        self - &(normal * 2. * self.dot(normal))
    }
}

impl_op_ex!(+ |a: &Tuple, b: &Tuple| -> Tuple {
    Tuple::from_vec(
        a.data
            .iter()
            .zip(b.data.iter())
            .map(|(x, y)| x + y)
            .collect(),
    )
});

impl_op_ex!(-|a: &Tuple, b: &Tuple| -> Tuple {
    Tuple::from_vec(
        a.data
            .iter()
            .zip(b.data.iter())
            .map(|(x, y)| *x - *y)
            .collect(),
    )
});

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

impl_op_ex!(*|a: &Tuple, b: f64| -> Tuple {
    Tuple::from_vec(a.data.iter().map(|x| *x * b).collect())
});

// function hadamard_product(c1, c2)
//   r ← c1.red * c2.red
//   g ← c1.green * c2.green
//   b ← c1.blue * c2.blue
//   return color(r, g, b)
// end function
impl_op_ex!(*|a: &Tuple, b: &Tuple| -> Tuple {
    Tuple::from_vec(
        a.data
            .iter()
            .zip(b.data.iter())
            .map(|(x, y)| *x * *y)
            .collect(),
    )
});

impl_op_ex!(/ |a: &Tuple, b: f64| -> Tuple {
    Tuple::from_vec(a.data.iter().map(|x| *x / b).collect())
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
impl AbsDiffEq for Tuple {
    type Epsilon = f64;

    fn default_epsilon() -> f64 {
        1e-5f64
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        for (x, y) in self.data.iter().zip(other.data.iter()) {
            if !x.abs_diff_eq(y, epsilon) {
                return false;
            }
        }
        true
    }
}

impl RelativeEq for Tuple {
    fn default_max_relative() -> f64 {
        f64::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        for (x, y) in self.data.iter().zip(other.data.iter()) {
            if !x.relative_eq(y, epsilon, max_relative) {
                return false;
            }
        }
        true
    }
}

impl UlpsEq for Tuple {
    fn default_max_ulps() -> u32 {
        f64::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        for (x, y) in self.data.iter().zip(other.data.iter()) {
            if !x.ulps_eq(y, epsilon, max_ulps) {
                return false;
            }
        }
        true
    }
}

// Feature: Tuples, Vectors, and Points
#[cfg(test)]
mod tests {
    use std::f64::consts::SQRT_2;

    use super::*;
    use approx::assert_relative_eq;

    // Scenario: A tuple with w=1.0 is a point
    //   Given a ← tuple(4.3, -4.2, 3.1, 1.0)
    //   Then a.x = 4.3
    //     And a.y = -4.2
    //     And a.z = 3.1
    //     And a.w = 1.0
    //     And a is a point
    //     And a is not a vector
    #[test]
    fn tuple_with_w_1_is_point() {
        let point = Tuple::from_vec(vec![4.3, -4.2, 3.1, 1.0]);
        assert!(point.is_point());
    }

    // Scenario: A tuple with w=0 is a vector
    //   Given a ← tuple(4.3, -4.2, 3.1, 0.0)
    //   Then a.x = 4.3
    //     And a.y = -4.2
    //     And a.z = 3.1
    //     And a.w = 0.0
    //     And a is not a point
    //     And a is a vector
    #[test]
    fn tuple_with_w_0_is_vector() {
        let vector = Tuple::from_vec(vec![4.3, -4.2, 3.1, 0.0]);
        assert!(vector.is_vector());
    }

    // Scenario: point() creates tuples with w=1
    //   Given p ← point(4, -4, 3)
    //   Then p = tuple(4, -4, 3, 1)
    #[test]
    fn point_creates_tuples_with_w_1() {
        let point = Tuple::point(4.0, -4.0, 3.0);
        assert!(point.is_point());
    }

    // Scenario: vector() creates tuples with w=0
    //   Given v ← vector(4, -4, 3)
    //   Then v = tuple(4, -4, 3, 0)
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

    // Scenario: Adding two tuples
    //   Given a1 ← tuple(3, -2, 5, 1)
    //     And a2 ← tuple(-2, 3, 1, 0)
    //    Then a1 + a2 = tuple(1, 1, 6, 1)
    #[test]
    fn adding_two_tuples() {
        let tuple_a = Tuple::from_vec(vec![3.0, -2.0, 5.0, 1.0]);
        let tuple_b = Tuple::from_vec(vec![-2.0, 3.0, 1.0, 0.0]);
        let expected = Tuple::from_vec(vec![1.0, 1.0, 6.0, 1.0]);
        assert_eq!(tuple_a + tuple_b, expected);
    }

    // Scenario: Subtracting two points
    //   Given p1 ← point(3, 2, 1)
    //     And p2 ← point(5, 6, 7)
    //   Then p1 - p2 = vector(-2, -4, -6)
    #[test]
    fn subtracting_two_points() {
        let point_a = Tuple::point(3.0, 2.0, 1.0);
        let point_b = Tuple::point(5.0, 6.0, 7.0);
        let expected = Tuple::vector(-2.0, -4.0, -6.0);
        assert_eq!(point_a - point_b, expected);
    }

    // Scenario: Subtracting a vector from a point
    //   Given p ← point(3, 2, 1)
    //     And v ← vector(5, 6, 7)
    //   Then p - v = point(-2, -4, -6)
    #[test]
    fn subtracting_vector_from_point() {
        let point = Tuple::point(3.0, 2.0, 1.0);
        let vector = Tuple::vector(5.0, 6.0, 7.0);
        let expected = Tuple::point(-2.0, -4.0, -6.0);
        assert_eq!(point - vector, expected);
    }

    // Scenario: Subtracting two vectors
    //   Given v1 ← vector(3, 2, 1)
    //     And v2 ← vector(5, 6, 7)
    //   Then v1 - v2 = vector(-2, -4, -6)
    #[test]
    fn subtracting_two_vectors() {
        let vector_a = Tuple::vector(3.0, 2.0, 1.0);
        let vector_b = Tuple::vector(5.0, 6.0, 7.0);
        let expected = Tuple::vector(-2.0, -4.0, -6.0);
        assert_eq!(vector_a - vector_b, expected);
    }

    // Scenario: Subtracting a vector from the zero vector
    //   Given zero ← vector(0, 0, 0)
    //     And v ← vector(1, -2, 3)
    //   Then zero - v = vector(-1, 2, -3)
    #[test]
    fn subtracting_vector_from_zero_vector() {
        let zero = Tuple::vector(0.0, 0.0, 0.0);
        let vector = Tuple::vector(1.0, -2.0, 3.0);
        let expected = Tuple::vector(-1.0, 2.0, -3.0);
        assert_eq!(zero - vector, expected);
    }

    // Scenario: Negating a tuple
    //   Given a ← tuple(1, -2, 3, -4)
    //   Then -a = tuple(-1, 2, -3, 4)
    #[test]
    fn negating_a_tuple() {
        let vector = Tuple::from_vec(vec![1.0, -2.0, 3.0, -4.0]);
        let expected = Tuple::from_vec(vec![-1.0, 2.0, -3.0, 4.0]);
        assert_eq!(-vector, expected);
    }

    // Scenario: Multiplying a tuple by a scalar
    //   Given a ← tuple(1, -2, 3, -4)
    //   Then a * 3.5 = tuple(3.5, -7, 10.5, -14)
    #[test]
    fn multiplying_tuple_by_scalar() {
        let vector = Tuple::from_vec(vec![1.0, -2.0, 3.0, -4.0]);
        let expected = Tuple::from_vec(vec![3.5, -7.0, 10.5, -14.0]);
        assert_eq!(vector * 3.5, expected);
    }

    // Scenario: Multiplying a tuple by a fraction
    //   Given a ← tuple(1, -2, 3, -4)
    //   Then a * 0.5 = tuple(0.5, -1, 1.5, -2)
    #[test]
    fn multiplying_tuple_by_fraction() {
        let vector = Tuple::from_vec(vec![1.0, -2.0, 3.0, -4.0]);
        let expected = Tuple::from_vec(vec![0.5, -1.0, 1.5, -2.0]);
        assert_eq!(vector * 0.5, expected);
    }

    // Scenario: Dividing a tuple by a scalar
    //   Given a ← tuple(1, -2, 3, -4)
    //   Then a / 2 = tuple(0.5, -1, 1.5, -2)
    #[test]
    fn dividing_tuple_by_scalar() {
        let vector = Tuple::from_vec(vec![1.0, -2.0, 3.0, -4.0]);
        let expected = Tuple::from_vec(vec![0.5, -1.0, 1.5, -2.0]);
        assert_eq!(vector / 2.0, expected);
    }

    // Scenario: Computing the magnitude of vector(1, 0, 0)
    //   Given v ← vector(1, 0, 0)
    //   Then magnitude(v) = 1
    #[test]
    fn computing_magnitude_of_vector_1_0_0() {
        let vector = Tuple::vector(1.0, 0.0, 0.0);
        assert_eq!(vector.magnitude(), 1.0);
    }

    // Scenario: Computing the magnitude of vector(0, 1, 0)
    //   Given v ← vector(0, 1, 0)
    //   Then magnitude(v) = 1
    #[test]
    fn computing_magnitude_of_vector_0_1_0() {
        let vector = Tuple::vector(0.0, 1.0, 0.0);
        assert_eq!(vector.magnitude(), 1.0);
    }

    // Scenario: Computing the magnitude of vector(0, 0, 1)
    //   Given v ← vector(0, 0, 1)
    //   Then magnitude(v) = 1
    #[test]
    fn computing_magnitude_of_vector_0_0_1() {
        let vector = Tuple::vector(0.0, 0.0, 1.0);
        assert_eq!(vector.magnitude(), 1.0);
    }

    // Scenario: Computing the magnitude of vector(1, 2, 3)
    //   Given v ← vector(1, 2, 3)
    //   Then magnitude(v) = √14
    #[test]
    fn computing_magnitude_of_vector_1_2_3() {
        let vector = Tuple::vector(1.0, 2.0, 3.0);
        assert_eq!(vector.magnitude(), 14.0_f64.sqrt());
    }

    // Scenario: Computing the magnitude of vector(-1, -2, -3)
    //   Given v ← vector(-1, -2, -3)
    //   Then magnitude(v) = √14
    #[test]
    fn computing_magnitude_of_vector_neg1_neg2_neg3() {
        let vector = Tuple::vector(-1.0, -2.0, -3.0);
        assert_eq!(vector.magnitude(), 14.0_f64.sqrt());
    }

    // Scenario: Normalizing vector(4, 0, 0) gives (1, 0, 0)
    //   Given v ← vector(4, 0, 0)
    //   Then normalize(v) = vector(1, 0, 0)
    #[test]
    fn normalizing_vector_4_0_0_gives_1_0_0() {
        let vector = Tuple::vector(4.0, 0.0, 0.0);
        let expected = Tuple::vector(1.0, 0.0, 0.0);
        assert_eq!(vector.normalize(), expected);
    }

    // Scenario: Normalizing vector(1, 2, 3)
    //   Given v ← vector(1, 2, 3)
    //                                   # vector(1/√14,   2/√14,   3/√14)
    //   Then normalize(v) = approximately vector(0.26726, 0.53452, 0.80178)
    #[test]
    fn normalizing_vector_1_2_3() {
        let vector = Tuple::vector(1.0, 2.0, 3.0);
        let sqrt_14 = 14.0_f64.sqrt();
        let expected = Tuple::vector(1.0 / sqrt_14, 2.0 / sqrt_14, 3.0 / sqrt_14);
        assert_eq!(vector.normalize(), expected);
    }

    // Scenario: The magnitude of a normalized vector
    //   Given v ← vector(1, 2, 3)
    //   When norm ← normalize(v)
    //   Then magnitude(norm) = 1
    #[test]
    fn magnitude_of_normalized_vector() {
        let vector = Tuple::vector(1.0, 2.0, 3.0);
        let normalized = vector.normalize();
        assert_eq!(normalized.magnitude(), 1.0);
    }

    // Scenario: The dot product of two tuples
    //   Given a ← vector(1, 2, 3)
    //     And b ← vector(2, 3, 4)
    //   Then dot(a, b) = 20
    #[test]
    fn dot_product_of_two_tuples() {
        let vector_a = Tuple::vector(1.0, 2.0, 3.0);
        let vector_b = Tuple::vector(2.0, 3.0, 4.0);
        assert_eq!(vector_a.dot(&vector_b), 20.0);
    }

    // Scenario: The cross product of two vectors
    //   Given a ← vector(1, 2, 3)
    //     And b ← vector(2, 3, 4)
    //   Then cross(a, b) = vector(-1, 2, -1)
    //     And cross(b, a) = vector(1, -2, 1)
    #[test]
    fn cross_product_of_two_vectors() {
        let vector_a = Tuple::vector(1.0, 2.0, 3.0);
        let vector_b = Tuple::vector(2.0, 3.0, 4.0);
        let expected_a = Tuple::vector(-1.0, 2.0, -1.0);
        let expected_b = Tuple::vector(1.0, -2.0, 1.0);
        assert_eq!(vector_a.cross(&vector_b), expected_a);
        assert_eq!(vector_b.cross(&vector_a), expected_b);
    }

    // Scenario: Colors are (red, green, blue) tuples
    //   Given c ← color(-0.5, 0.4, 1.7)
    //   Then c.red = -0.5
    //     And c.green = 0.4
    //     And c.blue = 1.7
    #[test]
    fn creating_color() {
        let color = Tuple::color(-0.5, 0.4, 1.7);
        assert_eq!(color.r(), -0.5);
        assert_eq!(color.g(), 0.4);
        assert_eq!(color.b(), 1.7);
    }

    // Scenario: Adding colors
    //   Given c1 ← color(0.9, 0.6, 0.75)
    //     And c2 ← color(0.7, 0.1, 0.25)
    //    Then c1 + c2 = color(1.6, 0.7, 1.0)
    #[test]
    fn adding_colors() {
        let color_a = Tuple::color(0.9, 0.6, 0.75);
        let color_b = Tuple::color(0.7, 0.1, 0.25);
        assert_eq!(color_a + color_b, Tuple::color(1.6, 0.7, 1.0));
    }

    // Scenario: Subtracting colors
    //   Given c1 ← color(0.9, 0.6, 0.75)
    //     And c2 ← color(0.7, 0.1, 0.25)
    //    Then c1 - c2 = color(0.2, 0.5, 0.5)
    #[test]
    fn subtracting_colors() {
        let color_a = Tuple::color(0.9, 0.6, 0.75);
        let color_b = Tuple::color(0.7, 0.1, 0.25);
        let expected = Tuple::color(0.2, 0.5, 0.5);
        let result = color_a - color_b;
        assert_relative_eq!(result, expected);
    }

    // Scenario: Multiplying a color by a scalar
    //   Given c ← color(0.2, 0.3, 0.4)
    //   Then c * 2 = color(0.4, 0.6, 0.8)
    #[test]
    fn multiplying_color_by_scalar() {
        let color = Tuple::color(0.2, 0.3, 0.4);
        assert_eq!(color * 2.0, Tuple::color(0.4, 0.6, 0.8));
    }

    // Scenario: Multiplying colors
    //   Given c1 ← color(1, 0.2, 0.4)
    //     And c2 ← color(0.9, 1, 0.1)
    //    Then c1 * c2 = color(0.9, 0.2, 0.04)
    #[test]
    fn multiplying_colors() {
        let color_a = Tuple::color(1.0, 0.2, 0.4);
        let color_b = Tuple::color(0.9, 1.0, 0.1);
        let expected = Tuple::color(0.9, 0.2, 0.04);
        let result = color_a * color_b;
        assert_relative_eq!(result, expected);
    }

    #[test]
    fn to_vector_makes_the_last_element_0() {
        let point = Tuple::point(1., 2., 3.);
        let expected = Tuple::vector(1., 2., 3.);
        assert_eq!(point.to_vector(), expected);
    }

    // Scenario: Reflecting a vector approaching at 45°
    //   Given v ← vector(1, -1, 0)
    //     And n ← vector(0, 1, 0)
    //   When r ← reflect(v, n)
    //   Then r = vector(1, 1, 0)
    #[test]
    fn reflect_vector_approaching_at_45_degrees() {
        let vector = Tuple::vector(1., -1., 0.);
        let normal = Tuple::vector(0., 1., 0.);
        let ray = vector.reflect(&normal);
        let expected = Tuple::vector(1., 1., 0.);
        assert_eq!(ray, expected);
    }

    // Scenario: Reflecting a vector off a slanted surface
    //   Given v ← vector(0, -1, 0)
    //     And n ← vector(√2/2, √2/2, 0)
    //   When r ← reflect(v, n)
    //   Then r = vector(1, 0, 0)
    #[test]
    fn reflect_vector_off_slanted_surface() {
        let frac_sqrt_2_2 = SQRT_2 / 2.;
        let vector = Tuple::vector(0., -1., 0.);
        let normal = Tuple::vector(frac_sqrt_2_2, frac_sqrt_2_2, 0.);
        let ray = vector.reflect(&normal);
        let expected = Tuple::vector(1., 0., 0.);
        assert_relative_eq!(ray, expected);
    }
}
