use crate::{Matrix, Tuple};

#[derive(PartialEq, Debug, Clone)]
pub struct Ray {
    pub origin: Tuple,
    pub direction: Tuple,
}

impl Ray {
    pub fn new(origin: Tuple, direction: Tuple) -> Ray {
        Ray { origin, direction }
    }

    // function position(ray, t)
    //   return ray.origin + ray.direction * t
    // end function
    pub fn position(&self, t: f64) -> Tuple {
        &self.origin + &self.direction * t
    }

    pub fn transform(&self, m: &Matrix) -> Ray {
        Ray::new(m * &self.origin, m * &self.direction)
    }
}

// Feature: Rays
#[cfg(test)]
mod tests {
    use super::*;

    // Scenario: Creating and querying a ray
    //   Given origin ← point(1, 2, 3)
    //     And direction ← vector(4, 5, 6)
    //   When r ← ray(origin, direction)
    //   Then r.origin = origin
    //     And r.direction = direction
    #[test]
    fn creating_and_querying_ray() {
        let ray = Ray::new(Tuple::point(1.0, 2.0, 3.0), Tuple::vector(4.0, 5.0, 6.0));
        assert_eq!(ray.origin, Tuple::point(1.0, 2.0, 3.0));
        assert_eq!(ray.direction, Tuple::vector(4.0, 5.0, 6.0));
    }

    // Scenario: Computing a point from a distance
    //   Given r ← ray(point(2, 3, 4), vector(1, 0, 0))
    //   Then position(r, 0) = point(2, 3, 4)
    //     And position(r, 1) = point(3, 3, 4)
    //     And position(r, -1) = point(1, 3, 4)
    //     And position(r, 2.5) = point(4.5, 3, 4)
    #[test]
    fn computing_point_from_distance() {
        let ray = Ray::new(Tuple::point(2.0, 3.0, 4.0), Tuple::vector(1.0, 0.0, 0.0));
        assert_eq!(ray.position(0.0), Tuple::point(2.0, 3.0, 4.0));
        assert_eq!(ray.position(1.0), Tuple::point(3.0, 3.0, 4.0));
        assert_eq!(ray.position(-1.0), Tuple::point(1.0, 3.0, 4.0));
        assert_eq!(ray.position(2.5), Tuple::point(4.5, 3.0, 4.0));
    }

    // Scenario: Translating a ray
    //   Given r ← ray(point(1, 2, 3), vector(0, 1, 0))
    //     And m ← translation(3, 4, 5)
    //   When r2 ← transform(r, m)
    //   Then r2.origin = point(4, 6, 8)
    //     And r2.direction = vector(0, 1, 0)
    #[test]
    fn translating_ray() {
        let ray = Ray::new(Tuple::point(1.0, 2.0, 3.0), Tuple::vector(0.0, 1.0, 0.0));
        let translation = Matrix::translation(3.0, 4.0, 5.0);
        let result = ray.transform(&translation);
        assert_eq!(result.origin, Tuple::point(4.0, 6.0, 8.0));
        assert_eq!(result.direction, Tuple::vector(0.0, 1.0, 0.0));
    }

    // Scenario: Scaling a ray
    //   Given r ← ray(point(1, 2, 3), vector(0, 1, 0))
    //     And m ← scaling(2, 3, 4)
    //   When r2 ← transform(r, m)
    //   Then r2.origin = point(2, 6, 12)
    //     And r2.direction = vector(0, 3, 0)
    #[test]
    fn scaling_ray() {
        let ray = Ray::new(Tuple::point(1.0, 2.0, 3.0), Tuple::vector(0.0, 1.0, 0.0));
        let scaling = Matrix::scaling(2.0, 3.0, 4.0);
        let result = ray.transform(&scaling);
        assert_eq!(result.origin, Tuple::point(2.0, 6.0, 12.0));
        assert_eq!(result.direction, Tuple::vector(0.0, 3.0, 0.0));
    }
}
