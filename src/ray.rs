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

    pub fn position(&self, t: f64) -> Tuple {
        &self.origin + &self.direction * t
    }

    pub fn transform(&self, m: &Matrix) -> Ray {
        Ray::new(m * &self.origin, m * &self.direction)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn creating_and_querying_ray() {
        let ray = Ray::new(Tuple::point(1.0, 2.0, 3.0), Tuple::vector(4.0, 5.0, 6.0));
        assert_eq!(ray.origin, Tuple::point(1.0, 2.0, 3.0));
        assert_eq!(ray.direction, Tuple::vector(4.0, 5.0, 6.0));
    }

    #[test]
    fn computing_point_from_distance() {
        let ray = Ray::new(Tuple::point(2.0, 3.0, 4.0), Tuple::vector(1.0, 0.0, 0.0));
        assert_eq!(ray.position(0.0), Tuple::point(2.0, 3.0, 4.0));
        assert_eq!(ray.position(1.0), Tuple::point(3.0, 3.0, 4.0));
        assert_eq!(ray.position(-1.0), Tuple::point(1.0, 3.0, 4.0));
        assert_eq!(ray.position(2.5), Tuple::point(4.5, 3.0, 4.0));
    }

    #[test]
    fn translating_ray() {
        let ray = Ray::new(Tuple::point(1.0, 2.0, 3.0), Tuple::vector(0.0, 1.0, 0.0));
        let translation = Matrix::translation(3.0, 4.0, 5.0);
        let result = ray.transform(&translation);
        assert_eq!(result.origin, Tuple::point(4.0, 6.0, 8.0));
        assert_eq!(result.direction, Tuple::vector(0.0, 1.0, 0.0));
    }

    #[test]
    fn scaling_ray() {
        let ray = Ray::new(Tuple::point(1.0, 2.0, 3.0), Tuple::vector(0.0, 1.0, 0.0));
        let scaling = Matrix::scaling(2.0, 3.0, 4.0);
        let result = ray.transform(&scaling);
        assert_eq!(result.origin, Tuple::point(2.0, 6.0, 12.0));
        assert_eq!(result.direction, Tuple::vector(0.0, 3.0, 0.0));
    }
}
