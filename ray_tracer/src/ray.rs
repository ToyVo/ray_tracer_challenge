use nalgebra_glm::{DMat4, DVec4};

#[derive(PartialEq, Debug, Clone)]
pub struct Ray {
    pub origin: DVec4,
    pub direction: DVec4,
}

impl Ray {
    pub fn new(origin: DVec4, direction: DVec4) -> Ray {
        Ray { origin, direction }
    }

    pub fn position(&self, t: f64) -> DVec4 {
        self.origin + self.direction * t
    }

    pub fn transform(&self, m: &DMat4) -> Ray {
        Ray::new(m * self.origin, m * self.direction)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra_glm::{scale, translate, vec3, vec4};

    #[test]
    fn creating_and_querying_ray() {
        let ray = Ray::new(vec4(1.0, 2.0, 3.0, 1.), vec4(4.0, 5.0, 6.0, 0.));
        assert_eq!(ray.origin, vec4(1.0, 2.0, 3.0, 1.));
        assert_eq!(ray.direction, vec4(4.0, 5.0, 6.0, 0.));
    }

    #[test]
    fn computing_point_from_distance() {
        let ray = Ray::new(vec4(2.0, 3.0, 4.0, 1.), vec4(1.0, 0.0, 0.0, 0.));
        assert_eq!(ray.position(0.0), vec4(2.0, 3.0, 4.0, 1.));
        assert_eq!(ray.position(1.0), vec4(3.0, 3.0, 4.0, 1.));
        assert_eq!(ray.position(-1.0), vec4(1.0, 3.0, 4.0, 1.));
        assert_eq!(ray.position(2.5), vec4(4.5, 3.0, 4.0, 1.));
    }

    #[test]
    fn translating_ray() {
        let ray = Ray::new(vec4(1.0, 2.0, 3.0, 1.), vec4(0.0, 1.0, 0.0, 0.));
        let translation = translate(&DMat4::identity(), &vec3(3.0, 4.0, 5.0));
        let result = ray.transform(&translation);
        assert_eq!(result.origin, vec4(4.0, 6.0, 8.0, 1.));
        assert_eq!(result.direction, vec4(0.0, 1.0, 0.0, 0.));
    }

    #[test]
    fn scaling_ray() {
        let ray = Ray::new(vec4(1.0, 2.0, 3.0, 1.), vec4(0.0, 1.0, 0.0, 0.));
        let scaling = scale(&DMat4::identity(), &vec3(2.0, 3.0, 4.0));
        let result = ray.transform(&scaling);
        assert_eq!(result.origin, vec4(2.0, 6.0, 12.0, 1.));
        assert_eq!(result.direction, vec4(0.0, 3.0, 0.0, 0.));
    }
}
