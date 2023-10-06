use crate::{Intersection, Material, Matrix, Ray, Shape, Transform, Tuple};

#[derive(PartialEq, Debug, Clone)]
pub struct Cone {
    transform: Matrix,
    material: Material,
    id: u32,
    maximum: f64,
    minimum: f64,
    closed: bool,
}

impl Cone {
    pub fn new(id: u32) -> Cone {
        Cone {
            transform: Matrix::identity(4),
            material: Material::new(),
            id,
            maximum: f64::INFINITY,
            minimum: f64::NEG_INFINITY,
            closed: false,
        }
    }

    fn check_caps(ray: &Ray, t: f64, y: f64) -> bool {
        let x = ray.origin.x() + t * ray.direction.x();
        let z = ray.origin.z() + t * ray.direction.z();
        (x.powi(2) + z.powi(2)) <= y.abs()
    }

    fn intersect_caps(&self, ray: &Ray) -> Vec<Intersection> {
        let mut intersections = vec![];
        if !self.closed || ray.direction.y().abs() < f64::EPSILON {
            return intersections;
        }

        let t = (self.minimum - ray.origin.y()) / ray.direction.y();
        if Self::check_caps(ray, t, self.minimum) {
            intersections.push(Intersection::new(t, Box::new(self.clone())));
        }

        let t = (self.maximum - ray.origin.y()) / ray.direction.y();
        if Self::check_caps(ray, t, self.maximum) {
            intersections.push(Intersection::new(t, Box::new(self.clone())));
        }

        intersections
    }
}

impl Shape for Cone {
    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let a = ray.direction.x().powi(2) - ray.direction.y().powi(2) + ray.direction.z().powi(2);
        let b = 2.0 * ray.origin.x() * ray.direction.x() - 2.0 * ray.origin.y() * ray.direction.y()
            + 2.0 * ray.origin.z() * ray.direction.z();

        if a.abs() < f64::EPSILON && b.abs() < f64::EPSILON {
            return self.intersect_caps(ray);
        }

        let c = ray.origin.x().powi(2) - ray.origin.y().powi(2) + ray.origin.z().powi(2);

        if a.abs() < f64::EPSILON && b.abs() >= f64::EPSILON {
            let t = -c / (2.0 * b);
            let mut intersections = vec![Intersection::new(t, Box::new(self.clone()))];
            intersections.append(&mut self.intersect_caps(ray));
            return intersections;
        }

        let discriminant = b.powi(2) - 4.0 * a * c;

        if discriminant < 0.0 {
            return vec![];
        }

        let disc = discriminant.sqrt();
        let t0 = (-b - disc) / (2.0 * a);
        let t1 = (-b + disc) / (2.0 * a);
        let (t0, t1) = if t0 > t1 { (t1, t0) } else { (t0, t1) };
        let shape = Box::new(self.clone());
        let mut intersections = vec![];

        let y0 = ray.origin.y() + t0 * ray.direction.y();
        if self.minimum < y0 && y0 < self.maximum {
            intersections.push(Intersection::new(t0, shape.clone()));
        }

        let y1 = ray.origin.y() + t1 * ray.direction.y();
        if self.minimum < y1 && y1 < self.maximum {
            intersections.push(Intersection::new(t1, shape));
        }

        intersections.append(&mut self.intersect_caps(ray));

        intersections
    }
    fn local_normal_at(&self, point: &Tuple) -> Tuple {
        let dist = point.x().powi(2) + point.z().powi(2);
        if dist < 1.0 && point.y() >= self.maximum - f64::EPSILON {
            Tuple::vector(0.0, 1.0, 0.0)
        } else if dist < 1.0 && point.y() <= self.minimum + f64::EPSILON {
            Tuple::vector(0.0, -1.0, 0.0)
        } else {
            let y = (point.x().powi(2) + point.z().powi(2)).sqrt()
                * if point.y() > 0.0 { -1.0 } else { 1.0 };
            Tuple::vector(point.x(), y, point.z())
        }
    }
    fn material(&self) -> &Material {
        &self.material
    }
    fn material_mut(&mut self) -> &mut Material {
        &mut self.material
    }
    fn id(&self) -> u32 {
        self.id
    }
}

impl Transform for Cone {
    fn transform(&self) -> &Matrix {
        &self.transform
    }
    fn transform_mut(&mut self) -> &mut Matrix {
        &mut self.transform
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn ray_intersects_cone_a() {
        let cone = Cone::new(0);
        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0).normalize());
        let intersections = cone.local_intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, 5.0);
        assert_eq!(intersections[1].t, 5.0);
    }

    #[test]
    fn ray_intersects_cone_b() {
        let cone = Cone::new(0);
        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(1.0, 1.0, 1.0).normalize());
        let intersections = cone.local_intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_relative_eq!(intersections[0].t, 8.66025, epsilon = 1e-5f64);
        assert_relative_eq!(intersections[1].t, 8.66025, epsilon = 1e-5f64);
    }

    #[test]
    fn ray_intersects_cone_c() {
        let cone = Cone::new(0);
        let ray = Ray::new(Tuple::point(1.0, 1.0, -5.0), Tuple::vector(-0.5, -1.0, 1.0).normalize());
        let intersections = cone.local_intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_relative_eq!(intersections[0].t, 4.55006, epsilon = 1e-5f64);
        assert_relative_eq!(intersections[1].t, 49.44994, epsilon = 1e-5f64);
    }

    #[test]
    fn ray_intersects_cone_parallel() {
        let cone = Cone::new(0);
        let ray = Ray::new(Tuple::point(0.0, 0.0, -1.0), Tuple::vector(0.0, 1.0, 1.0).normalize());
        let intersections = cone.local_intersect(&ray);
        assert_eq!(intersections.len(), 1);
        assert_relative_eq!(intersections[0].t, 0.35355, epsilon = 1e-5f64);
    }

    #[test]
    fn ray_intersects_cone_end_caps_a() {
        let mut cone = Cone::new(0);
        cone.minimum = -0.5;
        cone.maximum = 0.5;
        cone.closed = true;
        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 1.0, 0.0).normalize());
        let intersections = cone.local_intersect(&ray);
        assert_eq!(intersections.len(), 0);
    }

    #[test]
    fn ray_intersects_cone_end_caps_b() {
        let mut cone = Cone::new(0);
        cone.minimum = -0.5;
        cone.maximum = 0.5;
        cone.closed = true;
        let ray = Ray::new(Tuple::point(0.0, 0.0, -0.25), Tuple::vector(0.0, 1.0, 1.0).normalize());
        let intersections = cone.local_intersect(&ray);
        assert_eq!(intersections.len(), 2);
    }

    #[test]
    fn ray_intersects_cone_end_caps_c() {
        let mut cone = Cone::new(0);
        cone.minimum = -0.5;
        cone.maximum = 0.5;
        cone.closed = true;
        let ray = Ray::new(Tuple::point(0.0, 0.0, -0.25), Tuple::vector(0.0, 1.0, 0.0).normalize());
        let intersections = cone.local_intersect(&ray);
        assert_eq!(intersections.len(), 4);
    }

    #[test]
    fn normal_vector_on_cone_a() {
        let cone = Cone::new(0);
        let n = cone.local_normal_at(&Tuple::point(0.0, 0.0, 0.0));
        assert_eq!(n, Tuple::vector(0.0, 0.0, 0.0));
    }

    #[test]
    fn normal_vector_on_cone_b() {
        use std::f64::consts::SQRT_2;
        let cone = Cone::new(0);
        let n = cone.local_normal_at(&Tuple::point(1.0, 1.0, 1.0));
        assert_eq!(n, Tuple::vector(1.0, -SQRT_2, 1.0));
    }

    #[test]
    fn normal_vector_on_cone_c() {
        let cone = Cone::new(0);
        let n = cone.local_normal_at(&Tuple::point(-1.0, -1.0, 0.0));
        assert_eq!(n, Tuple::vector(-1.0, 1.0, 0.0));
    }
}
