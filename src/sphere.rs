use crate::{Intersection, Matrix, Ray, Tuple, Material};

#[derive(PartialEq, Debug)]
pub struct Sphere {
    pub transform: Matrix,
    pub material: Material,
}

impl Sphere {
    pub fn new() -> Sphere {
        Sphere {
            transform: Matrix::identity(4),
            material: Material::new(),
        }
    }

    pub fn intersects(&self, ray: &Ray) -> Vec<Intersection> {
        let ray = ray.transform(&self.transform.inverse());
        let radius = 1.0;
        let sphere_to_ray = &ray.origin - Tuple::point(0.0, 0.0, 0.0);
        let a = ray.direction.dot(&ray.direction);
        let b = 2.0 * ray.direction.dot(&sphere_to_ray);
        let c = sphere_to_ray.dot(&sphere_to_ray) - radius;
        let discriminant = b.powi(2) - 4.0 * a * c;
        if discriminant < 0.0 {
            vec![]
        } else {
            let t1 = (-b - discriminant.sqrt()) / (2.0 * a);
            let t2 = (-b + discriminant.sqrt()) / (2.0 * a);
            vec![Intersection::new(t1, self), Intersection::new(t2, self)]
        }
    }

    pub fn normal_at(&self, point: &Tuple) -> Tuple {
        let object_point = &self.transform.inverse() * point;
        let object_normal = object_point - Tuple::point(0.0, 0.0, 0.0);
        let world_normal = &self.transform.inverse().transpose() * &object_normal;
        world_normal.to_vector().normalize()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{PI,SQRT_2,FRAC_1_SQRT_2};

    #[test]
    fn ray_intersects_sphere_at_two_points() {
        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new();
        let intersections = shape.intersects(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, 4.0);
        assert_eq!(intersections[1].t, 6.0);
    }

    #[test]
    fn ray_intersects_sphere_at_tangent() {
        let ray = Ray::new(Tuple::point(0.0, 1.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new();
        let intersections = shape.intersects(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, 5.0);
        assert_eq!(intersections[1].t, 5.0);
    }

    #[test]
    fn ray_misses_sphere() {
        let ray = Ray::new(Tuple::point(0.0, 2.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new();
        let intersections = shape.intersects(&ray);
        assert_eq!(intersections.len(), 0);
    }

    #[test]
    fn ray_originates_inside_sphere() {
        let ray = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new();
        let intersections = shape.intersects(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, -1.0);
        assert_eq!(intersections[1].t, 1.0);
    }

    #[test]
    fn sphere_behind_ray() {
        let ray = Ray::new(Tuple::point(0.0, 0.0, 5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new();
        let intersections = shape.intersects(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, -6.0);
        assert_eq!(intersections[1].t, -4.0);
    }

    #[test]
    fn intersect_sets_object_on_intersection() {
        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new();
        let intersections = shape.intersects(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].object, &shape);
        assert_eq!(intersections[1].object, &shape);
    }

    #[test]
    fn sphere_default_transformation() {
        let shape = Sphere::new();
        assert_eq!(shape.transform, Matrix::identity(4));
    }

    #[test]
    fn changing_sphere_transformation() {
        let mut shape = Sphere::new();
        let translation = Matrix::translation(2.0, 3.0, 4.0);
        shape.transform = translation.clone();
        assert_eq!(shape.transform, translation);
    }

    #[test]
    fn intersecting_scaled_sphere_with_ray() {
        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let mut shape = Sphere::new();
        shape.transform=Matrix::scaling(2.0, 2.0, 2.0);
        let intersections = shape.intersects(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, 3.0);
        assert_eq!(intersections[1].t, 7.0);
    }

    #[test]
    fn intersecting_translated_sphere_with_ray() {
        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let mut shape = Sphere::new();
        shape.transform=Matrix::translation(5.0, 0.0, 0.0);
        let intersections = shape.intersects(&ray);
        assert_eq!(intersections.len(), 0);
    }

    #[test]
    fn normal_on_sphere_at_point_on_x_axis() {
        let shape = Sphere::new();
        let normal = shape.normal_at(&Tuple::point(1.0, 0.0, 0.0));
        assert_eq!(normal, Tuple::vector(1.0, 0.0, 0.0));
    }

    #[test]
    fn normal_on_sphere_at_point_on_y_axis() {
        let shape = Sphere::new();
        let normal = shape.normal_at(&Tuple::point(0.0, 1.0, 0.0));
        assert_eq!(normal, Tuple::vector(0.0, 1.0, 0.0));
    }

    #[test]
    fn normal_on_sphere_at_point_on_z_axis() {
        let shape = Sphere::new();
        let normal = shape.normal_at(&Tuple::point(0.0, 0.0, 1.0));
        assert_eq!(normal, Tuple::vector(0.0, 0.0, 1.0));
    }

    #[test]
    fn normal_on_sphere_at_nonaxial_point() {
        let shape = Sphere::new();
        let normal = shape.normal_at(&Tuple::point(
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
        ));
        assert_eq!(
            normal,
            Tuple::vector(
                3.0_f64.sqrt() / 3.0,
                3.0_f64.sqrt() / 3.0,
                3.0_f64.sqrt() / 3.0
            )
        );
    }

    #[test]
    fn normal_is_normalized() {
        let shape = Sphere::new();
        let normal = shape.normal_at(&Tuple::point(
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
        ));
        assert_eq!(normal, normal.normalize());
    }

    #[test]
    fn normal_on_translated_sphere() {
        let mut shape = Sphere::new();
        shape.transform=Matrix::translation(0.0, 1.0, 0.0);
        let normal = shape.normal_at(&Tuple::point(0.0, FRAC_1_SQRT_2 + 1., -FRAC_1_SQRT_2));
        assert!(normal.nearly_equals(&Tuple::vector(0.0, FRAC_1_SQRT_2, -FRAC_1_SQRT_2), 0.00001));
    }

    #[test]
    fn normal_on_transformed_sphere() {
        let mut shape = Sphere::new();
        let transform = Matrix::scaling(1.0, 0.5, 1.0) * Matrix::rotation_z(PI / 5.0);
        shape.transform= transform;
        let normal = shape.normal_at(&Tuple::point(
            0.0,
            SQRT_2 / 2.0,
            -SQRT_2 / 2.0,
        ));
        assert!(normal.nearly_equals(&Tuple::vector(0.0, 0.97014, -0.24254), 0.00001));
    }

    #[test]
    fn sphere_has_default_material() {
        let shape = Sphere::new();
        let material = shape.material;
        assert_eq!(material, Material::new());
    }

    #[test]
    fn sphere_may_be_assigned_material() {
        let mut shape = Sphere::new();
        let mut material = Material::new();
        material.ambient = 1.0;
        shape.material = material.clone();
        assert_eq!(shape.material, material);
    }
}