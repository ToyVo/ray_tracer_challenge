use crate::{Intersection, Material, Ray, Shape, Transform};
use nalgebra_glm::{vec4, DMat4, DVec4};

#[derive(PartialEq, Debug, Clone)]
pub struct Sphere {
    transform: DMat4,
    material: Material,
    id: u32,
}

impl Sphere {
    pub fn new(id: u32) -> Sphere {
        Sphere {
            transform: DMat4::identity(),
            material: Material::default(),
            id,
        }
    }

    pub fn glass(id: u32) -> Sphere {
        let mut sphere = Sphere::new(id);
        sphere.material_mut().transparency = 1.0;
        sphere.material_mut().refractive_index = 1.5;
        sphere
    }
}

impl Shape for Sphere {
    fn material(&self) -> &Material {
        &self.material
    }
    fn material_mut(&mut self) -> &mut Material {
        &mut self.material
    }
    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let radius = 1.0;
        let sphere_to_ray = ray.origin - vec4(0.0, 0.0, 0.0, 1.);
        let a = ray.direction.dot(&ray.direction);
        let b = 2.0 * ray.direction.dot(&sphere_to_ray);
        let c = sphere_to_ray.dot(&sphere_to_ray) - radius;
        let discriminant = b.powi(2) - 4.0 * a * c;
        let shape = Box::new(self.clone());
        if discriminant < 0.0 {
            vec![]
        } else {
            let t1 = (-b - discriminant.sqrt()) / (2.0 * a);
            let t2 = (-b + discriminant.sqrt()) / (2.0 * a);
            vec![
                Intersection::new(t1, shape.clone()),
                Intersection::new(t2, shape),
            ]
        }
    }
    fn local_normal_at(&self, point: &DVec4) -> DVec4 {
        point - vec4(0.0, 0.0, 0.0, 1.)
    }
    fn id(&self) -> u32 {
        self.id
    }
}

impl Transform for Sphere {
    fn transform(&self) -> &DMat4 {
        &self.transform
    }
    fn transform_mut(&mut self) -> &mut DMat4 {
        &mut self.transform
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use nalgebra_glm::{rotate_z, scale, translate, vec3, vec4, DMat4};
    use std::f64::consts::{FRAC_1_SQRT_2, PI, SQRT_2};

    use super::*;

    #[test]
    fn ray_intersects_sphere_at_two_points() {
        let ray = Ray::new(vec4(0.0, 0.0, -5.0, 1.), vec4(0.0, 0.0, 1.0, 0.));
        let shape = Sphere::new(0);
        let intersections = shape.intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, 4.0);
        assert_eq!(intersections[1].t, 6.0);
    }

    #[test]
    fn ray_intersects_sphere_at_tangent() {
        let ray = Ray::new(vec4(0.0, 1.0, -5.0, 1.), vec4(0.0, 0.0, 1.0, 0.));
        let shape = Sphere::new(0);
        let intersections = shape.intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, 5.0);
        assert_eq!(intersections[1].t, 5.0);
    }

    #[test]
    fn ray_misses_sphere() {
        let ray = Ray::new(vec4(0.0, 2.0, -5.0, 1.), vec4(0.0, 0.0, 1.0, 0.));
        let shape = Sphere::new(0);
        let intersections = shape.intersect(&ray);
        assert_eq!(intersections.len(), 0);
    }

    #[test]
    fn ray_originates_inside_sphere() {
        let ray = Ray::new(vec4(0.0, 0.0, 0.0, 1.), vec4(0.0, 0.0, 1.0, 0.));
        let shape = Sphere::new(0);
        let intersections = shape.intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, -1.0);
        assert_eq!(intersections[1].t, 1.0);
    }

    #[test]
    fn sphere_behind_ray() {
        let ray = Ray::new(vec4(0.0, 0.0, 5.0, 1.), vec4(0.0, 0.0, 1.0, 0.));
        let shape = Sphere::new(0);
        let intersections = shape.intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, -6.0);
        assert_eq!(intersections[1].t, -4.0);
    }

    #[test]
    fn intersect_sets_object_on_intersection() {
        let ray = Ray::new(vec4(0.0, 0.0, -5.0, 1.), vec4(0.0, 0.0, 1.0, 0.));
        let shape = Sphere::new(0);
        let intersections = shape.intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert!(<dyn Shape>::eq(&shape, &intersections[0].object));
        assert!(<dyn Shape>::eq(&shape, &intersections[1].object));
    }

    #[test]
    fn sphere_default_transformation() {
        let shape = Sphere::new(0);
        assert_eq!(shape.transform(), &DMat4::identity());
    }

    #[test]
    fn changing_sphere_transformation() {
        let mut shape = Sphere::new(0);
        *shape.transform_mut() = translate(&DMat4::identity(), &vec3(2.0, 3.0, 4.0));
        assert_eq!(
            shape.transform(),
            &translate(&DMat4::identity(), &vec3(2.0, 3.0, 4.0))
        );
    }

    #[test]
    fn intersecting_scaled_sphere_with_ray() {
        let ray = Ray::new(vec4(0.0, 0.0, -5.0, 1.), vec4(0.0, 0.0, 1.0, 0.));
        let mut shape = Sphere::new(0);
        *shape.transform_mut() = scale(&DMat4::identity(), &vec3(2.0, 2.0, 2.0));
        let intersections = shape.intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, 3.0);
        assert_eq!(intersections[1].t, 7.0);
    }

    #[test]
    fn intersecting_translated_sphere_with_ray() {
        let ray = Ray::new(vec4(0.0, 0.0, -5.0, 1.), vec4(0.0, 0.0, 1.0, 0.));
        let mut shape = Sphere::new(0);
        *shape.transform_mut() = translate(&DMat4::identity(), &vec3(5.0, 0.0, 0.0));
        let intersections = shape.intersect(&ray);
        assert_eq!(intersections.len(), 0);
    }

    #[test]
    fn normal_on_sphere_at_point_on_x_axis() {
        let shape = Sphere::new(0);
        let normal = shape.normal_at(&vec4(1.0, 0.0, 0.0, 1.));
        assert_eq!(normal, vec4(1.0, 0.0, 0.0, 0.));
    }

    #[test]
    fn normal_on_sphere_at_point_on_y_axis() {
        let shape = Sphere::new(0);
        let normal = shape.normal_at(&vec4(0.0, 1.0, 0.0, 1.));
        assert_eq!(normal, vec4(0.0, 1.0, 0.0, 0.));
    }

    #[test]
    fn normal_on_sphere_at_point_on_z_axis() {
        let shape = Sphere::new(0);
        let normal = shape.normal_at(&vec4(0.0, 0.0, 1.0, 1.));
        assert_eq!(normal, vec4(0.0, 0.0, 1.0, 0.));
    }

    #[test]
    fn normal_on_sphere_at_non_axial_point() {
        let shape = Sphere::new(0);
        let normal = shape.normal_at(&vec4(
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
            1.,
        ));
        assert_eq!(
            normal,
            vec4(
                3.0_f64.sqrt() / 3.0,
                3.0_f64.sqrt() / 3.0,
                3.0_f64.sqrt() / 3.0,
                0.
            )
        );
    }

    #[test]
    fn normal_is_normalized() {
        let shape = Sphere::new(0);
        let normal = shape.normal_at(&vec4(
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
            1.,
        ));
        assert_eq!(normal, normal.normalize());
    }

    #[test]
    fn normal_on_translated_sphere() {
        let mut shape = Sphere::new(0);
        *shape.transform_mut() = translate(&DMat4::identity(), &vec3(0.0, 1.0, 0.0));
        let normal = shape.normal_at(&vec4(0.0, FRAC_1_SQRT_2 + 1., -FRAC_1_SQRT_2, 1.));
        assert_relative_eq!(
            normal,
            vec4(0.0, FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.),
            epsilon = 1e-5f64
        );
    }

    #[test]
    fn normal_on_transformed_sphere() {
        let mut shape = Sphere::new(0);
        *shape.transform_mut() = scale(&DMat4::identity(), &vec3(1.0, 0.5, 1.0))
            * rotate_z(&DMat4::identity(), PI / 5.0);
        let normal = shape.normal_at(&vec4(0.0, SQRT_2 / 2.0, -SQRT_2 / 2.0, 0.));
        assert_relative_eq!(normal, vec4(0.0, 0.97014, -0.24254, 0.), epsilon = 1e-5f64);
    }

    #[test]
    fn sphere_has_default_material() {
        let shape = Sphere::new(0);
        assert_eq!(shape.material(), &Material::default());
    }

    #[test]
    fn sphere_may_be_assigned_material() {
        let mut shape = Sphere::new(0);
        let mut material = Material::default();
        material.ambient = 1.0;
        shape.material_mut().ambient = 1.0;
        assert_eq!(shape.material(), &material);
    }

    #[test]
    fn sphere_with_glassy_material() {
        let sphere = Sphere::glass(0);
        assert_eq!(sphere.material().transparency, 1.0);
        assert_eq!(sphere.material().refractive_index, 1.5);
    }
}
