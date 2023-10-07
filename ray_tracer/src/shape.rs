use crate::{Intersection, Material, Ray, Transform};
use dyn_clone::DynClone;
use nalgebra_glm::{inverse, DVec4};

pub trait Shape: std::fmt::Debug + DynClone + Transform {
    fn material(&self) -> &Material;
    fn material_mut(&mut self) -> &mut Material;
    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection>;
    fn local_normal_at(&self, point: &DVec4) -> DVec4;
    fn intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let ray = ray.transform(&inverse(self.transform()));
        self.local_intersect(&ray)
    }
    fn normal_at(&self, point: &DVec4) -> DVec4 {
        let local_point = inverse(self.transform()) * point;
        let local_normal = self.local_normal_at(&local_point);
        let mut world_normal = inverse(self.transform()).transpose() * local_normal;
        world_normal.w = 0.;
        world_normal.normalize()
    }
    fn id(&self) -> u32;
}

dyn_clone::clone_trait_object!(Shape);

impl PartialEq for dyn Shape {
    fn eq(&self, other: &dyn Shape) -> bool {
        self.id() == other.id()
    }
}

impl PartialEq<Box<dyn Shape>> for dyn Shape {
    fn eq(&self, other: &Box<dyn Shape>) -> bool {
        self.id() == other.id()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use nalgebra_glm::{rotate_z, scale, translate, vec3, vec4, DMat4};
    use std::f64::consts::{FRAC_1_SQRT_2, SQRT_2};

    #[derive(Debug, PartialEq, Clone)]
    struct TestShape {
        transform: DMat4,
        material: Material,
        id: u32,
    }

    static mut RAY: Option<Ray> = None;

    impl Transform for TestShape {
        fn transform(&self) -> &DMat4 {
            &self.transform
        }

        fn transform_mut(&mut self) -> &mut DMat4 {
            &mut self.transform
        }
    }

    impl Shape for TestShape {
        fn material(&self) -> &Material {
            &self.material
        }

        fn material_mut(&mut self) -> &mut Material {
            &mut self.material
        }

        fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
            unsafe { RAY = Some(ray.clone()) };
            vec![]
        }
        fn local_normal_at(&self, point: &DVec4) -> DVec4 {
            *point
        }

        fn id(&self) -> u32 {
            self.id
        }
    }

    #[test]
    fn intersecting_scaled_shape_with_ray() {
        let shape = TestShape {
            transform: scale(&DMat4::identity(), &vec3(2., 2., 2.)),
            material: Material::default(),
            id: 0,
        };
        let ray = Ray::new(vec4(0., 0., -5., 1.), vec4(0., 0., 1., 0.));
        shape.intersect(&ray);
        unsafe {
            assert_eq!(
                RAY,
                Some(Ray::new(vec4(0., 0., -2.5, 1.), vec4(0., 0., 0.5, 0.)))
            );
        }
    }

    #[test]
    fn intersecting_translated_shape_with_ray() {
        let shape = TestShape {
            transform: translate(&DMat4::identity(), &vec3(5., 0., 0.)),
            material: Material::default(),
            id: 0,
        };
        let ray = Ray::new(vec4(0., 0., -5., 1.), vec4(0., 0., 1., 0.));
        shape.intersect(&ray);
        unsafe {
            assert_eq!(
                RAY,
                Some(Ray::new(vec4(-5., 0., -5., 1.), vec4(0., 0., 1., 0.)))
            );
        }
    }

    #[test]
    fn computing_normal_on_translated_shape() {
        let shape = TestShape {
            transform: translate(&DMat4::identity(), &vec3(0., 1., 0.)),
            material: Material::default(),
            id: 0,
        };
        let normal = shape.normal_at(&vec4(0., FRAC_1_SQRT_2 + 1., -FRAC_1_SQRT_2, 1.));
        assert_relative_eq!(
            normal,
            vec4(0., FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.),
            epsilon = 1e-5f64
        );
    }

    #[test]
    fn computing_normal_on_transformed_shape() {
        let shape = TestShape {
            transform: scale(&DMat4::identity(), &vec3(1., 0.5, 1.))
                * rotate_z(&DMat4::identity(), std::f64::consts::PI / 5.),
            material: Material::default(),
            id: 0,
        };
        let normal = shape.normal_at(&vec4(0., SQRT_2 / 2., -SQRT_2 / 2., 1.));
        assert_relative_eq!(normal, vec4(0., 0.97014, -0.24254, 0.), epsilon = 1e-5f64);
    }
}
