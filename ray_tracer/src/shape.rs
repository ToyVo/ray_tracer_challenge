use crate::{Intersection, Material, Ray, Transform, Tuple};
use dyn_clone::DynClone;

pub trait Shape: std::fmt::Debug + DynClone + Transform {
    fn material(&self) -> &Material;
    fn material_mut(&mut self) -> &mut Material;
    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection>;
    fn local_normal_at(&self, point: &Tuple) -> Tuple;
    fn intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let ray = ray.transform(&self.transform().inverse());
        self.local_intersect(&ray)
    }
    fn normal_at(&self, point: &Tuple) -> Tuple {
        let local_point = &self.transform().inverse() * point;
        let local_normal = self.local_normal_at(&local_point);
        let world_normal = &self.transform().inverse().transpose() * &local_normal;
        world_normal.to_vector().normalize()
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
    use crate::Matrix;
    use std::f64::consts::SQRT_2;

    #[derive(Debug, PartialEq, Clone)]
    struct TestShape {
        transform: Matrix,
        material: Material,
        id: u32,
    }

    static mut RAY: Option<Ray> = None;

    impl Transform for TestShape {
        fn transform(&self) -> &Matrix {
            &self.transform
        }

        fn transform_mut(&mut self) -> &mut Matrix {
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
        fn local_normal_at(&self, point: &Tuple) -> Tuple {
            point.clone()
        }

        fn id(&self) -> u32 {
            self.id
        }
    }

    #[test]
    fn intersecting_scaled_shape_with_ray() {
        let shape = TestShape {
            transform: Matrix::scaling(2., 2., 2.),
            material: Material::new(),
            id: 0,
        };
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        shape.intersect(&ray);
        unsafe {
            assert_eq!(
                RAY,
                Some(Ray::new(
                    Tuple::point(0., 0., -2.5),
                    Tuple::vector(0., 0., 0.5)
                ))
            );
        }
    }

    #[test]
    fn intersecting_translated_shape_with_ray() {
        let shape = TestShape {
            transform: Matrix::translation(5., 0., 0.),
            material: Material::new(),
            id: 0,
        };
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        shape.intersect(&ray);
        unsafe {
            assert_eq!(
                RAY,
                Some(Ray::new(
                    Tuple::point(-5., 0., -5.),
                    Tuple::vector(0., 0., 1.)
                ))
            );
        }
    }

    #[test]
    fn computing_normal_on_translated_shape() {
        let shape = TestShape {
            transform: Matrix::translation(0., 1., 0.),
            material: Material::new(),
            id: 0,
        };
        let normal = shape.normal_at(&Tuple::point(0., 1.70711, -0.70711));
        assert!(normal.nearly_equals(&Tuple::vector(0., 0.70711, -0.70711), 1e-5f64));
    }

    #[test]
    fn computing_normal_on_transformed_shape() {
        let shape = TestShape {
            transform: Matrix::scaling(1., 0.5, 1.) * Matrix::rotation_z(std::f64::consts::PI / 5.),
            material: Material::new(),
            id: 0,
        };
        let normal = shape.normal_at(&Tuple::point(0., SQRT_2 / 2., -SQRT_2 / 2.));
        assert!(normal.nearly_equals(&Tuple::vector(0., 0.97014, -0.24254), 1e-5f64));
    }
}
