use crate::{Intersection, Material, Ray, Transform, Tuple};
use dyn_clone::DynClone;

pub trait Shape: std::fmt::Debug + DynClone + Transform {
    fn material(&self) -> &Material;
    fn material_mut(&mut self) -> &mut Material;
    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection>;
    fn local_normal_at(&self, point: &Tuple) -> Tuple;

    // function intersect(shape, ray)
    //   local_ray ← transform(ray, inverse(shape.transform))
    //   return local_intersect(shape, local_ray)
    // end function
    fn intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let ray = ray.transform(&self.transform().inverse());
        self.local_intersect(&ray)
    }

    // function normal_at(shape, point)
    //   local_point  ← inverse(shape.transform) * point
    //   local_normal ← local_normal_at(shape, local_point)
    //   world_normal ← transpose(inverse(shape.transform)) * local_normal
    //   world_normal.w ← 0
    //
    //   return normalize(world_normal)
    // end function
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

// Feature: Abstract Shapes
#[cfg(test)]
mod tests {
    use super::*;
    use crate::Matrix;
    use approx::assert_relative_eq;
    use std::f64::consts::{FRAC_1_SQRT_2, SQRT_2};

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

    // Scenario: The default transformation
    //   Given s ← test_shape()
    //   Then s.transform = identity_matrix

    // Scenario: Assigning a transformation
    //   Given s ← test_shape()
    //   When set_transform(s, translation(2, 3, 4))
    //   Then s.transform = translation(2, 3, 4)

    // Scenario: The default material
    //   Given s ← test_shape()
    //   When m ← s.material
    //   Then m = material()

    // Scenario: Assigning a material
    //   Given s ← test_shape()
    //     And m ← material()
    //     And m.ambient ← 1
    //   When s.material ← m
    //   Then s.material = m

    // Scenario: Intersecting a scaled shape with a ray
    //   Given r ← ray(point(0, 0, -5), vector(0, 0, 1))
    //     And s ← test_shape()
    //   When set_transform(s, scaling(2, 2, 2))
    //     And xs ← intersect(s, r)
    //   Then s.saved_ray.origin = point(0, 0, -2.5)
    //     And s.saved_ray.direction = vector(0, 0, 0.5)
    #[test]
    fn intersecting_scaled_shape_with_ray() {
        let shape = TestShape {
            transform: Matrix::scaling(2., 2., 2.),
            material: Material::default(),
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

    // Scenario: Intersecting a translated shape with a ray
    //   Given r ← ray(point(0, 0, -5), vector(0, 0, 1))
    //     And s ← test_shape()
    //   When set_transform(s, translation(5, 0, 0))
    //     And xs ← intersect(s, r)
    //   Then s.saved_ray.origin = point(-5, 0, -5)
    //     And s.saved_ray.direction = vector(0, 0, 1)
    #[test]
    fn intersecting_translated_shape_with_ray() {
        let shape = TestShape {
            transform: Matrix::translation(5., 0., 0.),
            material: Material::default(),
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

    // Scenario: Computing the normal on a translated shape
    //   Given s ← test_shape()
    //   When set_transform(s, translation(0, 1, 0))
    //     And n ← normal_at(s, point(0, 1.70711, -0.70711))
    //   Then n = vector(0, 0.70711, -0.70711)
    #[test]
    fn computing_normal_on_translated_shape() {
        let shape = TestShape {
            transform: Matrix::translation(0., 1., 0.),
            material: Material::default(),
            id: 0,
        };
        let normal = shape.normal_at(&Tuple::point(0., FRAC_1_SQRT_2 + 1., -FRAC_1_SQRT_2));
        assert_relative_eq!(
            normal,
            Tuple::vector(0., FRAC_1_SQRT_2, -FRAC_1_SQRT_2)
        );
    }

    // Scenario: Computing the normal on a transformed shape
    //   Given s ← test_shape()
    //     And m ← scaling(1, 0.5, 1) * rotation_z(π/5)
    //   When set_transform(s, m)
    //     And n ← normal_at(s, point(0, √2/2, -√2/2))
    //   Then n = vector(0, 0.97014, -0.24254)
    #[test]
    fn computing_normal_on_transformed_shape() {
        let shape = TestShape {
            transform: Matrix::scaling(1., 0.5, 1.) * Matrix::rotation_z(std::f64::consts::PI / 5.),
            material: Material::default(),
            id: 0,
        };
        let normal = shape.normal_at(&Tuple::point(0., SQRT_2 / 2., -SQRT_2 / 2.));
        assert_relative_eq!(
            normal,
            Tuple::vector(0., 0.97014, -0.24254)
        );
    }

    // Scenario: A shape has a parent attribute
    //   Given s ← test_shape()
    //   Then s.parent is nothing

    // Scenario: Converting a point from world to object space
    //   Given g1 ← group()
    //     And set_transform(g1, rotation_y(π/2))
    //     And g2 ← group()
    //     And set_transform(g2, scaling(2, 2, 2))
    //     And add_child(g1, g2)
    //     And s ← sphere()
    //     And set_transform(s, translation(5, 0, 0))
    //     And add_child(g2, s)
    //   When p ← world_to_object(s, point(-2, 0, -10))
    //   Then p = point(0, 0, -1)

    // Scenario: Converting a normal from object to world space
    //   Given g1 ← group()
    //     And set_transform(g1, rotation_y(π/2))
    //     And g2 ← group()
    //     And set_transform(g2, scaling(1, 2, 3))
    //     And add_child(g1, g2)
    //     And s ← sphere()
    //     And set_transform(s, translation(5, 0, 0))
    //     And add_child(g2, s)
    //   When n ← normal_to_world(s, vector(√3/3, √3/3, √3/3))
    //   Then n = vector(0.2857, 0.4286, -0.8571)

    // Scenario: Finding the normal on a child object
    //   Given g1 ← group()
    //     And set_transform(g1, rotation_y(π/2))
    //     And g2 ← group()
    //     And set_transform(g2, scaling(1, 2, 3))
    //     And add_child(g1, g2)
    //     And s ← sphere()
    //     And set_transform(s, translation(5, 0, 0))
    //     And add_child(g2, s)
    //   When n ← normal_at(s, point(1.7321, 1.1547, -5.5774))
    //   Then n = vector(0.2857, 0.4286, -0.8571)
}
