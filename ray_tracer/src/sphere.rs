use crate::{Intersection, Material, Matrix, Ray, Shape, Transform, Tuple};

#[derive(PartialEq, Debug, Clone)]
pub struct Sphere {
    transform: Matrix,
    material: Material,
    id: u32,
}

impl Sphere {
    pub fn new(id: u32) -> Sphere {
        Sphere {
            transform: Matrix::identity(4),
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

    // function intersect(sphere, ray)
    //   # the vector from the sphere's center, to the ray origin
    //   # remember: the sphere is centered at the world origin
    //   sphere_to_ray ← ray.origin - point(0, 0, 0)
    //
    //   a ← dot(ray.direction, ray.direction)
    //   b ← 2 * dot(ray.direction, sphere_to_ray)
    //   c ← dot(sphere_to_ray, sphere_to_ray) - 1
    //
    //   discriminant ← b² - 4 * a * c
    //
    //   if discriminant < 0 then
    //     return ()
    //   end if
    //
    //   t1 ← (-b - √(discriminant)) / (2 * a)
    //   t2 ← (-b + √(discriminant)) / (2 * a)
    //
    //   return (t1, t2)
    // end function
    // function intersect(sphere, ray)
    //   ray2 ← transform(ray, inverse(sphere.transform))
    //
    //   # ...
    // end function
    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let radius = 1.0;
        let sphere_to_ray = &ray.origin - Tuple::point(0.0, 0.0, 0.0);
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
    fn local_normal_at(&self, point: &Tuple) -> Tuple {
        point - Tuple::point(0.0, 0.0, 0.0)
    }
    fn id(&self) -> u32 {
        self.id
    }
}

impl Transform for Sphere {
    fn transform(&self) -> &Matrix {
        &self.transform
    }
    fn transform_mut(&mut self) -> &mut Matrix {
        &mut self.transform
    }
}

// Feature: Spheres
#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use std::f64::consts::{FRAC_1_SQRT_2, PI, SQRT_2};

    use super::*;

    // Scenario: A ray intersects a sphere at two points
    //   Given r ← ray(point(0, 0, -5), vector(0, 0, 1))
    //     And s ← sphere()
    //   When xs ← intersect(s, r)
    //   Then xs.count = 2
    //     And xs[0] = 4.0
    //     And xs[1] = 6.0
    #[test]
    fn ray_intersects_sphere_at_two_points() {
        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new(0);
        let intersections = shape.intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, 4.0);
        assert_eq!(intersections[1].t, 6.0);
    }

    // Scenario: A ray intersects a sphere at a tangent
    //   Given r ← ray(point(0, 1, -5), vector(0, 0, 1))
    //     And s ← sphere()
    //   When xs ← intersect(s, r)
    //   Then xs.count = 2
    //     And xs[0] = 5.0
    //     And xs[1] = 5.0
    #[test]
    fn ray_intersects_sphere_at_tangent() {
        let ray = Ray::new(Tuple::point(0.0, 1.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new(0);
        let intersections = shape.intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, 5.0);
        assert_eq!(intersections[1].t, 5.0);
    }

    // Scenario: A ray misses a sphere
    //   Given r ← ray(point(0, 2, -5), vector(0, 0, 1))
    //     And s ← sphere()
    //   When xs ← intersect(s, r)
    //   Then xs.count = 0
    #[test]
    fn ray_misses_sphere() {
        let ray = Ray::new(Tuple::point(0.0, 2.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new(0);
        let intersections = shape.intersect(&ray);
        assert_eq!(intersections.len(), 0);
    }

    // Scenario: A ray originates inside a sphere
    //   Given r ← ray(point(0, 0, 0), vector(0, 0, 1))
    //     And s ← sphere()
    //   When xs ← intersect(s, r)
    //   Then xs.count = 2
    //     And xs[0] = -1.0
    //     And xs[1] = 1.0
    #[test]
    fn ray_originates_inside_sphere() {
        let ray = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new(0);
        let intersections = shape.intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, -1.0);
        assert_eq!(intersections[1].t, 1.0);
    }

    // Scenario: A sphere is behind a ray
    //   Given r ← ray(point(0, 0, 5), vector(0, 0, 1))
    //     And s ← sphere()
    //   When xs ← intersect(s, r)
    //   Then xs.count = 2
    //     And xs[0] = -6.0
    //     And xs[1] = -4.0
    #[test]
    fn sphere_behind_ray() {
        let ray = Ray::new(Tuple::point(0.0, 0.0, 5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new(0);
        let intersections = shape.intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, -6.0);
        assert_eq!(intersections[1].t, -4.0);
    }

    // Scenario: Intersect sets the object on the intersection
    //   Given r ← ray(point(0, 0, -5), vector(0, 0, 1))
    //     And s ← sphere()
    //   When xs ← intersect(s, r)
    //   Then xs.count = 2
    //     And xs[0].object = s
    //     And xs[1].object = s
    #[test]
    fn intersect_sets_object_on_intersection() {
        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new(0);
        let intersections = shape.intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert!(<dyn Shape>::eq(&shape, &intersections[0].object));
        assert!(<dyn Shape>::eq(&shape, &intersections[1].object));
    }

    // Scenario: A sphere's default transformation
    //   Given s ← sphere()
    //   Then s.transform = identity_matrix
    #[test]
    fn sphere_default_transformation() {
        let shape = Sphere::new(0);
        assert_eq!(shape.transform(), &Matrix::identity(4));
    }

    // Scenario: Changing a sphere's transformation
    //   Given s ← sphere()
    //     And t ← translation(2, 3, 4)
    //   When set_transform(s, t)
    //   Then s.transform = t
    #[test]
    fn changing_sphere_transformation() {
        let mut shape = Sphere::new(0);
        *shape.transform_mut() = Matrix::translation(2.0, 3.0, 4.0);
        assert_eq!(shape.transform(), &Matrix::translation(2.0, 3.0, 4.0));
    }

    // Scenario: Intersecting a scaled sphere with a ray
    //   Given r ← ray(point(0, 0, -5), vector(0, 0, 1))
    //     And s ← sphere()
    //   When set_transform(s, scaling(2, 2, 2))
    //     And xs ← intersect(s, r)
    //   Then xs.count = 2
    //     And xs[0].t = 3
    //     And xs[1].t = 7
    #[test]
    fn intersecting_scaled_sphere_with_ray() {
        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let mut shape = Sphere::new(0);
        *shape.transform_mut() = Matrix::scaling(2.0, 2.0, 2.0);
        let intersections = shape.intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, 3.0);
        assert_eq!(intersections[1].t, 7.0);
    }

    // Scenario: Intersecting a translated sphere with a ray
    //   Given r ← ray(point(0, 0, -5), vector(0, 0, 1))
    //     And s ← sphere()
    //   When set_transform(s, translation(5, 0, 0))
    //     And xs ← intersect(s, r)
    //   Then xs.count = 0
    #[test]
    fn intersecting_translated_sphere_with_ray() {
        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let mut shape = Sphere::new(0);
        *shape.transform_mut() = Matrix::translation(5.0, 0.0, 0.0);
        let intersections = shape.intersect(&ray);
        assert_eq!(intersections.len(), 0);
    }

    // Scenario: The normal on a sphere at a point on the x axis
    //   Given s ← sphere()
    //   When n ← normal_at(s, point(1, 0, 0))
    //   Then n = vector(1, 0, 0)
    #[test]
    fn normal_on_sphere_at_point_on_x_axis() {
        let shape = Sphere::new(0);
        let normal = shape.normal_at(&Tuple::point(1.0, 0.0, 0.0));
        assert_eq!(normal, Tuple::vector(1.0, 0.0, 0.0));
    }

    // Scenario: The normal on a sphere at a point on the y axis
    //   Given s ← sphere()
    //   When n ← normal_at(s, point(0, 1, 0))
    //   Then n = vector(0, 1, 0)
    #[test]
    fn normal_on_sphere_at_point_on_y_axis() {
        let shape = Sphere::new(0);
        let normal = shape.normal_at(&Tuple::point(0.0, 1.0, 0.0));
        assert_eq!(normal, Tuple::vector(0.0, 1.0, 0.0));
    }

    // Scenario: The normal on a sphere at a point on the z axis
    //   Given s ← sphere()
    //   When n ← normal_at(s, point(0, 0, 1))
    //   Then n = vector(0, 0, 1)
    #[test]
    fn normal_on_sphere_at_point_on_z_axis() {
        let shape = Sphere::new(0);
        let normal = shape.normal_at(&Tuple::point(0.0, 0.0, 1.0));
        assert_eq!(normal, Tuple::vector(0.0, 0.0, 1.0));
    }

    // Scenario: The normal on a sphere at a nonaxial point
    //   Given s ← sphere()
    //   When n ← normal_at(s, point(√3/3, √3/3, √3/3))
    //   Then n = vector(√3/3, √3/3, √3/3)
    #[test]
    fn normal_on_sphere_at_nonaxial_point() {
        let shape = Sphere::new(0);
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
                3.0_f64.sqrt() / 3.0,
            )
        );
    }

    // Scenario: The normal is a normalized vector
    //   Given s ← sphere()
    //   When n ← normal_at(s, point(√3/3, √3/3, √3/3))
    //   Then n = normalize(n)
    #[test]
    fn normal_is_normalized() {
        let shape = Sphere::new(0);
        let normal = shape.normal_at(&Tuple::point(
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
        ));
        assert_eq!(normal, normal.normalize());
    }

    // Scenario: Computing the normal on a translated sphere
    //   Given s ← sphere()
    //     And set_transform(s, translation(0, 1, 0))
    //   When n ← normal_at(s, point(0, 1.70711, -0.70711))
    //   Then n = vector(0, 0.70711, -0.70711)
    #[test]
    fn normal_on_translated_sphere() {
        let mut shape = Sphere::new(0);
        *shape.transform_mut() = Matrix::translation(0.0, 1.0, 0.0);
        let normal = shape.normal_at(&Tuple::point(0.0, FRAC_1_SQRT_2 + 1., -FRAC_1_SQRT_2));
        assert_relative_eq!(
            normal,
            Tuple::vector(0.0, FRAC_1_SQRT_2, -FRAC_1_SQRT_2)
        );
    }

    // Scenario: Computing the normal on a transformed sphere
    //   Given s ← sphere()
    //     And m ← scaling(1, 0.5, 1) * rotation_z(π/5)
    //     And set_transform(s, m)
    //   When n ← normal_at(s, point(0, √2/2, -√2/2))
    //   Then n = vector(0, 0.97014, -0.24254)
    #[test]
    fn normal_on_transformed_sphere() {
        let mut shape = Sphere::new(0);
        *shape.transform_mut() = Matrix::scaling(1.0, 0.5, 1.0) * Matrix::rotation_z(PI / 5.0);
        let normal = shape.normal_at(&Tuple::point(0.0, SQRT_2 / 2.0, -SQRT_2 / 2.0));
        assert_relative_eq!(
            normal,
            Tuple::vector(0.0, 0.97014, -0.24254)
        );
    }

    // Scenario: A sphere has a default material
    //   Given s ← sphere()
    //   When m ← s.material
    //   Then m = material()
    #[test]
    fn sphere_has_default_material() {
        let shape = Sphere::new(0);
        assert_eq!(shape.material(), &Material::default());
    }

    // Scenario: A sphere may be assigned a material
    //   Given s ← sphere()
    //     And m ← material()
    //     And m.ambient ← 1
    //   When s.material ← m
    //   Then s.material = m
    #[test]
    fn sphere_may_be_assigned_material() {
        let mut shape = Sphere::new(0);
        let mut material = Material::default();
        material.ambient = 1.0;
        shape.material_mut().ambient = 1.0;
        assert_eq!(shape.material(), &material);
    }

    // Scenario: A helper for producing a sphere with a glassy material
    //   Given s ← glass_sphere()
    //   Then s.transform = identity_matrix
    //     And s.material.transparency = 1.0
    //     And s.material.refractive_index = 1.5
    #[test]
    fn sphere_with_glassy_material() {
        let sphere = Sphere::glass(0);
        assert_eq!(sphere.material().transparency, 1.0);
        assert_eq!(sphere.material().refractive_index, 1.5);
    }
}
