use crate::{Ray, Shape, Tuple};

#[derive(Clone, Debug)]
pub struct Intersection {
    pub t: f64,
    pub object: Box<dyn Shape>,
}

impl PartialEq for Intersection {
    fn eq(&self, other: &Self) -> bool {
        self.t == other.t
            && self.object.transform() == other.object.transform()
            && self.object.material() == other.object.material()
    }
}

#[derive(Debug)]
pub struct Computations {
    pub t: f64,
    pub object: Box<dyn Shape>,
    pub point: Tuple,
    pub eye_vector: Tuple,
    pub normal_vector: Tuple,
    pub inside: bool,
    pub over_point: Tuple,
    pub under_point: Tuple,
    pub reflect_vector: Tuple,
    pub n1: f64,
    pub n2: f64,
}

impl Intersection {
    pub fn new(t: f64, object: Box<dyn Shape>) -> Intersection {
        Intersection { t, object }
    }

    pub fn hit(intersections: &[Intersection]) -> Option<&Intersection> {
        let mut min = f64::INFINITY;
        let mut min_index = None;
        for (index, intersection) in intersections.iter().enumerate() {
            if intersection.t < min && intersection.t >= 0.0 {
                min = intersection.t;
                min_index = Some(index);
            }
        }
        match min_index {
            Some(index) => Some(&intersections[index]),
            None => None,
        }
    }

    // function prepare_computations(intersection, ray)
    //   # instantiate a data structure for storing some precomputed values
    //   comps ← new computations data structure
    //
    //   # copy the intersection's properties, for convenience
    //   comps.t       ← intersection.t
    //   comps.object  ← intersection.object
    //
    //   # precompute some useful values
    //   comps.point   ← position(ray, comps.t)
    //   comps.eyev    ← -ray.direction
    //   comps.normalv ← normal_at(comps.object, comps.point)
    //
    //   return comps
    // end function
    pub fn prepare_computations(
        &self,
        ray: &Ray,
        intersections: Option<&Vec<Intersection>>,
    ) -> Computations {
        let point = ray.position(self.t);
        let eye_vector = -ray.direction.clone();
        let mut normal_vector = self.object.normal_at(&point);
        let inside = normal_vector.dot(&eye_vector) < 0.0;
        if inside {
            normal_vector = -normal_vector;
        }
        let delta = &normal_vector * 1e-10f64;
        let over_point = &point + &delta;
        let under_point = &point - &delta;
        let reflect_vector = ray.direction.reflect(&normal_vector);
        let mut containers: Vec<Box<dyn Shape>> = vec![];
        let mut n1 = 1.;
        let mut n2 = self.object.material().refractive_index;
        if let Some(intersections) = intersections {
            for intersection in intersections {
                if intersection == self {
                    n1 = if containers.is_empty() {
                        1.
                    } else {
                        containers[containers.len() - 1].material().refractive_index
                    };
                };
                let mut found = false;
                for i in 0..containers.len() {
                    if containers[i].id() == intersection.object.id() {
                        containers.remove(i);
                        found = true;
                        break;
                    }
                }
                if !found {
                    containers.push(intersection.object.clone());
                }
                if intersection == self {
                    n2 = if containers.is_empty() {
                        1.
                    } else {
                        containers[containers.len() - 1].material().refractive_index
                    };
                    break;
                };
            }
        };
        Computations {
            t: self.t,
            object: self.object.clone(),
            point,
            eye_vector,
            normal_vector,
            inside,
            over_point,
            under_point,
            reflect_vector,
            n1,
            n2,
        }
    }
}

impl Computations {
    // function schlick(comps)
    //   # find the cosine of the angle between the eye and normal vectors
    //   cos ← dot(comps.eyev, comps.normalv)
    //
    //   # total internal reflection can only occur if n1 > n2
    //   if comps.n1 > comps.n2
    //     n ← comps.n1 / comps.n2
    //     sin2_t = n^2 * (1.0 - cos^2)
    //     return 1.0 if sin2_t > 1.0
    //
    //     # compute cosine of theta_t using trig identity
    //     cos_t ← sqrt(1.0 - sin2_t)
    //
    //     # when n1 > n2, use cos(theta_t) instead
    //     cos ← cos_t
    //   end if
    //
    //   r0 ← ((comps.n1 - comps.n2) / (comps.n1 + comps.n2))^2
    //   return r0 + (1 - r0) * (1 - cos)^5
    // end function
    // function schlick(comps)
    //   # find the cosine of the angle between the eye and normal vectors
    //   cos ← dot(comps.eyev, comps.normalv)
    //
    //   # total internal reflection can only occur if n1 > n2
    //   if comps.n1 > comps.n2
    //     n ← comps.n1 / comps.n2
    //     sin2_t = n^2 * (1.0 - cos^2)
    //     return 1.0 if sin2_t > 1.0
    //   end if
    //
    //   # return anything but 1.0 here, so that the test will fail
    //   # appropriately if something goes wrong.
    //   return 0.0
    // end function
    pub fn schlick(&self) -> f64 {
        let mut cos = self.eye_vector.dot(&self.normal_vector);
        if self.n1 > self.n2 {
            let n = self.n1 / self.n2;
            let sin2_t = n.powi(2) * (1. - cos.powi(2));
            if sin2_t > 1. {
                return 1.;
            }
            let cos_t = (1. - sin2_t).sqrt();
            cos = cos_t;
        }

        let r0 = ((self.n1 - self.n2) / (self.n1 + self.n2)).powi(2);
        r0 + (1. - r0) * (1. - cos).powi(5)
    }
}

// Feature: Intersections
#[cfg(test)]
mod tests {
    use crate::{Matrix, Plane, Ray, Sphere, Transform, Tuple};
    use approx::assert_relative_eq;
    use std::f64::consts::SQRT_2;

    use super::*;

    // Scenario: An intersection encapsulates t and object
    //   Given s ← sphere()
    //   When i ← intersection(3.5, s)
    //   Then i.t = 3.5
    //     And i.object = s
    #[test]
    fn intersection_encapsulates_t_and_object() {
        let shape = Sphere::new(0);
        let intersection = Intersection::new(3.5, Box::new(shape.clone()));
        assert_eq!(intersection.t, 3.5);
        assert!(<dyn Shape>::eq(&shape, &intersection.object));
    }

    // Scenario: Aggregating intersections
    //   Given s ← sphere()
    //     And i1 ← intersection(1, s)
    //     And i2 ← intersection(2, s)
    //   When xs ← intersections(i1, i2)
    //   Then xs.count = 2
    //     And xs[0].t = 1
    //     And xs[1].t = 2
    #[test]
    fn aggregating_intersections() {
        let shape = Box::new(Sphere::new(0));
        let i1 = Intersection::new(1.0, shape.clone());
        let i2 = Intersection::new(2.0, shape);
        let intersections = vec![i1, i2];
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, 1.0);
        assert_eq!(intersections[1].t, 2.0);
    }

    // Scenario: The hit, when all intersections have positive t
    //   Given s ← sphere()
    //     And i1 ← intersection(1, s)
    //     And i2 ← intersection(2, s)
    //     And xs ← intersections(i2, i1)
    //   When i ← hit(xs)
    //   Then i = i1
    #[test]
    fn when_all_intersections_have_positive_t() {
        let shape = Box::new(Sphere::new(0));
        let i1 = Intersection::new(1.0, shape.clone());
        let i2 = Intersection::new(2.0, shape);
        let intersections = vec![i1.clone(), i2];
        let hit = Intersection::hit(&intersections);
        assert_eq!(hit, Some(&i1));
    }

    // Scenario: The hit, when some intersections have negative t
    //   Given s ← sphere()
    //     And i1 ← intersection(-1, s)
    //     And i2 ← intersection(1, s)
    //     And xs ← intersections(i2, i1)
    //   When i ← hit(xs)
    //   Then i = i2
    #[test]
    fn when_some_intersections_have_negative_t() {
        let shape = Box::new(Sphere::new(0));
        let i1 = Intersection::new(-1.0, shape.clone());
        let i2 = Intersection::new(1.0, shape);
        let intersections = vec![i1, i2.clone()];
        let hit = Intersection::hit(&intersections);
        assert_eq!(hit, Some(&i2));
    }

    // Scenario: The hit, when all intersections have negative t
    //   Given s ← sphere()
    //     And i1 ← intersection(-2, s)
    //     And i2 ← intersection(-1, s)
    //     And xs ← intersections(i2, i1)
    //   When i ← hit(xs)
    //   Then i is nothing
    #[test]
    fn when_all_intersections_have_negative_t() {
        let shape = Box::new(Sphere::new(0));
        let i1 = Intersection::new(-2.0, shape.clone());
        let i2 = Intersection::new(-1.0, shape);
        let intersections = vec![i1, i2];
        let hit = Intersection::hit(&intersections);
        assert_eq!(hit, None);
    }

    // Scenario: The hit is always the lowest nonnegative intersection
    //   Given s ← sphere()
    //   And i1 ← intersection(5, s)
    //   And i2 ← intersection(7, s)
    //   And i3 ← intersection(-3, s)
    //   And i4 ← intersection(2, s)
    //   And xs ← intersections(i1, i2, i3, i4)
    // When i ← hit(xs)
    // Then i = i4
    #[test]
    fn hit_always_lowest_nonnegative_intersection() {
        let shape = Box::new(Sphere::new(0));
        let i1 = Intersection::new(5.0, shape.clone());
        let i2 = Intersection::new(7.0, shape.clone());
        let i3 = Intersection::new(-3.0, shape.clone());
        let i4 = Intersection::new(2.0, shape);
        let intersections = vec![i1, i2, i3, i4.clone()];
        let hit = Intersection::hit(&intersections);
        assert_eq!(hit, Some(&i4));
    }

    // Scenario: Precomputing the state of an intersection
    //   Given r ← ray(point(0, 0, -5), vector(0, 0, 1))
    //     And shape ← sphere()
    //     And i ← intersection(4, shape)
    //   When comps ← prepare_computations(i, r)
    //   Then comps.t = i.t
    //     And comps.object = i.object
    //     And comps.point = point(0, 0, -1)
    //     And comps.eyev = vector(0, 0, -1)
    //     And comps.normalv = vector(0, 0, -1)
    #[test]
    fn precomputing_state_of_intersection() {
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let shape = Sphere::new(0);
        let intersection = Intersection::new(4., Box::new(shape.clone()));
        let comps = intersection.prepare_computations(&ray, None);
        assert_eq!(comps.t, intersection.t);
        assert!(<dyn Shape>::eq(&shape, &intersection.object));
        assert!(<dyn Shape>::eq(&shape, &comps.object));
        assert_eq!(comps.point, Tuple::point(0., 0., -1.));
        assert_eq!(comps.eye_vector, Tuple::vector(0., 0., -1.));
        assert_eq!(comps.normal_vector, Tuple::vector(0., 0., -1.));
        assert!(!comps.inside);
    }

    // Scenario: The hit, when an intersection occurs on the inside
    //   Given r ← ray(point(0, 0, 0), vector(0, 0, 1))
    //     And shape ← sphere()
    //     And i ← intersection(1, shape)
    //   When comps ← prepare_computations(i, r)
    //   Then comps.point = point(0, 0, 1)
    //     And comps.eyev = vector(0, 0, -1)
    //     And comps.inside = true
    //       # normal would have been (0, 0, 1), but is inverted!
    //     And comps.normalv = vector(0, 0, -1)
    #[test]
    fn hit_when_intersection_inside() {
        let ray = Ray::new(Tuple::point(0., 0., 0.), Tuple::vector(0., 0., 1.));
        let shape = Box::new(Sphere::new(0));
        let intersection = Intersection::new(1., shape);
        let comps = intersection.prepare_computations(&ray, None);
        assert_eq!(comps.point, Tuple::point(0., 0., 1.));
        assert_eq!(comps.eye_vector, Tuple::vector(0., 0., -1.));
        assert_eq!(comps.normal_vector, Tuple::vector(0., 0., -1.));
        assert!(comps.inside);
    }

    // Scenario: The hit should offset the point
    //   Given r ← ray(point(0, 0, -5), vector(0, 0, 1))
    //     And shape ← sphere() with:
    //       | transform | translation(0, 0, 1) |
    //     And i ← intersection(5, shape)
    //   When comps ← prepare_computations(i, r)
    //   Then comps.over_point.z < -EPSILON/2
    //     And comps.point.z > comps.over_point.z
    #[test]
    fn hit_should_offset_point() {
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let mut shape = Sphere::new(0);
        *shape.transform_mut() = Matrix::translation(0., 0., 1.);
        let intersection = Intersection::new(5., Box::new(shape));
        let comps = intersection.prepare_computations(&ray, None);
        assert!(comps.over_point.z() < -f64::EPSILON / 2.);
        assert!(comps.point.z() > comps.over_point.z());
    }

    // Scenario: Precomputing the reflection vector
    //   Given shape ← plane()
    //     And r ← ray(point(0, 1, -1), vector(0, -√2/2, √2/2))
    //     And i ← intersection(√2, shape)
    //   When comps ← prepare_computations(i, r)
    //   Then comps.reflectv = vector(0, √2/2, √2/2)
    #[test]
    fn precomputing_reflection_vector() {
        let shape = Plane::new(0);
        let ray = Ray::new(
            Tuple::point(0., 1., -1.),
            Tuple::vector(0., -SQRT_2 / 2., SQRT_2 / 2.),
        );
        let intersection = Intersection::new(SQRT_2, Box::new(shape));
        let comps = intersection.prepare_computations(&ray, None);
        assert_eq!(
            comps.reflect_vector,
            Tuple::vector(0., SQRT_2 / 2., SQRT_2 / 2.)
        );
    }

    // Scenario Outline: Finding n1 and n2 at various intersections
    //   Given A ← glass_sphere() with:
    //       | transform                 | scaling(2, 2, 2) |
    //       | material.refractive_index | 1.5              |
    //     And B ← glass_sphere() with:
    //       | transform                 | translation(0, 0, -0.25) |
    //       | material.refractive_index | 2.0                      |
    //     And C ← glass_sphere() with:
    //       | transform                 | translation(0, 0, 0.25) |
    //       | material.refractive_index | 2.5                     |
    //     And r ← ray(point(0, 0, -4), vector(0, 0, 1))
    //     And xs ← intersections(2:A, 2.75:B, 3.25:C, 4.75:B, 5.25:C, 6:A)
    //   When comps ← prepare_computations(xs[<index>], r, xs)
    //   Then comps.n1 = <n1>
    //     And comps.n2 = <n2>
    //
    //   Examples:
    //     | index | n1  | n2  |
    //     | 0     | 1.0 | 1.5 |
    //     | 1     | 1.5 | 2.0 |
    //     | 2     | 2.0 | 2.5 |
    //     | 3     | 2.5 | 2.5 |
    //     | 4     | 2.5 | 1.5 |
    //     | 5     | 1.5 | 1.0 |
    #[test]
    fn finding_n1_and_n2_at_various_intersections() {
        let mut a = Box::new(Sphere::glass(0));
        *a.transform_mut() = Matrix::scaling(2., 2., 2.);
        a.material_mut().refractive_index = 1.5;
        let mut b = Box::new(Sphere::glass(1));
        *b.transform_mut() = Matrix::translation(0., 0., -0.25);
        b.material_mut().refractive_index = 2.;
        let mut c = Box::new(Sphere::glass(2));
        *c.transform_mut() = Matrix::translation(0., 0., 0.25);
        c.material_mut().refractive_index = 2.5;
        let ray = Ray::new(Tuple::point(0., 0., -4.), Tuple::vector(0., 0., 1.));
        let intersections = vec![
            Intersection::new(2., a.clone()),
            Intersection::new(2.75, b.clone()),
            Intersection::new(3.25, c.clone()),
            Intersection::new(4.75, b),
            Intersection::new(5.25, c),
            Intersection::new(6., a),
        ];
        let expected_n1 = [1., 1.5, 2., 2.5, 2.5, 1.5];
        let expected_n2 = [1.5, 2., 2.5, 2.5, 1.5, 1.];
        for (index, intersection) in intersections.iter().enumerate() {
            let comps = intersection.prepare_computations(&ray, Some(&intersections));
            assert_eq!(comps.n1, expected_n1[index]);
            assert_eq!(comps.n2, expected_n2[index]);
        }
    }

    // Scenario: The under point is offset below the surface
    //   Given r ← ray(point(0, 0, -5), vector(0, 0, 1))
    //     And shape ← glass_sphere() with:
    //       | transform | translation(0, 0, 1) |
    //     And i ← intersection(5, shape)
    //     And xs ← intersections(i)
    //   When comps ← prepare_computations(i, r, xs)
    //   Then comps.under_point.z > EPSILON/2
    //     And comps.point.z < comps.under_point.z
    #[test]
    fn under_point_is_offset_below_surface() {
        let mut shape = Sphere::glass(0);
        *shape.transform_mut() = Matrix::translation(0., 0., 1.);
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let intersection = Intersection::new(5., Box::new(shape));
        let intersections = vec![intersection.clone()];
        let comps = intersection.prepare_computations(&ray, Some(&intersections));
        assert!(comps.under_point.z() > f64::EPSILON / 2.);
        assert!(comps.point.z() < comps.under_point.z());
    }

    // Scenario: The Schlick approximation under total internal reflection
    //   Given shape ← glass_sphere()
    //     And r ← ray(point(0, 0, √2/2), vector(0, 1, 0))
    //     And xs ← intersections(-√2/2:shape, √2/2:shape)
    //   When comps ← prepare_computations(xs[1], r, xs)
    //     And reflectance ← schlick(comps)
    //   Then reflectance = 1.0
    #[test]
    fn schlick_approximation_under_total_reflection() {
        let shape = Box::new(Sphere::glass(0));
        let ray = Ray::new(Tuple::point(0., 0., SQRT_2 / 2.), Tuple::vector(0., 1., 0.));
        let intersections = vec![
            Intersection::new(-SQRT_2 / 2., shape.clone()),
            Intersection::new(SQRT_2 / 2., shape),
        ];
        let comps = intersections[1].prepare_computations(&ray, Some(&intersections));
        let reflectance = comps.schlick();
        assert_eq!(reflectance, 1.);
    }

    // Scenario: The Schlick approximation with a perpendicular viewing angle
    //   Given shape ← glass_sphere()
    //     And r ← ray(point(0, 0, 0), vector(0, 1, 0))
    //     And xs ← intersections(-1:shape, 1:shape)
    //   When comps ← prepare_computations(xs[1], r, xs)
    //     And reflectance ← schlick(comps)
    //   Then reflectance = 0.04
    #[test]
    fn schlick_approximation_with_perpendicular_viewing_angle() {
        let shape = Box::new(Sphere::glass(0));
        let ray = Ray::new(Tuple::point(0., 0., 0.), Tuple::vector(0., 1., 0.));
        let intersections = vec![
            Intersection::new(-1., shape.clone()),
            Intersection::new(1., shape),
        ];
        let comps = intersections[1].prepare_computations(&ray, Some(&intersections));
        let reflectance = comps.schlick();
        assert_relative_eq!(reflectance, 0.04);
    }

    // Scenario: The Schlick approximation with small angle and n2 > n1
    //   Given shape ← glass_sphere()
    //     And r ← ray(point(0, 0.99, -2), vector(0, 0, 1))
    //     And xs ← intersections(1.8589:shape)
    //   When comps ← prepare_computations(xs[0], r, xs)
    //     And reflectance ← schlick(comps)
    //   Then reflectance = 0.48873
    #[test]
    fn schlick_approximation_with_small_angle_and_n2_greater_than_n1() {
        let shape = Box::new(Sphere::glass(0));
        let ray = Ray::new(Tuple::point(0., 0.99, -2.), Tuple::vector(0., 0., 1.));
        let intersections = vec![Intersection::new(1.8589, shape)];
        let comps = intersections[0].prepare_computations(&ray, Some(&intersections));
        let reflectance = comps.schlick();
        assert!(reflectance - 0.48873 < 1e-5f64);
    }

    // Scenario: The hit, when an intersection occurs on the outside
    //   Given r ← ray(point(0, 0, -5), vector(0, 0, 1))
    //     And shape ← sphere()
    //     And i ← intersection(4, shape)
    //   When comps ← prepare_computations(i, r)
    //   Then comps.inside = false

    // Scenario: An intersection can encapsulate `u` and `v`
    //   Given s ← triangle(point(0, 1, 0), point(-1, 0, 0), point(1, 0, 0))
    //   When i ← intersection_with_uv(3.5, s, 0.2, 0.4)
    //   Then i.u = 0.2
    //     And i.v = 0.4
}
