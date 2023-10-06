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

    pub fn hit(intersections: &Vec<Intersection>) -> Option<&Intersection> {
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
        let delta = &normal_vector * (f64::EPSILON * 2.);
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

#[cfg(test)]
mod tests {
    use crate::{Matrix, Plane, Ray, Sphere, Transform, Tuple};
    use approx::assert_relative_eq;
    use std::f64::consts::SQRT_2;

    use super::*;

    #[test]
    fn intersection_encapsulates_t_and_object() {
        let shape = Sphere::new(0);
        let intersection = Intersection::new(3.5, Box::new(shape.clone()));
        assert_eq!(intersection.t, 3.5);
        assert!(<dyn Shape>::eq(&shape, &intersection.object));
    }

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

    #[test]
    fn when_all_intersections_have_positive_t() {
        let shape = Box::new(Sphere::new(0));
        let i1 = Intersection::new(1.0, shape.clone());
        let i2 = Intersection::new(2.0, shape);
        let intersections = vec![i1.clone(), i2];
        let hit = Intersection::hit(&intersections);
        assert_eq!(hit, Some(&i1));
    }

    #[test]
    fn when_some_intersections_have_negative_t() {
        let shape = Box::new(Sphere::new(0));
        let i1 = Intersection::new(-1.0, shape.clone());
        let i2 = Intersection::new(1.0, shape);
        let intersections = vec![i1, i2.clone()];
        let hit = Intersection::hit(&intersections);
        assert_eq!(hit, Some(&i2));
    }

    #[test]
    fn when_all_intersections_have_negative_t() {
        let shape = Box::new(Sphere::new(0));
        let i1 = Intersection::new(-2.0, shape.clone());
        let i2 = Intersection::new(-1.0, shape);
        let intersections = vec![i1, i2];
        let hit = Intersection::hit(&intersections);
        assert_eq!(hit, None);
    }

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
        let expected_n1 = vec![1., 1.5, 2., 2.5, 2.5, 1.5];
        let expected_n2 = vec![1.5, 2., 2.5, 2.5, 1.5, 1.];
        for (index, intersection) in intersections.iter().enumerate() {
            let comps = intersection.prepare_computations(&ray, Some(&intersections));
            assert_eq!(comps.n1, expected_n1[index]);
            assert_eq!(comps.n2, expected_n2[index]);
        }
    }

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
        assert_relative_eq!(reflectance, 0.04, epsilon = 1e-5f64);
    }

    #[test]
    fn schlick_approximation_with_small_angle_and_n2_greater_than_n1() {
        let shape = Box::new(Sphere::glass(0));
        let ray = Ray::new(Tuple::point(0., 0.99, -2.), Tuple::vector(0., 0., 1.));
        let intersections = vec![Intersection::new(1.8589, shape)];
        let comps = intersections[0].prepare_computations(&ray, Some(&intersections));
        let reflectance = comps.schlick();
        assert!(reflectance - 0.48873 < 1e-5f64);
    }
}
