use crate::{Ray, Sphere, Tuple};

#[derive(Clone, Debug, PartialEq)]
pub struct Intersection<'a> {
    pub t: f64,
    pub object: &'a Sphere,
}

pub struct Computations<'a> {
    pub t: f64,
    pub object: &'a Sphere,
    pub point: Tuple,
    pub eye_vector: Tuple,
    pub normal_vector: Tuple,
    pub inside: bool,
    pub over_point: Tuple,
}

impl<'a> Intersection<'a> {
    pub fn new(t: f64, object: &'a Sphere) -> Intersection<'a> {
        Intersection { t, object }
    }

    pub fn hit(intersections: &'a Vec<Intersection>) -> Option<&'a Intersection<'a>> {
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

    pub fn prepare_computations(&self, ray: &Ray) -> Computations<'a> {
        let point = ray.position(self.t);
        let eye_vector = -ray.direction.clone();
        let mut normal_vector = self.object.normal_at(&point);
        let inside = normal_vector.dot(&eye_vector) < 0.0;
        if inside {
            normal_vector = -normal_vector;
        }
        let over_point = &point + &normal_vector * 1e-11;
        Computations {
            t: self.t,
            object: self.object,
            point,
            eye_vector,
            normal_vector,
            inside,
            over_point,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{Matrix, Ray, Tuple};

    use super::*;

    #[test]
    fn intersection_encapsulates_t_and_object() {
        let shape = Sphere::new();
        let intersection = Intersection::new(3.5, &shape);
        assert_eq!(intersection.t, 3.5);
        assert_eq!(intersection.object, &shape);
    }

    #[test]
    fn aggregating_intersections() {
        let shape = Sphere::new();
        let i1 = Intersection::new(1.0, &shape);
        let i2 = Intersection::new(2.0, &shape);
        let intersections = vec![i1, i2];
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, 1.0);
        assert_eq!(intersections[1].t, 2.0);
    }

    #[test]
    fn when_all_intersections_have_positive_t() {
        let shape = Sphere::new();
        let i1 = Intersection::new(1.0, &shape);
        let i2 = Intersection::new(2.0, &shape);
        let intersections = vec![i1.clone(), i2];
        let hit = Intersection::hit(&intersections);
        assert_eq!(hit, Some(&i1));
    }

    #[test]
    fn when_some_intersections_have_negative_t() {
        let shape = Sphere::new();
        let i1 = Intersection::new(-1.0, &shape);
        let i2 = Intersection::new(1.0, &shape);
        let intersections = vec![i1, i2.clone()];
        let hit = Intersection::hit(&intersections);
        assert_eq!(hit, Some(&i2));
    }

    #[test]
    fn when_all_intersections_have_negative_t() {
        let shape = Sphere::new();
        let i1 = Intersection::new(-2.0, &shape);
        let i2 = Intersection::new(-1.0, &shape);
        let intersections = vec![i1, i2];
        let hit = Intersection::hit(&intersections);
        assert_eq!(hit, None);
    }

    #[test]
    fn hit_always_lowest_nonnegative_intersection() {
        let shape = Sphere::new();
        let i1 = Intersection::new(5.0, &shape);
        let i2 = Intersection::new(7.0, &shape);
        let i3 = Intersection::new(-3.0, &shape);
        let i4 = Intersection::new(2.0, &shape);
        let intersections = vec![i1, i2, i3, i4.clone()];
        let hit = Intersection::hit(&intersections);
        assert_eq!(hit, Some(&i4));
    }

    #[test]
    fn precomputing_state_of_intersection() {
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let shape = Sphere::new();
        let intersection = Intersection::new(4., &shape);
        let comps = intersection.prepare_computations(&ray);
        assert_eq!(comps.t, intersection.t);
        assert_eq!(comps.object, intersection.object);
        assert_eq!(comps.point, Tuple::point(0., 0., -1.));
        assert_eq!(comps.eye_vector, Tuple::vector(0., 0., -1.));
        assert_eq!(comps.normal_vector, Tuple::vector(0., 0., -1.));
        assert!(!comps.inside);
    }

    #[test]
    fn hit_when_intersection_inside() {
        let ray = Ray::new(Tuple::point(0., 0., 0.), Tuple::vector(0., 0., 1.));
        let shape = Sphere::new();
        let intersection = Intersection::new(1., &shape);
        let comps = intersection.prepare_computations(&ray);
        assert_eq!(comps.point, Tuple::point(0., 0., 1.));
        assert_eq!(comps.eye_vector, Tuple::vector(0., 0., -1.));
        assert_eq!(comps.normal_vector, Tuple::vector(0., 0., -1.));
        assert!(comps.inside);
    }

    #[test]
    fn hit_should_offset_point() {
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let mut shape = Sphere::new();
        shape.transform = Matrix::translation(0., 0., 1.);
        let intersection = Intersection::new(5., &shape);
        let comps = intersection.prepare_computations(&ray);
        assert!(comps.over_point.z() < -f64::EPSILON / 2.);
        assert!(comps.point.z() > comps.over_point.z());
    }
}
