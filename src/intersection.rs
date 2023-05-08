use crate::Sphere;

#[derive(Clone, Debug, PartialEq)]
pub struct Intersection<'a> {
    pub t: f64,
    pub object: &'a Sphere,
}

impl<'a> Intersection<'a> {
    pub fn new(t: f64, object: &'a Sphere) -> Intersection<'a> {
        Intersection { t, object }
    }

    pub fn hit(xs: &'a Vec<Intersection>) -> Option<&'a Intersection<'a>> {
        let mut min = std::f64::INFINITY;
        let mut min_index = None;
        for (i, x) in xs.iter().enumerate() {
            if x.t < min && x.t >= 0.0 {
                min = x.t;
                min_index = Some(i);
            }
        }
        match min_index {
            Some(i) => Some(&xs[i]),
            None => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intersection_encapsulates_t_and_object() {
        let s = Sphere::new();
        let i = Intersection::new(3.5, &s);
        assert_eq!(i.t, 3.5);
        assert_eq!(i.object, &s);
    }

    #[test]
    fn aggregating_intersections() {
        let s = Sphere::new();
        let i1 = Intersection::new(1.0, &s);
        let i2 = Intersection::new(2.0, &s);
        let xs = vec![i1, i2];
        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].t, 1.0);
        assert_eq!(xs[1].t, 2.0);
    }

    #[test]
    fn when_all_intersections_have_positive_t() {
        let s = Sphere::new();
        let i1 = Intersection::new(1.0, &s);
        let i2 = Intersection::new(2.0, &s);
        let xs = vec![i1.clone(), i2];
        let i = Intersection::hit(&xs);
        assert_eq!(i, Some(&i1));
    }

    #[test]
    fn when_some_intersections_have_negative_t() {
        let s = Sphere::new();
        let i1 = Intersection::new(-1.0, &s);
        let i2 = Intersection::new(1.0, &s);
        let xs = vec![i1, i2.clone()];
        let i = Intersection::hit(&xs);
        assert_eq!(i, Some(&i2));
    }

    #[test]
    fn when_all_intersections_have_negative_t() {
        let s = Sphere::new();
        let i1 = Intersection::new(-2.0, &s);
        let i2 = Intersection::new(-1.0, &s);
        let xs = vec![i1, i2];
        let i = Intersection::hit(&xs);
        assert_eq!(i, None);
    }

    #[test]
    fn hit_always_lowest_nonnegative_intersection() {
        let s = Sphere::new();
        let i1 = Intersection::new(5.0, &s);
        let i2 = Intersection::new(7.0, &s);
        let i3 = Intersection::new(-3.0, &s);
        let i4 = Intersection::new(2.0, &s);
        let xs = vec![i1, i2, i3, i4.clone()];
        let i = Intersection::hit(&xs);
        assert_eq!(i, Some(&i4));
    }
}
