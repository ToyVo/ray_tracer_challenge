use crate::{Intersection, Material, Matrix, Ray, Shape, Transform, Tuple};

#[derive(PartialEq, Debug, Clone)]
pub struct Cylinder {
    transform: Matrix,
    material: Material,
    id: u32,
    maximum: f64,
    minimum: f64,
    closed: bool,
}

impl Cylinder {
    pub fn new(id: u32) -> Cylinder {
        Cylinder {
            transform: Matrix::identity(4),
            material: Material::default(),
            id,
            maximum: f64::INFINITY,
            minimum: f64::NEG_INFINITY,
            closed: false,
        }
    }

    // # a helper function to reduce duplication.
    // # checks to see if the intersection at `t` is within a radius
    // # of 1 (the radius of your cylinders) from the y axis.
    // function check_cap(ray, t)
    //   x ← ray.origin.x + t * ray.direction.x
    //   z ← ray.origin.z + t * ray.direction.z
    //
    //   return (x² + z²) <= 1
    // end
    fn check_caps(ray: &Ray, t: f64) -> bool {
        let x = ray.origin.x() + t * ray.direction.x();
        let z = ray.origin.z() + t * ray.direction.z();
        (x.powi(2) + z.powi(2)) <= 1.0
    }

    // function intersect_caps(cyl, ray, xs)
    //   # caps only matter if the cylinder is closed, and might possibly be
    //   # intersected by the ray.
    //   if cyl is not closed or ray.direction.y is close to zero
    //     return
    //   end if
    //
    //   # check for an intersection with the lower end cap by intersecting
    //   # the ray with the plane at y=cyl.minimum
    //   t ← (cyl.minimum - ray.origin.y) / ray.direction.y
    //   if check_cap(ray, t)
    //     add intersection(t, cyl) to xs
    //   end if
    //
    //   # check for an intersection with the upper end cap by intersecting
    //   # the ray with the plane at y=cyl.maximum
    //   t ← (cyl.maximum - ray.origin.y) / ray.direction.y
    //   if check_cap(ray, t)
    //     add intersection(t, cyl) to xs
    //   end if
    // end function
    fn intersect_caps(&self, ray: &Ray) -> Vec<Intersection> {
        let mut intersections = vec![];
        if !self.closed || ray.direction.y().abs() < f64::EPSILON {
            return intersections;
        }

        let t = (self.minimum - ray.origin.y()) / ray.direction.y();
        if Self::check_caps(ray, t) {
            intersections.push(Intersection::new(t, Box::new(self.clone())));
        }

        let t = (self.maximum - ray.origin.y()) / ray.direction.y();
        if Self::check_caps(ray, t) {
            intersections.push(Intersection::new(t, Box::new(self.clone())));
        }

        intersections
    }
}

impl Shape for Cylinder {
    fn material(&self) -> &Material {
        &self.material
    }
    fn material_mut(&mut self) -> &mut Material {
        &mut self.material
    }

    // function local_intersect(cylinder, ray)
    //   a ← ray.direction.x² + ray.direction.z²
    //
    //   # ray is parallel to the y axis
    //   return () if a is approximately zero
    //
    //   b ← 2 * ray.origin.x * ray.direction.x +
    //        2 * ray.origin.z * ray.direction.z
    //   c ← ray.origin.x² + ray.origin.z² - 1
    //
    //   disc ← b² - 4 * a * c
    //
    //   # ray does not intersect the cylinder
    //   return () if disc < 0
    //
    //   # this is just a placeholder, to ensure the tests
    //   # pass that expect the ray to miss.
    //   return ( intersection(1, cylinder) )
    //
    //   t0 ← (-b - √(disc)) / (2 * a)
    //   t1 ← (-b + √(disc)) / (2 * a)
    //
    //   return ( intersection(t0, cylinder), intersection(t1, cylinder) )
    // end function
    // t0 ← (-b - √(disc)) / (2 * a)
    // t1 ← (-b + √(disc)) / (2 * a)
    // if t0 > t1 then swap(t0, t1)
    //
    // xs = ()
    //
    // y0 ← ray.origin.y + t0 * ray.direction.y
    // if cylinder.minimum < y0 and y0 < cylinder.maximum
    //   add intersection(t0, cylinder) to xs
    // end if
    //
    // y1 ← ray.origin.y + t1 * ray.direction.y
    // if cylinder.minimum < y1 and y1 < cylinder.maximum
    //   add intersection(t1, cylinder) to xs
    // end if
    //
    // return xs
    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let mut intersections = vec![];
        let a = ray.direction.x().powi(2) + ray.direction.z().powi(2);
        let b = 2.0 * ray.origin.x() * ray.direction.x() + 2.0 * ray.origin.z() * ray.direction.z();
        let c = ray.origin.x().powi(2) + ray.origin.z().powi(2) - 1.0;

        if a.abs() >= f64::EPSILON {
            let discriminant = b.powi(2) - 4.0 * a * c;
            if discriminant < 0.0 {
                return intersections;
            }

            let disc = discriminant.sqrt();
            let t0 = (-b - disc) / (2.0 * a);
            let t1 = (-b + disc) / (2.0 * a);
            let (t0, t1) = if t0 > t1 { (t1, t0) } else { (t0, t1) };
            let shape = Box::new(self.clone());

            let y0 = ray.origin.y() + t0 * ray.direction.y();
            if self.minimum < y0 && y0 < self.maximum {
                intersections.push(Intersection::new(t0, shape.clone()));
            }

            let y1 = ray.origin.y() + t1 * ray.direction.y();
            if self.minimum < y1 && y1 < self.maximum {
                intersections.push(Intersection::new(t1, shape));
            }
        }

        intersections.append(&mut self.intersect_caps(ray));

        intersections
    }

    // function local_normal_at(cylinder, point)
    //   # compute the square of the distance from the y axis
    //   dist ← point.x² + point.z²
    //
    //   if dist < 1 and point.y >= cylinder.maximum - EPSILON
    //     return vector(0, 1, 0)
    //
    //   else if dist < 1 and point.y <= cylinder.minimum + EPSILON
    //     return vector(0, -1, 0)
    //
    //   else
    //     return vector(point.x, 0, point.z)
    //   end if
    // end function
    fn local_normal_at(&self, point: &Tuple) -> Tuple {
        let dist = point.x().powi(2) + point.z().powi(2);
        if dist < 1.0 && point.y() >= self.maximum - f64::EPSILON {
            Tuple::vector(0.0, 1.0, 0.0)
        } else if dist < 1.0 && point.y() <= self.minimum + f64::EPSILON {
            Tuple::vector(0.0, -1.0, 0.0)
        } else {
            Tuple::vector(point.x(), 0.0, point.z())
        }
    }
    fn id(&self) -> u32 {
        self.id
    }
}

impl Transform for Cylinder {
    fn transform(&self) -> &Matrix {
        &self.transform
    }
    fn transform_mut(&mut self) -> &mut Matrix {
        &mut self.transform
    }
}

// Feature: Cylinders
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    // Scenario Outline: A ray misses a cylinder
    //   Given cyl ← cylinder()
    //     And direction ← normalize(<direction>)
    //     And r ← ray(<origin>, direction)
    //   When xs ← local_intersect(cyl, r)
    //   Then xs.count = 0
    //
    //   Examples:
    //     | origin          | direction       |
    //     | point(1, 0, 0)  | vector(0, 1, 0) |
    #[test]
    fn ray_misses_cylinder_a() {
        let cylinder = Cylinder::new(0);
        let ray = Ray::new(
            Tuple::point(1.0, 0.0, 0.0),
            Tuple::vector(0.0, 1.0, 0.0).normalize(),
        );
        let intersections = cylinder.local_intersect(&ray);
        assert_eq!(intersections.len(), 0);
    }

    // Scenario Outline: A ray misses a cylinder
    //   Given cyl ← cylinder()
    //     And direction ← normalize(<direction>)
    //     And r ← ray(<origin>, direction)
    //   When xs ← local_intersect(cyl, r)
    //   Then xs.count = 0
    //
    //   Examples:
    //     | origin          | direction       |
    //     | point(0, 0, 0)  | vector(0, 1, 0) |
    #[test]
    fn ray_misses_cylinder_b() {
        let cylinder = Cylinder::new(0);
        let ray = Ray::new(
            Tuple::point(0.0, 0.0, 0.0),
            Tuple::vector(0.0, 1.0, 0.0).normalize(),
        );
        let intersections = cylinder.local_intersect(&ray);
        assert_eq!(intersections.len(), 0);
    }

    // Scenario Outline: A ray misses a cylinder
    //   Given cyl ← cylinder()
    //     And direction ← normalize(<direction>)
    //     And r ← ray(<origin>, direction)
    //   When xs ← local_intersect(cyl, r)
    //   Then xs.count = 0
    //
    //   Examples:
    //     | origin          | direction       |
    //     | point(0, 0, -5) | vector(1, 1, 1) |
    #[test]
    fn ray_misses_cylinder_c() {
        let cylinder = Cylinder::new(0);
        let ray = Ray::new(
            Tuple::point(0.0, 0.0, -5.0),
            Tuple::vector(1.0, 1.0, 1.0).normalize(),
        );
        let intersections = cylinder.local_intersect(&ray);
        assert_eq!(intersections.len(), 0);
    }

    // Scenario Outline: A ray strikes a cylinder
    //   Given cyl ← cylinder()
    //     And direction ← normalize(<direction>)
    //     And r ← ray(<origin>, direction)
    //   When xs ← local_intersect(cyl, r)
    //   Then xs.count = 2
    //     And xs[0].t = <t0>
    //     And xs[1].t = <t1>
    //
    //   Examples:
    //     | origin            | direction         | t0      | t1      |
    //     | point(1, 0, -5)   | vector(0, 0, 1)   | 5       | 5       |
    #[test]
    fn ray_hits_cylinder_a() {
        let cylinder = Cylinder::new(0);
        let ray = Ray::new(
            Tuple::point(1.0, 0.0, -5.0),
            Tuple::vector(0.0, 0.0, 1.0).normalize(),
        );
        let intersections = cylinder.local_intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, 5.0);
        assert_eq!(intersections[1].t, 5.0);
    }

    // Scenario Outline: A ray strikes a cylinder
    //   Given cyl ← cylinder()
    //     And direction ← normalize(<direction>)
    //     And r ← ray(<origin>, direction)
    //   When xs ← local_intersect(cyl, r)
    //   Then xs.count = 2
    //     And xs[0].t = <t0>
    //     And xs[1].t = <t1>
    //
    //   Examples:
    //     | origin            | direction         | t0      | t1      |
    //     | point(0, 0, -5)   | vector(0, 0, 1)   | 4       | 6       |
    #[test]
    fn ray_hits_cylinder_b() {
        let cylinder = Cylinder::new(0);
        let ray = Ray::new(
            Tuple::point(0.0, 0.0, -5.0),
            Tuple::vector(0.0, 0.0, 1.0).normalize(),
        );
        let intersections = cylinder.local_intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, 4.0);
        assert_eq!(intersections[1].t, 6.0);
    }

    // Scenario Outline: A ray strikes a cylinder
    //   Given cyl ← cylinder()
    //     And direction ← normalize(<direction>)
    //     And r ← ray(<origin>, direction)
    //   When xs ← local_intersect(cyl, r)
    //   Then xs.count = 2
    //     And xs[0].t = <t0>
    //     And xs[1].t = <t1>
    //
    //   Examples:
    //     | origin            | direction         | t0      | t1      |
    //     | point(0.5, 0, -5) | vector(0.1, 1, 1) | 6.80798 | 7.08872 |
    #[test]
    fn ray_hits_cylinder_c() {
        let cylinder = Cylinder::new(0);
        let ray = Ray::new(
            Tuple::point(0.5, 0.0, -5.0),
            Tuple::vector(0.1, 1.0, 1.0).normalize(),
        );
        let intersections = cylinder.local_intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_relative_eq!(intersections[0].t, 6.80798, epsilon = 1e-5f64);
        assert_relative_eq!(intersections[1].t, 7.08872, epsilon = 1e-5f64);
    }

    // Scenario Outline: Normal vector on a cylinder
    //   Given cyl ← cylinder()
    //   When n ← local_normal_at(cyl, <point>)
    //   Then n = <normal>
    //
    //   Examples:
    //     | point           | normal           |
    //     | point(1, 0, 0)  | vector(1, 0, 0)  |
    #[test]
    fn normal_vector_on_cylinder_a() {
        let cylinder = Cylinder::new(0);
        let normal = cylinder.local_normal_at(&Tuple::point(1.0, 0.0, 0.0));
        assert_eq!(normal, Tuple::vector(1.0, 0.0, 0.0));
    }

    // Scenario Outline: Normal vector on a cylinder
    //   Given cyl ← cylinder()
    //   When n ← local_normal_at(cyl, <point>)
    //   Then n = <normal>
    //
    //   Examples:
    //     | point           | normal           |
    //     | point(0, 5, -1) | vector(0, 0, -1) |
    #[test]
    fn normal_vector_on_cylinder_b() {
        let cylinder = Cylinder::new(0);
        let normal = cylinder.local_normal_at(&Tuple::point(0.0, 5.0, -1.0));
        assert_eq!(normal, Tuple::vector(0.0, 0.0, -1.0));
    }

    // Scenario Outline: Normal vector on a cylinder
    //   Given cyl ← cylinder()
    //   When n ← local_normal_at(cyl, <point>)
    //   Then n = <normal>
    //
    //   Examples:
    //     | point           | normal           |
    //     | point(0, -2, 1) | vector(0, 0, 1)  |
    #[test]
    fn normal_vector_on_cylinder_c() {
        let cylinder = Cylinder::new(0);
        let normal = cylinder.local_normal_at(&Tuple::point(0.0, -2.0, 1.0));
        assert_eq!(normal, Tuple::vector(0.0, 0.0, 1.0));
    }

    // Scenario Outline: Normal vector on a cylinder
    //   Given cyl ← cylinder()
    //   When n ← local_normal_at(cyl, <point>)
    //   Then n = <normal>
    //
    //   Examples:
    //     | point           | normal           |
    //     | point(-1, 1, 0) | vector(-1, 0, 0) |
    #[test]
    fn normal_vector_on_cylinder_d() {
        let cylinder = Cylinder::new(0);
        let normal = cylinder.local_normal_at(&Tuple::point(-1.0, 1.0, 0.0));
        assert_eq!(normal, Tuple::vector(-1.0, 0.0, 0.0));
    }

    // Scenario: The default minimum and maximum for a cylinder
    //   Given cyl ← cylinder()
    //   Then cyl.minimum = -infinity
    //     And cyl.maximum = infinity
    #[test]
    fn default_minimum_and_maximum_for_cylinder() {
        let cylinder = Cylinder::new(0);
        assert_eq!(cylinder.minimum, f64::NEG_INFINITY);
        assert_eq!(cylinder.maximum, f64::INFINITY);
    }

    // Scenario Outline: Intersecting a constrained cylinder
    //   Given cyl ← cylinder()
    //     And cyl.minimum ← 1
    //     And cyl.maximum ← 2
    //     And direction ← normalize(<direction>)
    //     And r ← ray(<point>, direction)
    //   When xs ← local_intersect(cyl, r)
    //   Then xs.count = <count>
    //
    //   Examples:
    //     |   | point             | direction         | count |
    //     | 1 | point(0, 1.5, 0)  | vector(0.1, 1, 0) | 0     |
    #[test]
    fn intersecting_constrained_cylinder_a() {
        let mut cylinder = Cylinder::new(0);
        cylinder.minimum = 1.0;
        cylinder.maximum = 2.0;
        let ray = Ray::new(Tuple::point(0.0, 1.5, 0.0), Tuple::vector(0.1, 1.0, 0.0));
        let intersections = cylinder.local_intersect(&ray);
        assert_eq!(intersections.len(), 0);
    }

    // Scenario Outline: Intersecting a constrained cylinder
    //   Given cyl ← cylinder()
    //     And cyl.minimum ← 1
    //     And cyl.maximum ← 2
    //     And direction ← normalize(<direction>)
    //     And r ← ray(<point>, direction)
    //   When xs ← local_intersect(cyl, r)
    //   Then xs.count = <count>
    //
    //   Examples:
    //     |   | point             | direction         | count |
    //     | 2 | point(0, 3, -5)   | vector(0, 0, 1)   | 0     |
    #[test]
    fn intersecting_constrained_cylinder_b() {
        let mut cylinder = Cylinder::new(0);
        cylinder.minimum = 1.0;
        cylinder.maximum = 2.0;
        let ray = Ray::new(Tuple::point(0.0, 3.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let intersections = cylinder.local_intersect(&ray);
        assert_eq!(intersections.len(), 0);
    }

    // Scenario Outline: Intersecting a constrained cylinder
    //   Given cyl ← cylinder()
    //     And cyl.minimum ← 1
    //     And cyl.maximum ← 2
    //     And direction ← normalize(<direction>)
    //     And r ← ray(<point>, direction)
    //   When xs ← local_intersect(cyl, r)
    //   Then xs.count = <count>
    //
    //   Examples:
    //     |   | point             | direction         | count |
    //     | 3 | point(0, 0, -5)   | vector(0, 0, 1)   | 0     |
    #[test]
    fn intersecting_constrained_cylinder_c() {
        let mut cylinder = Cylinder::new(0);
        cylinder.minimum = 1.0;
        cylinder.maximum = 2.0;
        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let intersections = cylinder.local_intersect(&ray);
        assert_eq!(intersections.len(), 0);
    }

    // Scenario Outline: Intersecting a constrained cylinder
    //   Given cyl ← cylinder()
    //     And cyl.minimum ← 1
    //     And cyl.maximum ← 2
    //     And direction ← normalize(<direction>)
    //     And r ← ray(<point>, direction)
    //   When xs ← local_intersect(cyl, r)
    //   Then xs.count = <count>
    //
    //   Examples:
    //     |   | point             | direction         | count |
    //     | 4 | point(0, 2, -5)   | vector(0, 0, 1)   | 0     |
    #[test]
    fn intersecting_constrained_cylinder_d() {
        let mut cylinder = Cylinder::new(0);
        cylinder.minimum = 1.0;
        cylinder.maximum = 2.0;
        let ray = Ray::new(Tuple::point(0.0, 2.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let intersections = cylinder.local_intersect(&ray);
        assert_eq!(intersections.len(), 0);
    }

    // Scenario Outline: Intersecting a constrained cylinder
    //   Given cyl ← cylinder()
    //     And cyl.minimum ← 1
    //     And cyl.maximum ← 2
    //     And direction ← normalize(<direction>)
    //     And r ← ray(<point>, direction)
    //   When xs ← local_intersect(cyl, r)
    //   Then xs.count = <count>
    //
    //   Examples:
    //     |   | point             | direction         | count |
    //     | 5 | point(0, 1, -5)   | vector(0, 0, 1)   | 0     |
    #[test]
    fn intersecting_constrained_cylinder_e() {
        let mut cylinder = Cylinder::new(0);
        cylinder.minimum = 1.0;
        cylinder.maximum = 2.0;
        let ray = Ray::new(Tuple::point(0.0, 1.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let intersections = cylinder.local_intersect(&ray);
        assert_eq!(intersections.len(), 0);
    }

    // Scenario Outline: Intersecting a constrained cylinder
    //   Given cyl ← cylinder()
    //     And cyl.minimum ← 1
    //     And cyl.maximum ← 2
    //     And direction ← normalize(<direction>)
    //     And r ← ray(<point>, direction)
    //   When xs ← local_intersect(cyl, r)
    //   Then xs.count = <count>
    //
    //   Examples:
    //     |   | point             | direction         | count |
    //     | 6 | point(0, 1.5, -2) | vector(0, 0, 1)   | 2     |
    #[test]
    fn intersecting_constrained_cylinder_f() {
        let mut cylinder = Cylinder::new(0);
        cylinder.minimum = 1.0;
        cylinder.maximum = 2.0;
        let ray = Ray::new(Tuple::point(0.0, 1.5, -2.0), Tuple::vector(0.0, 0.0, 1.0));
        let intersections = cylinder.local_intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, 1.0);
        assert_eq!(intersections[1].t, 3.0);
    }

    // Scenario: The default closed value for a cylinder
    //   Given cyl ← cylinder()
    //   Then cyl.closed = false
    #[test]
    fn default_closed_value_for_cylinder() {
        let cylinder = Cylinder::new(0);
        assert!(!cylinder.closed);
    }

    // Scenario Outline: Intersecting the caps of a closed cylinder
    //   Given cyl ← cylinder()
    //     And cyl.minimum ← 1
    //     And cyl.maximum ← 2
    //     And cyl.closed ← true
    //     And direction ← normalize(<direction>)
    //     And r ← ray(<point>, direction)
    //   When xs ← local_intersect(cyl, r)
    //   Then xs.count = <count>
    //
    //   Examples:
    //     |   | point            | direction        | count |
    //     | 1 | point(0, 3, 0)   | vector(0, -1, 0) | 2     |
    #[test]
    fn intersecting_caps_of_closed_cylinder_a() {
        let mut cylinder = Cylinder::new(0);
        cylinder.minimum = 1.0;
        cylinder.maximum = 2.0;
        cylinder.closed = true;
        let ray = Ray::new(Tuple::point(0.0, 3.0, 0.0), Tuple::vector(0.0, -1.0, 0.0));
        let intersections = cylinder.local_intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, 2.0);
        assert_eq!(intersections[0].object.id(), cylinder.id());
        assert_eq!(intersections[1].t, 1.0);
        assert_eq!(intersections[1].object.id(), cylinder.id());
    }

    // Scenario Outline: Intersecting the caps of a closed cylinder
    //   Given cyl ← cylinder()
    //     And cyl.minimum ← 1
    //     And cyl.maximum ← 2
    //     And cyl.closed ← true
    //     And direction ← normalize(<direction>)
    //     And r ← ray(<point>, direction)
    //   When xs ← local_intersect(cyl, r)
    //   Then xs.count = <count>
    //
    //   Examples:
    //     |   | point            | direction        | count |
    //     | 2 | point(0, 3, -2)  | vector(0, -1, 2) | 2     |
    #[test]
    fn intersecting_caps_of_closed_cylinder_b() {
        let mut cylinder = Cylinder::new(0);
        cylinder.minimum = 1.0;
        cylinder.maximum = 2.0;
        cylinder.closed = true;
        let ray = Ray::new(Tuple::point(0.0, 3.0, -2.0), Tuple::vector(0.0, -1.0, 2.0));
        let intersections = cylinder.local_intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, 1.5);
        assert_eq!(intersections[0].object.id(), cylinder.id());
        assert_eq!(intersections[1].t, 1.0);
        assert_eq!(intersections[1].object.id(), cylinder.id());
    }

    // Scenario Outline: Intersecting the caps of a closed cylinder
    //   Given cyl ← cylinder()
    //     And cyl.minimum ← 1
    //     And cyl.maximum ← 2
    //     And cyl.closed ← true
    //     And direction ← normalize(<direction>)
    //     And r ← ray(<point>, direction)
    //   When xs ← local_intersect(cyl, r)
    //   Then xs.count = <count>
    //
    //   Examples:
    //     |   | point            | direction        | count |
    //     | 3 | point(0, 4, -2)  | vector(0, -1, 1) | 2     | # corner case
    #[test]
    fn intersecting_caps_of_closed_cylinder_c() {
        let mut cylinder = Cylinder::new(0);
        cylinder.minimum = 1.0;
        cylinder.maximum = 2.0;
        cylinder.closed = true;
        let ray = Ray::new(Tuple::point(0.0, 4.0, -2.0), Tuple::vector(0.0, -1.0, 1.0));
        let intersections = cylinder.local_intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, 3.0);
        assert_eq!(intersections[0].object.id(), cylinder.id());
        assert_eq!(intersections[1].t, 2.0);
        assert_eq!(intersections[1].object.id(), cylinder.id());
    }

    // Scenario Outline: Intersecting the caps of a closed cylinder
    //   Given cyl ← cylinder()
    //     And cyl.minimum ← 1
    //     And cyl.maximum ← 2
    //     And cyl.closed ← true
    //     And direction ← normalize(<direction>)
    //     And r ← ray(<point>, direction)
    //   When xs ← local_intersect(cyl, r)
    //   Then xs.count = <count>
    //
    //   Examples:
    //     |   | point            | direction        | count |
    //     | 4 | point(0, 0, -2)  | vector(0, 1, 2)  | 2     |
    #[test]
    fn intersecting_caps_of_closed_cylinder_d() {
        let mut cylinder = Cylinder::new(0);
        cylinder.minimum = 1.0;
        cylinder.maximum = 2.0;
        cylinder.closed = true;
        let ray = Ray::new(Tuple::point(0.0, 0.0, -2.0), Tuple::vector(0.0, 1.0, 2.0));
        let intersections = cylinder.local_intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, 1.5);
        assert_eq!(intersections[0].object.id(), cylinder.id());
        assert_eq!(intersections[1].t, 1.0);
        assert_eq!(intersections[1].object.id(), cylinder.id());
    }

    // Scenario Outline: Intersecting the caps of a closed cylinder
    //   Given cyl ← cylinder()
    //     And cyl.minimum ← 1
    //     And cyl.maximum ← 2
    //     And cyl.closed ← true
    //     And direction ← normalize(<direction>)
    //     And r ← ray(<point>, direction)
    //   When xs ← local_intersect(cyl, r)
    //   Then xs.count = <count>
    //
    //   Examples:
    //     |   | point            | direction        | count |
    //     | 5 | point(0, -1, -2) | vector(0, 1, 1)  | 2     | # corner case
    #[test]
    fn intersecting_caps_of_closed_cylinder_e() {
        let mut cylinder = Cylinder::new(0);
        cylinder.minimum = 1.0;
        cylinder.maximum = 2.0;
        cylinder.closed = true;
        let ray = Ray::new(Tuple::point(0.0, -1.0, -2.0), Tuple::vector(0.0, 1.0, 1.0));
        let intersections = cylinder.local_intersect(&ray);
        assert_eq!(intersections.len(), 2);
        assert_eq!(intersections[0].t, 2.0);
        assert_eq!(intersections[0].object.id(), cylinder.id());
        assert_eq!(intersections[1].t, 3.0);
        assert_eq!(intersections[1].object.id(), cylinder.id());
    }

    // Scenario Outline: The normal vector on a cylinder's end caps
    //   Given cyl ← cylinder()
    //     And cyl.minimum ← 1
    //     And cyl.maximum ← 2
    //     And cyl.closed ← true
    //   When n ← local_normal_at(cyl, <point>)
    //   Then n = <normal>
    //
    //   Examples:
    //     | point            | normal           |
    //     | point(0, 1, 0)   | vector(0, -1, 0) |
    #[test]
    fn normal_vector_on_cylinder_end_caps_a() {
        let mut cylinder = Cylinder::new(0);
        cylinder.minimum = 1.0;
        cylinder.maximum = 2.0;
        cylinder.closed = true;
        let normal = cylinder.local_normal_at(&Tuple::point(0.0, 1.0, 0.0));
        assert_eq!(normal, Tuple::vector(0.0, -1.0, 0.0));
    }

    // Scenario Outline: The normal vector on a cylinder's end caps
    //   Given cyl ← cylinder()
    //     And cyl.minimum ← 1
    //     And cyl.maximum ← 2
    //     And cyl.closed ← true
    //   When n ← local_normal_at(cyl, <point>)
    //   Then n = <normal>
    //
    //   Examples:
    //     | point            | normal           |
    //     | point(0.5, 1, 0) | vector(0, -1, 0) |
    #[test]
    fn normal_vector_on_cylinder_end_caps_b() {
        let mut cylinder = Cylinder::new(0);
        cylinder.minimum = 1.0;
        cylinder.maximum = 2.0;
        cylinder.closed = true;
        let normal = cylinder.local_normal_at(&Tuple::point(0.5, 1.0, 0.0));
        assert_eq!(normal, Tuple::vector(0.0, -1.0, 0.0));
    }

    // Scenario Outline: The normal vector on a cylinder's end caps
    //   Given cyl ← cylinder()
    //     And cyl.minimum ← 1
    //     And cyl.maximum ← 2
    //     And cyl.closed ← true
    //   When n ← local_normal_at(cyl, <point>)
    //   Then n = <normal>
    //
    //   Examples:
    //     | point            | normal           |
    //     | point(0, 1, 0.5) | vector(0, -1, 0) |
    #[test]
    fn normal_vector_on_cylinder_end_caps_c() {
        let mut cylinder = Cylinder::new(0);
        cylinder.minimum = 1.0;
        cylinder.maximum = 2.0;
        cylinder.closed = true;
        let normal = cylinder.local_normal_at(&Tuple::point(0.0, 1.0, 0.5));
        assert_eq!(normal, Tuple::vector(0.0, -1.0, 0.0));
    }

    // Scenario Outline: The normal vector on a cylinder's end caps
    //   Given cyl ← cylinder()
    //     And cyl.minimum ← 1
    //     And cyl.maximum ← 2
    //     And cyl.closed ← true
    //   When n ← local_normal_at(cyl, <point>)
    //   Then n = <normal>
    //
    //   Examples:
    //     | point            | normal           |
    //     | point(0, 2, 0)   | vector(0, 1, 0)  |
    #[test]
    fn normal_vector_on_cylinder_end_caps_d() {
        let mut cylinder = Cylinder::new(0);
        cylinder.minimum = 1.0;
        cylinder.maximum = 2.0;
        cylinder.closed = true;
        let normal = cylinder.local_normal_at(&Tuple::point(0.0, 2.0, 0.0));
        assert_eq!(normal, Tuple::vector(0.0, 1.0, 0.0));
    }

    // Scenario Outline: The normal vector on a cylinder's end caps
    //   Given cyl ← cylinder()
    //     And cyl.minimum ← 1
    //     And cyl.maximum ← 2
    //     And cyl.closed ← true
    //   When n ← local_normal_at(cyl, <point>)
    //   Then n = <normal>
    //
    //   Examples:
    //     | point            | normal           |
    //     | point(0.5, 2, 0) | vector(0, 1, 0)  |
    #[test]
    fn normal_vector_on_cylinder_end_caps_e() {
        let mut cylinder = Cylinder::new(0);
        cylinder.minimum = 1.0;
        cylinder.maximum = 2.0;
        cylinder.closed = true;
        let normal = cylinder.local_normal_at(&Tuple::point(0.5, 2.0, 0.0));
        assert_eq!(normal, Tuple::vector(0.0, 1.0, 0.0));
    }

    // Scenario Outline: The normal vector on a cylinder's end caps
    //   Given cyl ← cylinder()
    //     And cyl.minimum ← 1
    //     And cyl.maximum ← 2
    //     And cyl.closed ← true
    //   When n ← local_normal_at(cyl, <point>)
    //   Then n = <normal>
    //
    //   Examples:
    //     | point            | normal           |
    //     | point(0, 2, 0.5) | vector(0, 1, 0)  |
    #[test]
    fn normal_vector_on_cylinder_end_caps_f() {
        let mut cylinder = Cylinder::new(0);
        cylinder.minimum = 1.0;
        cylinder.maximum = 2.0;
        cylinder.closed = true;
        let normal = cylinder.local_normal_at(&Tuple::point(0.0, 2.0, 0.5));
        assert_eq!(normal, Tuple::vector(0.0, 1.0, 0.0));
    }
}
