use crate::{Intersection, Material, Matrix, Ray, Shape, Transform, Tuple};

#[derive(PartialEq, Debug, Clone)]
pub struct Cube {
    transform: Matrix,
    material: Material,
    id: u32,
}

impl Cube {
    pub fn new(id: u32) -> Cube {
        Cube {
            transform: Matrix::identity(4),
            material: Material::default(),
            id,
        }
    }

    // function check_axis(origin, direction)
    //   tmin_numerator = (-1 - origin)
    //   tmax_numerator = (1 - origin)
    //
    //   if abs(direction) >= EPSILON
    //     tmin ← tmin_numerator / direction
    //     tmax ← tmax_numerator / direction
    //   else
    //     tmin ← tmin_numerator * INFINITY
    //     tmax ← tmax_numerator * INFINITY
    //   end if
    //
    //   if tmin > tmax then swap(tmin, tmax)
    //
    //   return tmin, tmax
    // end function
    fn check_axis(origin: f64, direction: f64) -> (f64, f64) {
        let t_min_numerator = -1.0 - origin;
        let t_max_numerator = 1.0 - origin;

        let (t_min, t_max) = if direction.abs() >= f64::EPSILON {
            (t_min_numerator / direction, t_max_numerator / direction)
        } else {
            (
                t_min_numerator * f64::INFINITY,
                t_max_numerator * f64::INFINITY,
            )
        };

        if t_min > t_max {
            (t_max, t_min)
        } else {
            (t_min, t_max)
        }
    }
}

impl Shape for Cube {
    fn material(&self) -> &Material {
        &self.material
    }
    fn material_mut(&mut self) -> &mut Material {
        &mut self.material
    }

    // function local_intersect(cube, ray)
    //   xtmin, xtmax ← check_axis(ray.origin.x, ray.direction.x)
    //   ytmin, ytmax ← check_axis(ray.origin.y, ray.direction.y)
    //   ztmin, ztmax ← check_axis(ray.origin.z, ray.direction.z)
    //
    //   tmin ← max(xtmin, ytmin, ztmin)
    //   tmax ← min(xtmax, ytmax, ztmax)
    //
    //   return () if tmin > tmax
    //
    //   return ( intersection(tmin, cube), intersection(tmax, cube) )
    // end function
    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let (x_min, x_max) = Self::check_axis(ray.origin.x(), ray.direction.x());
        let (y_min, y_max) = Self::check_axis(ray.origin.y(), ray.direction.y());
        let (z_min, z_max) = Self::check_axis(ray.origin.z(), ray.direction.z());
        let t_min = x_min.max(y_min).max(z_min);
        let t_max = x_max.min(y_max).min(z_max);
        if t_min > t_max {
            return vec![];
        }
        let shape = Box::new(self.clone());
        vec![
            Intersection::new(t_min, shape.clone()),
            Intersection::new(t_max, shape),
        ]
    }

    // function local_normal_at(cube, point)
    //   maxc ← max(abs(point.x), abs(point.y), abs(point.z))
    //
    //   if maxc = abs(point.x) then
    //     return vector(point.x, 0, 0)
    //   else if maxc = abs(point.y) then
    //     return vector(0, point.y, 0)
    //   end if
    //
    //   return vector(0, 0, point.z)
    // end function
    fn local_normal_at(&self, point: &Tuple) -> Tuple {
        let max = point.x().abs().max(point.y().abs()).max(point.z().abs());
        if max == point.x().abs() {
            Tuple::vector(point.x(), 0.0, 0.0)
        } else if max == point.y().abs() {
            Tuple::vector(0.0, point.y(), 0.0)
        } else {
            Tuple::vector(0.0, 0.0, point.z())
        }
    }
    fn id(&self) -> u32 {
        self.id
    }
}

impl Transform for Cube {
    fn transform(&self) -> &Matrix {
        &self.transform
    }
    fn transform_mut(&mut self) -> &mut Matrix {
        &mut self.transform
    }
}

// Feature: Cubes
#[cfg(test)]
mod tests {
    use super::*;

    // Scenario Outline: A ray intersects a cube
    //   Given c ← cube()
    //     And r ← ray(<origin>, <direction>)
    //   When xs ← local_intersect(c, r)
    //   Then xs.count = 2
    //     And xs[0].t = <t1>
    //     And xs[1].t = <t2>
    //
    //   Examples:
    //     |        | origin            | direction        | t1 | t2 |
    //     | +x     | point(5, 0.5, 0)  | vector(-1, 0, 0) |  4 |  6 |
    #[test]
    fn ray_intersects_cube_pos_x() {
        let cube = Cube::new(0);
        let ray = Ray::new(Tuple::point(5.0, 0.5, 0.0), Tuple::vector(-1.0, 0.0, 0.0));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, 4.0);
        assert_eq!(intersection[1].t, 6.0);
    }

    // Scenario Outline: A ray intersects a cube
    //   Given c ← cube()
    //     And r ← ray(<origin>, <direction>)
    //   When xs ← local_intersect(c, r)
    //   Then xs.count = 2
    //     And xs[0].t = <t1>
    //     And xs[1].t = <t2>
    //
    //   Examples:
    //     |        | origin            | direction        | t1 | t2 |
    //     | -x     | point(-5, 0.5, 0) | vector(1, 0, 0)  |  4 |  6 |
    #[test]
    fn ray_intersects_cube_neg_x() {
        let cube = Cube::new(0);
        let ray = Ray::new(Tuple::point(-5.0, 0.5, 0.0), Tuple::vector(1.0, 0.0, 0.0));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, 4.0);
        assert_eq!(intersection[1].t, 6.0);
    }

    // Scenario Outline: A ray intersects a cube
    //   Given c ← cube()
    //     And r ← ray(<origin>, <direction>)
    //   When xs ← local_intersect(c, r)
    //   Then xs.count = 2
    //     And xs[0].t = <t1>
    //     And xs[1].t = <t2>
    //
    //   Examples:
    //     |        | origin            | direction        | t1 | t2 |
    //     | +y     | point(0.5, 5, 0)  | vector(0, -1, 0) |  4 |  6 |
    #[test]
    fn ray_intersects_cube_pos_y() {
        let cube = Cube::new(0);
        let ray = Ray::new(Tuple::point(0.5, 5.0, 0.0), Tuple::vector(0.0, -1.0, 0.0));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, 4.0);
        assert_eq!(intersection[1].t, 6.0);
    }

    // Scenario Outline: A ray intersects a cube
    //   Given c ← cube()
    //     And r ← ray(<origin>, <direction>)
    //   When xs ← local_intersect(c, r)
    //   Then xs.count = 2
    //     And xs[0].t = <t1>
    //     And xs[1].t = <t2>
    //
    //   Examples:
    //     |        | origin            | direction        | t1 | t2 |
    //     | -y     | point(0.5, -5, 0) | vector(0, 1, 0)  |  4 |  6 |
    #[test]
    fn ray_intersects_cube_neg_y() {
        let cube = Cube::new(0);
        let ray = Ray::new(Tuple::point(0.5, -5.0, 0.0), Tuple::vector(0.0, 1.0, 0.0));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, 4.0);
        assert_eq!(intersection[1].t, 6.0);
    }

    // Scenario Outline: A ray intersects a cube
    //   Given c ← cube()
    //     And r ← ray(<origin>, <direction>)
    //   When xs ← local_intersect(c, r)
    //   Then xs.count = 2
    //     And xs[0].t = <t1>
    //     And xs[1].t = <t2>
    //
    //   Examples:
    //     |        | origin            | direction        | t1 | t2 |
    //     | +z     | point(0.5, 0, 5)  | vector(0, 0, -1) |  4 |  6 |
    #[test]
    fn ray_intersects_cube_pos_z() {
        let cube = Cube::new(0);
        let ray = Ray::new(Tuple::point(0.5, 0.0, 5.0), Tuple::vector(0.0, 0.0, -1.0));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, 4.0);
        assert_eq!(intersection[1].t, 6.0);
    }

    // Scenario Outline: A ray intersects a cube
    //   Given c ← cube()
    //     And r ← ray(<origin>, <direction>)
    //   When xs ← local_intersect(c, r)
    //   Then xs.count = 2
    //     And xs[0].t = <t1>
    //     And xs[1].t = <t2>
    //
    //   Examples:
    //     |        | origin            | direction        | t1 | t2 |
    //     | -z     | point(0.5, 0, -5) | vector(0, 0, 1)  |  4 |  6 |
    #[test]
    fn ray_intersects_cube_neg_z() {
        let cube = Cube::new(0);
        let ray = Ray::new(Tuple::point(0.5, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, 4.0);
        assert_eq!(intersection[1].t, 6.0);
    }

    // Scenario Outline: A ray intersects a cube
    //   Given c ← cube()
    //     And r ← ray(<origin>, <direction>)
    //   When xs ← local_intersect(c, r)
    //   Then xs.count = 2
    //     And xs[0].t = <t1>
    //     And xs[1].t = <t2>
    //
    //   Examples:
    //     |        | origin            | direction        | t1 | t2 |
    //     | inside | point(0, 0.5, 0)  | vector(0, 0, 1)  | -1 |  1 |
    #[test]
    fn ray_intersects_cube_inside() {
        let cube = Cube::new(0);
        let ray = Ray::new(Tuple::point(0.0, 0.5, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, -1.0);
        assert_eq!(intersection[1].t, 1.0);
    }

    // Scenario Outline: A ray misses a cube
    //   Given c ← cube()
    //     And r ← ray(<origin>, <direction>)
    //   When xs ← local_intersect(c, r)
    //   Then xs.count = 0
    //
    //   Examples:
    //     | origin           | direction                      |
    //     | point(-2, 0, 0)  | vector(0.2673, 0.5345, 0.8018) |
    #[test]
    fn ray_misses_cube_a() {
        let cube = Cube::new(0);
        let ray = Ray::new(
            Tuple::point(-2.0, 0.0, 0.0),
            Tuple::vector(0.2673, 0.5345, 0.8018),
        );
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 0);
    }

    // Scenario Outline: A ray misses a cube
    //   Given c ← cube()
    //     And r ← ray(<origin>, <direction>)
    //   When xs ← local_intersect(c, r)
    //   Then xs.count = 0
    //
    //   Examples:
    //     | origin           | direction                      |
    //     | point(0, -2, 0)  | vector(0.8018, 0.2673, 0.5345) |
    #[test]
    fn ray_misses_cube_b() {
        let cube = Cube::new(0);
        let ray = Ray::new(
            Tuple::point(0.0, -2.0, 0.0),
            Tuple::vector(0.8018, 0.2673, 0.5345),
        );
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 0);
    }

    // Scenario Outline: A ray misses a cube
    //   Given c ← cube()
    //     And r ← ray(<origin>, <direction>)
    //   When xs ← local_intersect(c, r)
    //   Then xs.count = 0
    //
    //   Examples:
    //     | origin           | direction                      |
    //     | point(0, 0, -2)  | vector(0.5345, 0.8018, 0.2673) |
    #[test]
    fn ray_misses_cube_c() {
        let cube = Cube::new(0);
        let ray = Ray::new(
            Tuple::point(0.0, 0.0, -2.0),
            Tuple::vector(0.5345, 0.8018, 0.2673),
        );
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 0);
    }

    // Scenario Outline: A ray misses a cube
    //   Given c ← cube()
    //     And r ← ray(<origin>, <direction>)
    //   When xs ← local_intersect(c, r)
    //   Then xs.count = 0
    //
    //   Examples:
    //     | origin           | direction                      |
    //     | point(2, 0, 2)   | vector(0, 0, -1)               |
    #[test]
    fn ray_misses_cube_d() {
        let cube = Cube::new(0);
        let ray = Ray::new(Tuple::point(2.0, 0.0, 2.0), Tuple::vector(0.0, 0.0, -1.0));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 0);
    }

    // Scenario Outline: A ray misses a cube
    //   Given c ← cube()
    //     And r ← ray(<origin>, <direction>)
    //   When xs ← local_intersect(c, r)
    //   Then xs.count = 0
    //
    //   Examples:
    //     | origin           | direction                      |
    //     | point(0, 2, 2)   | vector(0, -1, 0)               |
    #[test]
    fn ray_misses_cube_e() {
        let cube = Cube::new(0);
        let ray = Ray::new(Tuple::point(0.0, 2.0, 2.0), Tuple::vector(0.0, -1.0, 0.0));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 0);
    }

    // Scenario Outline: A ray misses a cube
    //   Given c ← cube()
    //     And r ← ray(<origin>, <direction>)
    //   When xs ← local_intersect(c, r)
    //   Then xs.count = 0
    //
    //   Examples:
    //     | origin           | direction                      |
    //     | point(2, 2, 0)   | vector(-1, 0, 0)               |
    #[test]
    fn ray_misses_cube_f() {
        let cube = Cube::new(0);
        let ray = Ray::new(Tuple::point(2.0, 2.0, 0.0), Tuple::vector(-1.0, 0.0, 0.0));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 0);
    }

    // Scenario Outline: The normal on the surface of a cube
    //   Given c ← cube()
    //     And p ← <point>
    //   When normal ← local_normal_at(c, p)
    //   Then normal = <normal>
    //
    //   Examples:
    //     | point                | normal           |
    //     | point(1, 0.5, -0.8)  | vector(1, 0, 0)  |
    #[test]
    fn normal_on_surface_of_cube_pos_x() {
        let cube = Cube::new(0);
        let point = Tuple::point(1.0, 0.5, -0.8);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, Tuple::vector(1.0, 0.0, 0.0));
    }

    // Scenario Outline: The normal on the surface of a cube
    //   Given c ← cube()
    //     And p ← <point>
    //   When normal ← local_normal_at(c, p)
    //   Then normal = <normal>
    //
    //   Examples:
    //     | point                | normal           |
    //     | point(-1, -0.2, 0.9) | vector(-1, 0, 0) |
    #[test]
    fn normal_on_surface_of_cube_neg_x() {
        let cube = Cube::new(0);
        let point = Tuple::point(-1.0, -0.2, 0.9);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, Tuple::vector(-1.0, 0.0, 0.0));
    }

    // Scenario Outline: The normal on the surface of a cube
    //   Given c ← cube()
    //     And p ← <point>
    //   When normal ← local_normal_at(c, p)
    //   Then normal = <normal>
    //
    //   Examples:
    //     | point                | normal           |
    //     | point(-0.4, 1, -0.1) | vector(0, 1, 0)  |
    #[test]
    fn normal_on_surface_of_cube_pos_y() {
        let cube = Cube::new(0);
        let point = Tuple::point(-0.4, 1.0, -0.1);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, Tuple::vector(0.0, 1.0, 0.0));
    }

    // Scenario Outline: The normal on the surface of a cube
    //   Given c ← cube()
    //     And p ← <point>
    //   When normal ← local_normal_at(c, p)
    //   Then normal = <normal>
    //
    //   Examples:
    //     | point                | normal           |
    //     | point(0.3, -1, -0.7) | vector(0, -1, 0) |
    #[test]
    fn normal_on_surface_of_cube_neg_y() {
        let cube = Cube::new(0);
        let point = Tuple::point(0.3, -1.0, -0.7);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, Tuple::vector(0.0, -1.0, 0.0));
    }

    // Scenario Outline: The normal on the surface of a cube
    //   Given c ← cube()
    //     And p ← <point>
    //   When normal ← local_normal_at(c, p)
    //   Then normal = <normal>
    //
    //   Examples:
    //     | point                | normal           |
    //     | point(-0.6, 0.3, 1)  | vector(0, 0, 1)  |
    #[test]
    fn normal_on_surface_of_cube_pos_z() {
        let cube = Cube::new(0);
        let point = Tuple::point(-0.6, 0.3, 1.0);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, Tuple::vector(0.0, 0.0, 1.0));
    }

    // Scenario Outline: The normal on the surface of a cube
    //   Given c ← cube()
    //     And p ← <point>
    //   When normal ← local_normal_at(c, p)
    //   Then normal = <normal>
    //
    //   Examples:
    //     | point                | normal           |
    //     | point(0.4, 0.4, -1)  | vector(0, 0, -1) |
    #[test]
    fn normal_on_surface_of_cube_neg_z() {
        let cube = Cube::new(0);
        let point = Tuple::point(0.4, 0.4, -1.0);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, Tuple::vector(0.0, 0.0, -1.0));
    }

    // Scenario Outline: The normal on the surface of a cube
    //   Given c ← cube()
    //     And p ← <point>
    //   When normal ← local_normal_at(c, p)
    //   Then normal = <normal>
    //
    //   Examples:
    //     | point                | normal           |
    //     | point(1, 1, 1)       | vector(1, 0, 0)  |
    #[test]
    fn normal_on_corner_of_cube_pos() {
        let cube = Cube::new(0);
        let point = Tuple::point(1.0, 1.0, 1.0);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, Tuple::vector(1.0, 0.0, 0.0));
    }

    // Scenario Outline: The normal on the surface of a cube
    //   Given c ← cube()
    //     And p ← <point>
    //   When normal ← local_normal_at(c, p)
    //   Then normal = <normal>
    //
    //   Examples:
    //     | point                | normal           |
    //     | point(-1, -1, -1)    | vector(-1, 0, 0) |
    #[test]
    fn normal_on_corner_of_cube_neg() {
        let cube = Cube::new(0);
        let point = Tuple::point(-1.0, -1.0, -1.0);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, Tuple::vector(-1.0, 0.0, 0.0));
    }
}
