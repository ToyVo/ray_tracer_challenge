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
            material: Material::new(),
            id,
        }
    }

    fn check_axis(origin: f64, direction: f64) -> (f64, f64) {
        let tmin_numerator = -1.0 - origin;
        let tmax_numerator = 1.0 - origin;

        let (tmin, tmax) = if direction.abs() >= f64::EPSILON {
            (tmin_numerator / direction, tmax_numerator / direction)
        } else {
            (
                tmin_numerator * f64::INFINITY,
                tmax_numerator * f64::INFINITY,
            )
        };

        if tmin > tmax {
            (tmax, tmin)
        } else {
            (tmin, tmax)
        }
    }
}

impl Shape for Cube {
    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let (xtmin, xtmax) = Self::check_axis(ray.origin.x(), ray.direction.x());
        let (ytmin, ytmax) = Self::check_axis(ray.origin.y(), ray.direction.y());
        let (ztmin, ztmax) = Self::check_axis(ray.origin.z(), ray.direction.z());
        let tmin = xtmin.max(ytmin).max(ztmin);
        let tmax = xtmax.min(ytmax).min(ztmax);
        if tmin > tmax {
            return vec![];
        }
        let shape = Box::new(self.clone());
        vec![
            Intersection::new(tmin, shape.clone()),
            Intersection::new(tmax, shape),
        ]
    }
    fn local_normal_at(&self, point: &Tuple) -> Tuple {
        let maxc = point.x().abs().max(point.y().abs()).max(point.z().abs());
        if maxc == point.x().abs() {
            Tuple::vector(point.x(), 0.0, 0.0)
        } else if maxc == point.y().abs() {
            Tuple::vector(0.0, point.y(), 0.0)
        } else {
            Tuple::vector(0.0, 0.0, point.z())
        }
    }
    fn material(&self) -> &Material {
        &self.material
    }
    fn material_mut(&mut self) -> &mut Material {
        &mut self.material
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ray_intersects_cube_pos_x() {
        let cube = Cube::new(0);
        let ray = Ray::new(Tuple::point(5.0, 0.5, 0.0), Tuple::vector(-1.0, 0.0, 0.0));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, 4.0);
        assert_eq!(intersection[1].t, 6.0);
    }

    #[test]
    fn ray_intersects_cube_neg_x() {
        let cube = Cube::new(0);
        let ray = Ray::new(Tuple::point(-5.0, 0.5, 0.0), Tuple::vector(1.0, 0.0, 0.0));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, 4.0);
        assert_eq!(intersection[1].t, 6.0);
    }

    #[test]
    fn ray_intersects_cube_pos_y() {
        let cube = Cube::new(0);
        let ray = Ray::new(Tuple::point(0.5, 5.0, 0.0), Tuple::vector(0.0, -1.0, 0.0));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, 4.0);
        assert_eq!(intersection[1].t, 6.0);
    }

    #[test]
    fn ray_intersects_cube_neg_y() {
        let cube = Cube::new(0);
        let ray = Ray::new(Tuple::point(0.5, -5.0, 0.0), Tuple::vector(0.0, 1.0, 0.0));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, 4.0);
        assert_eq!(intersection[1].t, 6.0);
    }

    #[test]
    fn ray_intersects_cube_pos_z() {
        let cube = Cube::new(0);
        let ray = Ray::new(Tuple::point(0.5, 0.0, 5.0), Tuple::vector(0.0, 0.0, -1.0));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, 4.0);
        assert_eq!(intersection[1].t, 6.0);
    }

    #[test]
    fn ray_intersects_cube_neg_z() {
        let cube = Cube::new(0);
        let ray = Ray::new(Tuple::point(0.5, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, 4.0);
        assert_eq!(intersection[1].t, 6.0);
    }

    #[test]
    fn ray_intersects_cube_inside() {
        let cube = Cube::new(0);
        let ray = Ray::new(Tuple::point(0.0, 0.5, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, -1.0);
        assert_eq!(intersection[1].t, 1.0);
    }

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

    #[test]
    fn ray_misses_cube_d() {
        let cube = Cube::new(0);
        let ray = Ray::new(Tuple::point(2.0, 0.0, 2.0), Tuple::vector(0.0, 0.0, -1.0));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 0);
    }

    #[test]
    fn ray_misses_cube_e() {
        let cube = Cube::new(0);
        let ray = Ray::new(Tuple::point(0.0, 2.0, 2.0), Tuple::vector(0.0, -1.0, 0.0));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 0);
    }

    #[test]
    fn ray_misses_cube_f() {
        let cube = Cube::new(0);
        let ray = Ray::new(Tuple::point(2.0, 2.0, 0.0), Tuple::vector(-1.0, 0.0, 0.0));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 0);
    }

    #[test]
    fn normal_on_surface_of_cube_pos_x() {
        let cube = Cube::new(0);
        let point = Tuple::point(1.0, 0.5, -0.8);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, Tuple::vector(1.0, 0.0, 0.0));
    }

    #[test]
    fn normal_on_surface_of_cube_neg_x() {
        let cube = Cube::new(0);
        let point = Tuple::point(-1.0, -0.2, 0.9);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, Tuple::vector(-1.0, 0.0, 0.0));
    }

    #[test]
    fn normal_on_surface_of_cube_pos_y() {
        let cube = Cube::new(0);
        let point = Tuple::point(-0.4, 1.0, -0.1);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, Tuple::vector(0.0, 1.0, 0.0));
    }

    #[test]
    fn normal_on_surface_of_cube_neg_y() {
        let cube = Cube::new(0);
        let point = Tuple::point(0.3, -1.0, -0.7);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, Tuple::vector(0.0, -1.0, 0.0));
    }

    #[test]
    fn normal_on_surface_of_cube_pos_z() {
        let cube = Cube::new(0);
        let point = Tuple::point(-0.6, 0.3, 1.0);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, Tuple::vector(0.0, 0.0, 1.0));
    }

    #[test]
    fn normal_on_surface_of_cube_neg_z() {
        let cube = Cube::new(0);
        let point = Tuple::point(0.4, 0.4, -1.0);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, Tuple::vector(0.0, 0.0, -1.0));
    }

    #[test]
    fn normal_on_corner_of_cube_pos() {
        let cube = Cube::new(0);
        let point = Tuple::point(1.0, 1.0, 1.0);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, Tuple::vector(1.0, 0.0, 0.0));
    }

    #[test]
    fn normal_on_corner_of_cube_neg() {
        let cube = Cube::new(0);
        let point = Tuple::point(-1.0, -1.0, -1.0);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, Tuple::vector(-1.0, 0.0, 0.0));
    }
}
