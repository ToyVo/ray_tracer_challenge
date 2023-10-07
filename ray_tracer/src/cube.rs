use crate::{Intersection, Material, Ray, Shape, Transform};
use nalgebra_glm::{vec4, DMat4, DVec4};

#[derive(PartialEq, Debug, Clone)]
pub struct Cube {
    transform: DMat4,
    material: Material,
    id: u32,
}

impl Cube {
    pub fn new(id: u32) -> Cube {
        Cube {
            transform: DMat4::identity(),
            material: Material::default(),
            id,
        }
    }

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
    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let (x_min, x_max) = Self::check_axis(ray.origin.x, ray.direction.x);
        let (y_min, y_max) = Self::check_axis(ray.origin.y, ray.direction.y);
        let (z_min, z_max) = Self::check_axis(ray.origin.z, ray.direction.z);
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
    fn local_normal_at(&self, point: &DVec4) -> DVec4 {
        let max = point.x.abs().max(point.y.abs()).max(point.z.abs());
        if max == point.x.abs() {
            vec4(point.x, 0.0, 0.0, 0.)
        } else if max == point.y.abs() {
            vec4(0.0, point.y, 0.0, 0.)
        } else {
            vec4(0.0, 0.0, point.z, 0.)
        }
    }
    fn id(&self) -> u32 {
        self.id
    }
}

impl Transform for Cube {
    fn transform(&self) -> &DMat4 {
        &self.transform
    }
    fn transform_mut(&mut self) -> &mut DMat4 {
        &mut self.transform
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra_glm::vec4;

    #[test]
    fn ray_intersects_cube_pos_x() {
        let cube = Cube::new(0);
        let ray = Ray::new(vec4(5.0, 0.5, 0.0, 1.), vec4(-1.0, 0.0, 0.0, 0.));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, 4.0);
        assert_eq!(intersection[1].t, 6.0);
    }

    #[test]
    fn ray_intersects_cube_neg_x() {
        let cube = Cube::new(0);
        let ray = Ray::new(vec4(-5.0, 0.5, 0.0, 1.), vec4(1.0, 0.0, 0.0, 0.));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, 4.0);
        assert_eq!(intersection[1].t, 6.0);
    }

    #[test]
    fn ray_intersects_cube_pos_y() {
        let cube = Cube::new(0);
        let ray = Ray::new(vec4(0.5, 5.0, 0.0, 1.), vec4(0.0, -1.0, 0.0, 0.));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, 4.0);
        assert_eq!(intersection[1].t, 6.0);
    }

    #[test]
    fn ray_intersects_cube_neg_y() {
        let cube = Cube::new(0);
        let ray = Ray::new(vec4(0.5, -5.0, 0.0, 1.), vec4(0.0, 1.0, 0.0, 0.));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, 4.0);
        assert_eq!(intersection[1].t, 6.0);
    }

    #[test]
    fn ray_intersects_cube_pos_z() {
        let cube = Cube::new(0);
        let ray = Ray::new(vec4(0.5, 0.0, 5.0, 1.), vec4(0.0, 0.0, -1.0, 0.));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, 4.0);
        assert_eq!(intersection[1].t, 6.0);
    }

    #[test]
    fn ray_intersects_cube_neg_z() {
        let cube = Cube::new(0);
        let ray = Ray::new(vec4(0.5, 0.0, -5.0, 1.), vec4(0.0, 0.0, 1.0, 0.));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, 4.0);
        assert_eq!(intersection[1].t, 6.0);
    }

    #[test]
    fn ray_intersects_cube_inside() {
        let cube = Cube::new(0);
        let ray = Ray::new(vec4(0.0, 0.5, 0.0, 1.), vec4(0.0, 0.0, 1.0, 0.));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 2);
        assert_eq!(intersection[0].t, -1.0);
        assert_eq!(intersection[1].t, 1.0);
    }

    #[test]
    fn ray_misses_cube_a() {
        let cube = Cube::new(0);
        let ray = Ray::new(vec4(-2.0, 0.0, 0.0, 1.), vec4(0.2673, 0.5345, 0.8018, 0.));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 0);
    }

    #[test]
    fn ray_misses_cube_b() {
        let cube = Cube::new(0);
        let ray = Ray::new(vec4(0.0, -2.0, 0.0, 1.), vec4(0.8018, 0.2673, 0.5345, 0.));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 0);
    }

    #[test]
    fn ray_misses_cube_c() {
        let cube = Cube::new(0);
        let ray = Ray::new(vec4(0.0, 0.0, -2.0, 1.), vec4(0.5345, 0.8018, 0.2673, 0.));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 0);
    }

    #[test]
    fn ray_misses_cube_d() {
        let cube = Cube::new(0);
        let ray = Ray::new(vec4(2.0, 0.0, 2.0, 1.), vec4(0.0, 0.0, -1.0, 0.));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 0);
    }

    #[test]
    fn ray_misses_cube_e() {
        let cube = Cube::new(0);
        let ray = Ray::new(vec4(0.0, 2.0, 2.0, 1.), vec4(0.0, -1.0, 0.0, 0.));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 0);
    }

    #[test]
    fn ray_misses_cube_f() {
        let cube = Cube::new(0);
        let ray = Ray::new(vec4(2.0, 2.0, 0.0, 1.), vec4(-1.0, 0.0, 0.0, 0.));
        let intersection = cube.local_intersect(&ray);
        assert_eq!(intersection.len(), 0);
    }

    #[test]
    fn normal_on_surface_of_cube_pos_x() {
        let cube = Cube::new(0);
        let point = vec4(1.0, 0.5, -0.8, 1.);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, vec4(1.0, 0.0, 0.0, 0.));
    }

    #[test]
    fn normal_on_surface_of_cube_neg_x() {
        let cube = Cube::new(0);
        let point = vec4(-1.0, -0.2, 0.9, 1.);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, vec4(-1.0, 0.0, 0.0, 0.));
    }

    #[test]
    fn normal_on_surface_of_cube_pos_y() {
        let cube = Cube::new(0);
        let point = vec4(-0.4, 1.0, -0.1, 1.);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, vec4(0.0, 1.0, 0.0, 0.));
    }

    #[test]
    fn normal_on_surface_of_cube_neg_y() {
        let cube = Cube::new(0);
        let point = vec4(0.3, -1.0, -0.7, 1.);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, vec4(0.0, -1.0, 0.0, 0.));
    }

    #[test]
    fn normal_on_surface_of_cube_pos_z() {
        let cube = Cube::new(0);
        let point = vec4(-0.6, 0.3, 1.0, 1.);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, vec4(0.0, 0.0, 1.0, 0.));
    }

    #[test]
    fn normal_on_surface_of_cube_neg_z() {
        let cube = Cube::new(0);
        let point = vec4(0.4, 0.4, -1.0, 1.);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, vec4(0.0, 0.0, -1.0, 0.));
    }

    #[test]
    fn normal_on_corner_of_cube_pos() {
        let cube = Cube::new(0);
        let point = vec4(1.0, 1.0, 1.0, 1.);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, vec4(1.0, 0.0, 0.0, 0.));
    }

    #[test]
    fn normal_on_corner_of_cube_neg() {
        let cube = Cube::new(0);
        let point = vec4(-1.0, -1.0, -1.0, 1.);
        let normal = cube.local_normal_at(&point);
        assert_eq!(normal, vec4(-1.0, 0.0, 0.0, 0.));
    }
}
