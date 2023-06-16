use crate::{Canvas, Matrix, Ray, Tuple, World};

pub struct Camera {
    h_size: usize,
    v_size: usize,
    pub transform: Matrix,
    pixel_size: f64,
    half_width: f64,
    half_height: f64,
}

impl Camera {
    pub fn new(h_size: usize, v_size: usize, field_of_view: f64) -> Camera {
        let half_view = (field_of_view / 2.).tan();
        let aspect = h_size as f64 / v_size as f64;
        let (half_width, half_height) = if aspect >= 1. {
            (half_view, half_view / aspect)
        } else {
            (half_view * aspect, half_view)
        };
        let pixel_size = (half_width * 2.) / h_size as f64;
        Camera {
            h_size,
            v_size,
            transform: Matrix::identity(4),
            pixel_size,
            half_width,
            half_height,
        }
    }

    pub fn ray_for_pixel(&self, px: usize, py: usize) -> Ray {
        let x_offset = (px as f64 + 0.5) * self.pixel_size;
        let y_offset = (py as f64 + 0.5) * self.pixel_size;
        let world_x = self.half_width - x_offset;
        let world_y = self.half_height - y_offset;
        let pixel = &self.transform.inverse() * Tuple::point(world_x, world_y, -1.);
        let origin = &self.transform.inverse() * Tuple::point(0., 0., 0.);
        let direction = (&pixel - &origin).normalize();
        Ray::new(origin, direction)
    }

    pub fn render(&self, world: &World) -> Canvas {
        let mut image = Canvas::new(self.h_size, self.v_size, Tuple::color(0., 0., 0.));
        for y in 0..self.v_size {
            for x in 0..self.h_size {
                let ray = self.ray_for_pixel(x, y);
                let color = world.color_at(&ray);
                image.write_pixel(x, y, color);
            }
        }
        image
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{PI, SQRT_2};

    use crate::view_transform;

    use super::*;

    #[test]
    fn constructing_camera() {
        let h_size = 160;
        let v_size = 120;
        let field_of_view = PI / 2.;
        let camera = Camera::new(h_size, v_size, field_of_view);
        assert_eq!(camera.h_size, 160);
        assert_eq!(camera.v_size, 120);
        assert_eq!(camera.transform, Matrix::identity(4));
    }

    #[test]
    fn pixel_size_horizontal_canvas() {
        let camera = Camera::new(200, 125, PI / 2.);
        assert!(camera.pixel_size - 0.01 < 1e-3f64);
    }

    #[test]
    fn pixel_size_vertical_canvas() {
        let camera = Camera::new(125, 200, PI / 2.);
        assert!(camera.pixel_size - 0.01 < 1e-3f64);
    }

    #[test]
    fn ray_through_center_of_canvas() {
        let camera = Camera::new(201, 101, PI / 2.);
        let ray = camera.ray_for_pixel(100, 50);
        assert_eq!(ray.origin, Tuple::point(0., 0., 0.));
        assert!(ray
            .direction
            .nearly_equals(&Tuple::vector(0., 0., -1.), 1e-3f64));
    }

    #[test]
    fn ray_through_corner_of_canvas() {
        let camera = Camera::new(201, 101, PI / 2.);
        let ray = camera.ray_for_pixel(0, 0);
        assert_eq!(ray.origin, Tuple::point(0., 0., 0.));
        assert!(ray
            .direction
            .nearly_equals(&Tuple::vector(0.66519, 0.33259, -0.66851), 1e-3f64));
    }

    #[test]
    fn ray_when_camera_is_transformed() {
        let mut camera = Camera::new(201, 101, PI / 2.);
        camera.transform = Matrix::rotation_y(PI / 4.) * Matrix::translation(0., -2., 5.);
        let ray = camera.ray_for_pixel(100, 50);
        assert_eq!(ray.origin, Tuple::point(0., 2., -5.));
        assert!(ray
            .direction
            .nearly_equals(&Tuple::vector(SQRT_2 / 2., 0., -SQRT_2 / 2.), 1e-3f64));
    }

    #[test]
    fn rendering_world_with_camera() {
        let world = World::default();
        let mut camera = Camera::new(11, 11, PI / 2.);
        let from = Tuple::point(0., 0., -5.);
        let to = Tuple::point(0., 0., 0.);
        let up = Tuple::vector(0., 1., 0.);
        camera.transform = view_transform(from, to, up);
        let image = camera.render(&world);
        assert!(image
            .pixel_at(5, 5)
            .nearly_equals(&Tuple::color(0.38066, 0.47583, 0.2855), 1e-3f64));
    }
}
