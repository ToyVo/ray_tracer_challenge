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

    // function ray_for_pixel(camera, px, py)
    //   # the offset from the edge of the canvas to the pixel's center
    //   xoffset ← (px + 0.5) * camera.pixel_size
    //   yoffset ← (py + 0.5) * camera.pixel_size
    //
    //   # the untransformed coordinates of the pixel in world space.
    //   # (remember that the camera looks toward -z, so +x is to the *left*.)
    //   world_x ← camera.half_width - xoffset
    //   world_y ← camera.half_height - yoffset
    //
    //   # using the camera matrix, transform the canvas point and the origin,
    //   # and then compute the ray's direction vector.
    //   # (remember that the canvas is at z=-1)
    //   pixel ← inverse(camera.transform) * point(world_x, world_y, -1)
    //   origin ← inverse(camera.transform) * point(0, 0, 0)
    //   direction ← normalize(pixel - origin)
    //
    //   return ray(origin, direction)
    // end function
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

    // function render(camera, world)
    //   image ← canvas(camera.hsize, camera.vsize)
    //
    //   for y ← 0 to camera.vsize - 1
    //     for x ← 0 to camera.hsize - 1
    //       ray ← ray_for_pixel(camera, x, y)
    //       color ← color_at(world, ray)
    //       write_pixel(image, x, y, color)
    //     end for
    //   end for
    //
    //   return image
    // end function
    pub fn render(&self, world: &World) -> Canvas {
        let mut image = Canvas::new(self.h_size, self.v_size, Tuple::color(0., 0., 0.));
        for y in 0..self.v_size {
            for x in 0..self.h_size {
                let ray = self.ray_for_pixel(x, y);
                let color = world.color_at(&ray, 5);
                image.write_pixel(x, y, color);
            }
        }
        image
    }
}

// Feature: Camera
#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use std::f64::consts::{PI, SQRT_2};

    use crate::view_transform;

    use super::*;

    // Scenario: Constructing a camera
    //   Given hsize ← 160
    //     And vsize ← 120
    //     And field_of_view ← π/2
    //   When c ← camera(hsize, vsize, field_of_view)
    //   Then c.hsize = 160
    //     And c.vsize = 120
    //     And c.field_of_view = π/2
    //     And c.transform = identity_matrix
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

    // Scenario: The pixel size for a horizontal canvas
    //   Given c ← camera(200, 125, π/2)
    //   Then c.pixel_size = 0.01
    #[test]
    fn pixel_size_horizontal_canvas() {
        let camera = Camera::new(200, 125, PI / 2.);
        assert_relative_eq!(camera.pixel_size, 0.01)
    }

    // Scenario: The pixel size for a vertical canvas
    //   Given c ← camera(125, 200, π/2)
    //   Then c.pixel_size = 0.01
    #[test]
    fn pixel_size_vertical_canvas() {
        let camera = Camera::new(125, 200, PI / 2.);
        assert_relative_eq!(camera.pixel_size, 0.01)
    }

    // Scenario: Constructing a ray through the center of the canvas
    //   Given c ← camera(201, 101, π/2)
    //   When r ← ray_for_pixel(c, 100, 50)
    //   Then r.origin = point(0, 0, 0)
    //     And r.direction = vector(0, 0, -1)
    #[test]
    fn ray_through_center_of_canvas() {
        let camera = Camera::new(201, 101, PI / 2.);
        let ray = camera.ray_for_pixel(100, 50);
        assert_eq!(ray.origin, Tuple::point(0., 0., 0.));
        assert_relative_eq!(ray.direction, Tuple::vector(0., 0., -1.));
    }

    // Scenario: Constructing a ray through a corner of the canvas
    //   Given c ← camera(201, 101, π/2)
    //   When r ← ray_for_pixel(c, 0, 0)
    //   Then r.origin = point(0, 0, 0)
    //     And r.direction = vector(0.66519, 0.33259, -0.66851)
    #[test]
    fn ray_through_corner_of_canvas() {
        let camera = Camera::new(201, 101, PI / 2.);
        let ray = camera.ray_for_pixel(0, 0);
        assert_eq!(ray.origin, Tuple::point(0., 0., 0.));
        assert_relative_eq!(
            ray.direction,
            Tuple::vector(0.66519, 0.33259, -0.66851)
        );
    }

    // Scenario: Constructing a ray when the camera is transformed
    //   Given c ← camera(201, 101, π/2)
    //   When c.transform ← rotation_y(π/4) * translation(0, -2, 5)
    //     And r ← ray_for_pixel(c, 100, 50)
    //   Then r.origin = point(0, 2, -5)
    //     And r.direction = vector(√2/2, 0, -√2/2)
    #[test]
    fn ray_when_camera_is_transformed() {
        let mut camera = Camera::new(201, 101, PI / 2.);
        camera.transform = Matrix::rotation_y(PI / 4.) * Matrix::translation(0., -2., 5.);
        let ray = camera.ray_for_pixel(100, 50);
        assert_eq!(ray.origin, Tuple::point(0., 2., -5.));
        assert_relative_eq!(
            ray.direction,
            Tuple::vector(SQRT_2 / 2., 0., -SQRT_2 / 2.)
        );
    }

    // Scenario: Rendering a world with a camera
    //   Given w ← default_world()
    //     And c ← camera(11, 11, π/2)
    //     And from ← point(0, 0, -5)
    //     And to ← point(0, 0, 0)
    //     And up ← vector(0, 1, 0)
    //     And c.transform ← view_transform(from, to, up)
    //   When image ← render(c, w)
    //   Then pixel_at(image, 5, 5) = color(0.38066, 0.47583, 0.2855)
    #[test]
    fn rendering_world_with_camera() {
        let world = World::default();
        let mut camera = Camera::new(11, 11, PI / 2.);
        let from = Tuple::point(0., 0., -5.);
        let to = Tuple::point(0., 0., 0.);
        let up = Tuple::vector(0., 1., 0.);
        camera.transform = view_transform(from, to, up);
        let image = camera.render(&world);
        assert_relative_eq!(
            *image.pixel_at(5, 5),
            Tuple::color(0.38066, 0.47583, 0.2855),
            epsilon = 1e-5f64
        );
    }
}
