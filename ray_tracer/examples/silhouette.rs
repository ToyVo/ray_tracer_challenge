use ray_tracer::{Canvas, Intersection, Light, Ray, Shape, SolidPattern, Sphere, Tuple};

// # start the ray at z = -5
// ray_origin ← point(0, 0, -5)
//
// # put the wall at z = 10
// wall_z ← 10
//
// wall_size ← 7.0
//
// canvas_pixels ← 100
//
// pixel_size ← wall_size / canvas_pixels
//
// half ← wall_size / 2
//
// canvas ← canvas(canvas_pixels, canvas_pixels)
//
// sphere.material ← material()
// sphere.material.color ← color(1, 0.2, 1)
//
// light_position ← point(-10, 10, -10)
// light_color    ← color(1, 1, 1)
// light          ← point_light(light_position, light_color)
//
// # for each row of pixels in the canvas
// for y ← 0 to canvas_pixels - 1
//
//   # compute the world y coordinate (top = +half, bottom = -half)
//   world_y ← half - pixel_size * y
//
//   # for each pixel in the row
//   for x ← 0 to canvas_pixels - 1
//
//     # compute the world x coordinate (left = -half, right = half)
//     world_x ← -half + pixel_size * x
//
//     # describe the point on the wall that the ray will target
//     position ← point(world_x, world_y, wall_z)
//
//     ray ← ray(ray_origin, normalize(position - ray_origin))
//     xs ← intersect(sphere, ray)
//
//     if hit(xs) is defined
//       point  ← position(ray, hit.t)
//       normal ← normal_at(hit.object, point)
//       eye    ← -ray.direction
//
//       color ← lighting(hit.object.material, light, point, eye, normal)
//
//       write_pixel(canvas, x, y, color)
//     end if
//
//   end for
//
// end for
fn main() {
    let size = 256;
    let ray_origin = Tuple::point(0., 0., -5.);
    let wall_z = 10.;
    let wall_size = 7.;
    let pixel_size = wall_size / size as f64;
    let mut sphere = Sphere::new(0);
    sphere.material_mut().pattern = Box::new(SolidPattern::new(Tuple::color(1., 0.2, 1.)));
    let mut canvas = Canvas::new(size, size, Tuple::color(0.0, 0.0, 0.0));
    let light = Light::new(Tuple::point(-10., 10., -10.), Tuple::color(1., 1., 1.));

    for x in 0..size {
        let world_x = -wall_size / 2. + pixel_size * x as f64;
        for y in 0..size {
            let world_y = wall_size / 2. - pixel_size * y as f64;
            let target = Tuple::point(world_x, world_y, wall_z);
            let ray = Ray::new(ray_origin.clone(), (target - &ray_origin).normalize());
            let intersections = sphere.intersect(&ray);
            if let Some(hit) = Intersection::hit(&intersections) {
                let point = ray.position(hit.t);
                let normal = hit.object.normal_at(&point);
                let eye = -ray.direction;
                let color = hit.object.material().lighting(
                    hit.object.as_ref(),
                    &light,
                    &point,
                    &eye,
                    &normal,
                    false,
                );
                canvas.write_pixel(x, y, color);
            }
        }
    }

    canvas.write_ppm("silhouette.ppm").unwrap();
}
