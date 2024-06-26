use nalgebra_glm::{vec3, vec4};
use ray_tracer::{Canvas, Intersection, Light, Ray, Shape, SolidPattern, Sphere};

fn main() {
    let size = 256;
    let ray_origin = vec4(0., 0., -5., 1.);
    let wall_z = 10.;
    let wall_size = 7.;
    let pixel_size = wall_size / size as f64;
    let mut sphere = Sphere::new(0);
    sphere.material_mut().pattern = Box::new(SolidPattern::new(vec3(1., 0.2, 1.)));
    let mut canvas = Canvas::new(size, size, vec3(0.0, 0.0, 0.0));
    let light = Light::new(vec4(-10., 10., -10., 1.), vec3(1., 1., 1.));

    for x in 0..size {
        let world_x = -wall_size / 2. + pixel_size * x as f64;
        for y in 0..size {
            let world_y = wall_size / 2. - pixel_size * y as f64;
            let target = vec4(world_x, world_y, wall_z, 1.);
            let ray = Ray::new(ray_origin, (target - ray_origin).normalize());
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
