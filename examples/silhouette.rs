use ray_tracer_challenge::{Canvas, Ray, Tuple, Sphere, Intersection, Light};
use std::fs::File;
use std::io::Write;

fn main() {
    let size = 256;
    let ray_origin = Tuple::point(0., 0., -5.);
    let wall_z = 10.;
    let wall_size = 7.;
    let pixel_size = wall_size / size as f64;
    let mut sphere = Sphere::new();
    sphere.material.color = Tuple::color(1., 0.2, 1.);
    let mut canvas = Canvas::new(size, size, Tuple::color(0.0, 0.0, 0.0));
    let light = Light::new(Tuple::point(-10., 10., -10.), Tuple::color(1., 1., 1.));

    for x in 0..size {
        let world_x = -wall_size / 2. + pixel_size * x as f64;
        for y in 0..size {
            let world_y = wall_size / 2. - pixel_size * y as f64;
            let target = Tuple::point(world_x, world_y, wall_z);
            let r = Ray::new(ray_origin.clone(), (target - &ray_origin).normalize());
            let xs = sphere.intersects(&r);
            if let Some(hit) = Intersection::hit(&xs) {
                let point = r.position(hit.t);
                let normal = hit.object.normal_at(&point);
                let eye = -r.direction;
                let color = hit.object.material.lighting(&light, &point, &eye, &normal);
                canvas.set(x,y,color);
            }
        }
    }

    // write the canvas to a file
    let mut file = File::create("silhouette.ppm").unwrap();
    file.write_all(canvas.to_string().as_bytes()).unwrap();
}
