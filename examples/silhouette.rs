use ray_tracer_challenge::{Canvas, Ray, Tuple, Sphere, Intersection};
use std::fs::File;
use std::io::Write;

fn main() {
    let size = 100;
    let ray_origin = Tuple::point(0., 0., -5.);
    let wall_z = 10.;
    let wall_size = 7.;
    let pixel_size = wall_size / size as f64;
    let sphere = Sphere::new();
    let mut canvas = Canvas::new(size, size, Tuple::color(0.0, 0.0, 0.0));

    for x in 0..size {
        let world_x = -wall_size / 2. + pixel_size * x as f64;
        for y in 0..size {
            let world_y = wall_size / 2. - pixel_size * y as f64;
            let target = Tuple::point(world_x, world_y, wall_z);
            let r = Ray::new(ray_origin.clone(), (target - &ray_origin).normalize());
            let xs = sphere.intersects(&r);
            if let Some(t) = Intersection::hit(&xs) {
                canvas.set(x,y,Tuple::color(1.,0.,0.));
            }
        }
    }

    // write the canvas to a file
    let mut file = File::create("silhouette.ppm").unwrap();
    file.write_all(canvas.to_string().as_bytes()).unwrap();
}
