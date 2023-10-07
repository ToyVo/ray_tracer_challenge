use nalgebra_glm::{scale, translate, vec3, vec4, DMat4};
use ray_tracer::Canvas;
use std::f64::consts::PI;

fn main() {
    let size = 512;
    let start = vec4(0., 1., 0., 1.);
    let transform = translate(
        &scale(
            &DMat4::identity(),
            &vec3(size as f64 * 3. / 8., size as f64 * 3. / 8., 1.),
        ),
        &vec3(size as f64 / 2., size as f64 / 2., 0.),
    );
    let mut canvas = Canvas::new(size, size, vec3(0.0, 0.0, 0.0));

    for hour in 0..12 {
        let point = nalgebra_glm::rotate_z_vec4(&start, (hour as f64 * PI) / 6.0);
        let result = transform * point;
        println!("{:?}", result);
        let x = result.x.round() as usize;
        let y = result.y.round() as usize;
        canvas.write_pixel(x, y, vec3(1.0, 0.0, 0.0));
    }

    canvas.write_ppm("clock.ppm").unwrap();
}
