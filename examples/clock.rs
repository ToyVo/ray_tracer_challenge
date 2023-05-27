use ray_tracer_challenge::{Canvas, Matrix, Tuple};
use std::f64::consts::PI;
use std::fs::File;
use std::io::Write;

fn main() {
    let size = 512;
    let start = Tuple::point(0., 1., 0.);
    let transform = Matrix::scaling(size as f64 * 3. / 8., size as f64 * 3. / 8., 1.).translate(
        size as f64 / 2.,
        size as f64 / 2.,
        0.,
    );
    let mut canvas = Canvas::new(size, size, Tuple::color(0.0, 0.0, 0.0));

    for i in 0..12 {
        let point = start.rotate_z((i as f64 * PI) / 6.0);
        let result = &transform * point;
        println!("{:?}", result);
        let x = result.x().round() as usize;
        let y = result.y().round() as usize;
        canvas.set(y, x, Tuple::color(1.0, 0.0, 0.0));
    }

    // write the canvas to a file
    let mut file = File::create("clock.ppm").unwrap();
    file.write_all(canvas.to_string().as_bytes()).unwrap();
}
