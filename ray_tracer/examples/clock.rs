use std::f64::consts::PI;

use ray_tracer::{Canvas, Matrix, Tuple};

fn main() {
    let size = 512;
    let start = Tuple::point(0., 1., 0.);
    let transform = Matrix::scaling(size as f64 * 3. / 8., size as f64 * 3. / 8., 1.).translate(
        size as f64 / 2.,
        size as f64 / 2.,
        0.,
    );
    let mut canvas = Canvas::new(size, size, Tuple::color(0.0, 0.0, 0.0));

    for hour in 0..12 {
        let point = start.rotate_z((hour as f64 * PI) / 6.0);
        let result = &transform * point;
        println!("{:?}", result);
        let x = result.x().round() as usize;
        let y = result.y().round() as usize;
        canvas.write_pixel(x, y, Tuple::color(1.0, 0.0, 0.0));
    }

    canvas.write_ppm("clock.ppm").unwrap();
}
