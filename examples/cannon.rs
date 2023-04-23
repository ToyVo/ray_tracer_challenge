use ray_tracer_challenge::{Matrix, Tuple};
use std::fs::File;
use std::io::Write;

struct Projectile {
    position: Tuple<f64>,
    velocity: Tuple<f64>,
}

struct Environment {
    gravity: Tuple<f64>,
    wind: Tuple<f64>,
}

fn tick(env: &Environment, proj: &Projectile) -> Projectile {
    Projectile {
        position: proj.position.clone() + proj.velocity.clone(),
        velocity: proj.velocity.clone() + env.gravity.clone() + env.wind.clone(),
    }
}

fn main() {
    let mut proj = Projectile {
        position: Tuple::point(0.0, 1.0, 0.0),
        velocity: Tuple::vector(1.0, 1.8, 0.0).normalize() * 11.25,
    };
    let env = Environment {
        gravity: Tuple::vector(0.0, -0.1, 0.0),
        wind: Tuple::vector(-0.01, 0.0, 0.0),
    };

    let width = 900;
    let height = 550;
    let mut positions = Vec::new();
    loop {
        proj = tick(&env, &proj);
        let x = proj.position.x().round() as usize;
        let y = height - proj.position.y().round() as usize;
        if y >= height || x >= width {
            break;
        }
        println!("x: {}, y: {}", x, y);
        positions.push((x, y));
    }

    println!("Took {} ticks", positions.len());
    let mut canvas = Matrix::new(height, width, Tuple::color(0.0, 0.0, 0.0));

    positions.iter().for_each(|(x, y)| {
        canvas.set(*y, *x, Tuple::color(1.0, 0.0, 0.0));
    });

    // write the canvas to a file
    let mut file = File::create("cannon.ppm").unwrap();
    file.write_all(canvas.to_string().as_bytes()).unwrap();
}
