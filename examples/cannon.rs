use ray_tracer_challenge::{Canvas, Tuple};
use std::fs::File;
use std::io::Write;

struct Projectile {
    position: Tuple,
    velocity: Tuple,
}

struct Environment {
    gravity: Tuple,
    wind: Tuple,
}

fn tick(env: &Environment, proj: &Projectile) -> Projectile {
    Projectile {
        position: &proj.position + &proj.velocity,
        velocity: &proj.velocity + &env.gravity + &env.wind,
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
    let mut ticks = 0;
    let mut canvas = Canvas::new(height, width, Tuple::color(0.0, 0.0, 0.0));

    loop {
        proj = tick(&env, &proj);
        let x = proj.position.x().round() as usize;
        let y = height - proj.position.y().round() as usize;
        if y >= height || x >= width {
            break;
        }
        println!("x: {}, y: {}", x, y);
        canvas.set(y, x, Tuple::color(1.0, 0.0, 0.0));
        ticks += 1;
    }

    println!("Took {} ticks", ticks);

    // write the canvas to a file
    let mut file = File::create("cannon.ppm").unwrap();
    file.write_all(canvas.to_string().as_bytes()).unwrap();
}
