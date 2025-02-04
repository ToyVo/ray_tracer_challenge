use ray_tracer::{Canvas, Tuple};

struct Projectile {
    position: Tuple,
    velocity: Tuple,
}

struct Environment {
    gravity: Tuple,
    wind: Tuple,
}

// function tick(env, proj)
//   position ← proj.position + proj.velocity
//   velocity ← proj.velocity + env.gravity + env.wind
//   return projectile(position, velocity)
// end function
fn tick(env: &Environment, proj: &Projectile) -> Projectile {
    Projectile {
        position: &proj.position + &proj.velocity,
        velocity: &proj.velocity + &env.gravity + &env.wind,
    }
}

fn main() {
    // # projectile starts one unit above the origin.
    // # velocity is normalized to 1 unit/tick.
    // p ← projectile(point(0, 1, 0), normalize(vector(1, 1, 0)))
    let mut proj = Projectile {
        position: Tuple::point(0.0, 1.0, 0.0),
        velocity: Tuple::vector(1.0, 1.8, 0.0).normalize() * 11.25,
    };

    // # gravity -0.1 unit/tick, and wind is -0.01 unit/tick.
    // e ← environment(vector(0, -0.1, 0), vector(-0.01, 0, 0))
    let env = Environment {
        gravity: Tuple::vector(0.0, -0.1, 0.0),
        wind: Tuple::vector(-0.01, 0.0, 0.0),
    };
    let width = 900;
    let height = 550;
    let mut ticks = 0;
    let mut canvas = Canvas::new(width, height, Tuple::color(0.0, 0.0, 0.0));

    loop {
        proj = tick(&env, &proj);
        let x = proj.position.x().round() as usize;
        let y = height - proj.position.y().round() as usize;
        if y >= height || x >= width {
            break;
        }
        println!("x: {}, y: {}", x, y);
        canvas.write_pixel(x, y, Tuple::color(1.0, 0.0, 0.0));
        ticks += 1;
    }

    println!("Took {} ticks", ticks);
    canvas.write_ppm("cannon.ppm").unwrap();
}
