use nalgebra_glm::{vec3, vec4, DVec4};
use ray_tracer::Canvas;

struct Projectile {
    position: DVec4,
    velocity: DVec4,
}

struct Environment {
    gravity: DVec4,
    wind: DVec4,
}

fn tick(env: &Environment, proj: &Projectile) -> Projectile {
    Projectile {
        position: proj.position + proj.velocity,
        velocity: proj.velocity + env.gravity + env.wind,
    }
}

fn main() {
    let mut proj = Projectile {
        position: vec4(0.0, 1.0, 0.0, 1.),
        velocity: vec4(1.0, 1.8, 0.0, 0.).normalize() * 11.25,
    };
    let env = Environment {
        gravity: vec4(0.0, -0.1, 0.0, 0.),
        wind: vec4(-0.01, 0.0, 0.0, 0.),
    };
    let width = 900;
    let height = 550;
    let mut ticks = 0;
    let mut canvas = Canvas::new(width, height, vec3(0.0, 0.0, 0.0));

    loop {
        proj = tick(&env, &proj);
        let x = proj.position.x.round() as usize;
        let y = height - proj.position.y.round() as usize;
        if y >= height || x >= width {
            break;
        }
        println!("x: {}, y: {}", x, y);
        canvas.write_pixel(x, y, vec3(1.0, 0.0, 0.0));
        ticks += 1;
    }

    println!("Took {} ticks", ticks);
    canvas.write_ppm("cannon.ppm").unwrap();
}
