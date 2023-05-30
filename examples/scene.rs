use ray_tracer_challenge::{Sphere, Matrix, Tuple, World, Light, Camera, view_transform};
use std::f64::consts::PI;

fn main() {
    let mut floor = Sphere::new();
    floor.transform = Matrix::scaling(10., 0.01, 10.);
    floor.material.color = Tuple::color(1., 0.9, 0.9);
    floor.material.specular = 0.;

    let mut left_wall = Sphere::new();
    left_wall.transform = Matrix::translation(0., 0., 5.)
        * Matrix::rotation_y(-PI / 4.)
        * Matrix::rotation_x(PI / 2.)
        * Matrix::scaling(10., 0.01, 10.);
    left_wall.material = floor.material.clone();

    let mut right_wall = Sphere::new();
    right_wall.transform = Matrix::translation(0., 0., 5.)
        * Matrix::rotation_y(PI / 4.)
        * Matrix::rotation_x(PI / 2.)
        * Matrix::scaling(10., 0.01, 10.);
    right_wall.material = floor.material.clone();

    let mut middle = Sphere::new();
    middle.transform = Matrix::translation(-0.5, 1., 0.5);
    middle.material.color = Tuple::color(0.1, 1., 0.5);
    middle.material.diffuse = 0.7;
    middle.material.specular = 0.3;

    let mut right = Sphere::new();
    right.transform = Matrix::translation(1.5, 0.5, -0.5) * Matrix::scaling(0.5, 0.5, 0.5);
    right.material.color = Tuple::color(0.5, 1., 0.1);
    right.material.diffuse = 0.7;
    right.material.specular = 0.3;

    let mut left = Sphere::new();
    left.transform = Matrix::translation(-1.5, 0.33, -0.75) * Matrix::scaling(0.33, 0.33, 0.33);
    left.material.color = Tuple::color(1., 0.8, 0.1);
    left.material.diffuse = 0.7;
    left.material.specular = 0.3;

    let mut world = World::new();
    world.objects = vec![floor, left_wall, right_wall, middle, right, left];
    world.lights = vec![Light::new(Tuple::point(-10., 10., -10.), Tuple::color(1., 1., 1.))];

    let mut camera = Camera::new(256, 256, PI / 3.);
    camera.transform = view_transform(
        Tuple::point(0., 1.5, -5.),
        Tuple::point(0., 1., 0.),
        Tuple::vector(0., 1., 0.),
    );

    let canvas = camera.render(&world);
    canvas.write_ppm("scene.ppm").unwrap();
}