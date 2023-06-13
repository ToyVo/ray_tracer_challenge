use std::f64::consts::PI;
use std::rc::Rc;

use ray_tracer_challenge::{Camera, Light, Sphere, Tuple, view_transform, World, Shape, Plane, Matrix};

fn main() {
    let mut floor = Plane::new();
    *floor.transform_mut() = Matrix::scaling(10., 0.01, 10.);
    floor.material_mut().color = Tuple::color(1., 0.9, 0.9);
    floor.material_mut().specular = 0.;

    let mut left_wall = Plane::new();
    *left_wall.transform_mut() = Matrix::translation(0., 0., 5.) * Matrix::rotation_y(-PI / 4.) * Matrix::rotation_x(PI / 2.) * Matrix::scaling(10., 0.01, 10.);
    left_wall.material_mut().color = Tuple::color(1., 0.9, 0.9);
    left_wall.material_mut().specular = 0.;

    let mut right_wall = Plane::new();
    *right_wall.transform_mut() = Matrix::translation(0., 0., 5.) * Matrix::rotation_y(PI / 4.) * Matrix::rotation_x(PI / 2.) * Matrix::scaling(10., 0.01, 10.);
    right_wall.material_mut().color = Tuple::color(1., 0.9, 0.9);
    right_wall.material_mut().specular = 0.;

    let mut middle = Sphere::new();
    *middle.transform_mut() = Matrix::translation(-0.5, 1., 0.5);
    middle.material_mut().color = Tuple::color(0.1, 1., 0.5);
    middle.material_mut().diffuse = 0.7;
    middle.material_mut().specular = 0.3;


    let mut right = Sphere::new();
    *right.transform_mut() = Matrix::translation(1.5, 0.5, -0.5) * Matrix::scaling(0.5, 0.5, 0.5);
    right.material_mut().color = Tuple::color(0.5, 1., 0.1);
    right.material_mut().diffuse = 0.7;
    right.material_mut().specular = 0.3;

    let mut left = Sphere::new();
    *left.transform_mut() = Matrix::translation(-1.5, 0.33, -0.75) * Matrix::scaling(0.33, 0.33, 0.33);
    left.material_mut().color = Tuple::color(1., 0.8, 0.1);
    left.material_mut().diffuse = 0.7;
    left.material_mut().specular = 0.3;

    let mut world = World::new();
    world.objects = vec![Rc::new(floor), Rc::new(left_wall), Rc::new(right_wall), Rc::new(middle), Rc::new(right), Rc::new(left)];
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
