use std::f64::consts::PI;

use ray_tracer::{
    view_transform, BlendedPattern, Camera, CheckeredPattern, Light, Matrix, PerturbedPattern,
    Plane, Shape, SolidPattern, Sphere, StripePattern, Transform, Tuple, World,
};

fn main() {
    let mut world = World::new();
    let mut s1 = StripePattern::new(
        Box::new(SolidPattern::new(Tuple::color(1., 1., 1.))),
        Box::new(SolidPattern::new(Tuple::color(0.75, 0.75, 0.75))),
    );
    *s1.transform_mut() = Matrix::scaling(0.1, 0.1, 0.1);
    let mut s2 = StripePattern::new(
        Box::new(SolidPattern::new(Tuple::color(0.5, 0.5, 0.5))),
        Box::new(SolidPattern::new(Tuple::color(0.25, 0.25, 0.25))),
    );
    *s2.transform_mut() = Matrix::scaling(0.1, 0.1, 0.1) * Matrix::rotation_y(PI / 2.);

    let mut floor = Plane::new(world.counter.increment());
    *floor.transform_mut() = Matrix::scaling(10., 0.01, 10.);
    floor.material_mut().pattern = Box::new(PerturbedPattern::new(Box::new(BlendedPattern::new(
        Box::new(s1),
        Box::new(s2),
    ))));
    floor.material_mut().specular = 0.;

    let mut middle = Sphere::new(world.counter.increment());
    *middle.transform_mut() = Matrix::translation(-0.5, 1., 0.5);
    middle.material_mut().pattern = Box::new(CheckeredPattern::new(
        Box::new(SolidPattern::new(Tuple::color(0.1, 1., 0.5))),
        Box::new(SolidPattern::new(Tuple::color(0.2, 0.5, 1.))),
    ));
    *middle.material_mut().pattern.transform_mut() = Matrix::scaling(0.2, 0.2, 0.2);
    middle.material_mut().diffuse = 0.7;
    middle.material_mut().specular = 0.3;

    let mut right = Sphere::new(world.counter.increment());
    *right.transform_mut() = Matrix::translation(1.5, 0.5, -0.5) * Matrix::scaling(0.5, 0.5, 0.5);
    right.material_mut().pattern = Box::new(SolidPattern::new(Tuple::color(0.5, 1., 0.1)));
    right.material_mut().diffuse = 0.7;
    right.material_mut().specular = 0.3;

    let mut left = Sphere::new(world.counter.increment());
    *left.transform_mut() =
        Matrix::translation(-1.5, 0.33, -0.75) * Matrix::scaling(0.33, 0.33, 0.33);
    left.material_mut().pattern = Box::new(SolidPattern::new(Tuple::color(1., 0.8, 0.1)));
    left.material_mut().diffuse = 0.7;
    left.material_mut().specular = 0.3;

    world.objects = vec![
        Box::new(floor),
        Box::new(middle),
        Box::new(right),
        Box::new(left),
    ];
    world.lights = vec![Light::new(
        Tuple::point(-10., 10., -10.),
        Tuple::color(1., 1., 1.),
    )];

    let mut camera = Camera::new(256, 256, PI / 3.);
    camera.transform = view_transform(
        Tuple::point(0., 1.5, -5.),
        Tuple::point(0., 1., 0.),
        Tuple::vector(0., 1., 0.),
    );

    let canvas = camera.render(&world);
    canvas.write_ppm("scene.ppm").unwrap();
}
