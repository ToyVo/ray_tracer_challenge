pub use camera::Camera;
pub use canvas::Canvas;
pub use counter::Counter;
pub use cube::Cube;
pub use intersection::{Computations, Intersection};
pub use light::Light;
pub use material::Material;
pub use matrix::Matrix;
pub use pattern::{
    BlendedPattern, CheckeredPattern, GradientPattern, Pattern, PerturbedPattern, RingPattern,
    SolidPattern, StripePattern,
};
pub use plane::Plane;
pub use ray::Ray;
pub use shape::Shape;
pub use sphere::Sphere;
pub use transformations::view_transform;
pub use tuple::Tuple;
pub use world::World;

#[cfg(test)]
pub use pattern::tests::TestPattern;

mod camera;
mod canvas;
mod counter;
mod cube;
mod intersection;
mod light;
mod material;
mod matrix;
mod pattern;
mod plane;
mod ray;
mod shape;
mod sphere;
mod transformations;
mod tuple;
mod world;

pub trait Transform {
    fn transform(&self) -> &Matrix;
    fn transform_mut(&mut self) -> &mut Matrix;
}
