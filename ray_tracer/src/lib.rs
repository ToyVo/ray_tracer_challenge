pub use camera::Camera;
pub use canvas::Canvas;
pub use cone::Cone;
pub use counter::Counter;
pub use cube::Cube;
pub use cylinder::Cylinder;
pub use intersection::{Computations, Intersection};
pub use light::Light;
pub use material::Material;
pub use pattern::{
    BlendedPattern, CheckeredPattern, GradientPattern, Pattern, PerturbedPattern, RingPattern,
    SolidPattern, StripePattern,
};
pub use plane::Plane;
pub use ray::Ray;
pub use shape::Shape;
pub use sphere::Sphere;
pub use transformations::view_transform;
pub use world::World;

#[cfg(test)]
pub use pattern::tests::TestPattern;

mod camera;
mod canvas;
mod cone;
mod counter;
mod cube;
mod cylinder;
mod intersection;
mod light;
mod material;
mod pattern;
mod plane;
mod ray;
mod shape;
mod sphere;
mod transformations;
mod world;

use nalgebra_glm::DMat4;

pub trait Transform {
    fn transform(&self) -> &DMat4;
    fn transform_mut(&mut self) -> &mut DMat4;
}
