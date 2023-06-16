pub use camera::*;
pub use canvas::*;
pub use intersection::*;
pub use light::*;
pub use material::*;
pub use matrix::*;
pub use pattern::*;
pub use plane::*;
pub use ray::*;
pub use shape::*;
pub use sphere::*;
pub use transformations::*;
pub use tuple::*;
pub use world::*;

mod camera;
mod canvas;
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
