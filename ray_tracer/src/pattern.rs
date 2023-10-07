use crate::{Shape, Transform};
use bracket_noise::prelude::*;
use dyn_clone::DynClone;
use nalgebra_glm::{inverse, vec4, DMat4, DVec3, DVec4};

pub trait Pattern: std::fmt::Debug + DynClone + Transform {
    fn pattern_at(&self, point: &DVec4) -> DVec3;
    fn colors(&self) -> Vec<&DVec3>;
    fn pattern_at_object(&self, object: &dyn Shape, world_point: &DVec4) -> DVec3 {
        let target_point = inverse(object.transform()) * world_point;
        let pattern_point = inverse(self.transform()) * target_point;
        self.pattern_at(&pattern_point)
    }
    fn pattern_at_pattern(&self, pattern: &dyn Pattern, world_point: &DVec4) -> DVec3 {
        let target_point = inverse(pattern.transform()) * world_point;
        let pattern_point = inverse(self.transform()) * target_point;
        self.pattern_at(&pattern_point)
    }
}

dyn_clone::clone_trait_object!(Pattern);

impl PartialEq for dyn Pattern {
    fn eq(&self, other: &Self) -> bool {
        self.colors() == other.colors()
    }
}

#[derive(Debug, Clone)]
pub struct StripePattern {
    pattern_a: Box<dyn Pattern>,
    pattern_b: Box<dyn Pattern>,
    transform: DMat4,
}

impl StripePattern {
    pub fn new(pattern_a: Box<dyn Pattern>, pattern_b: Box<dyn Pattern>) -> StripePattern {
        StripePattern {
            pattern_a,
            pattern_b,
            transform: DMat4::identity(),
        }
    }
}

impl Pattern for StripePattern {
    fn pattern_at(&self, point: &DVec4) -> DVec3 {
        if point.x.floor() as i32 % 2 == 0 {
            self.pattern_a.pattern_at_pattern(self, point)
        } else {
            self.pattern_b.pattern_at_pattern(self, point)
        }
    }

    fn colors(&self) -> Vec<&DVec3> {
        let mut colors = self.pattern_a.colors();
        colors.append(&mut self.pattern_b.colors());
        colors
    }
}

impl Transform for StripePattern {
    fn transform(&self) -> &DMat4 {
        &self.transform
    }

    fn transform_mut(&mut self) -> &mut DMat4 {
        &mut self.transform
    }
}

#[derive(Debug, Clone)]
pub struct SolidPattern {
    color: DVec3,
    transform: DMat4,
}

impl SolidPattern {
    pub fn new(color: DVec3) -> SolidPattern {
        SolidPattern {
            color,
            transform: DMat4::identity(),
        }
    }
}

impl Pattern for SolidPattern {
    fn pattern_at(&self, _point: &DVec4) -> DVec3 {
        self.color
    }

    fn colors(&self) -> Vec<&DVec3> {
        vec![&self.color]
    }
}

impl Transform for SolidPattern {
    fn transform(&self) -> &DMat4 {
        &self.transform
    }

    fn transform_mut(&mut self) -> &mut DMat4 {
        &mut self.transform
    }
}

#[derive(Debug, Clone)]
pub struct GradientPattern {
    color_a: DVec3,
    color_b: DVec3,
    transform: DMat4,
}

impl GradientPattern {
    pub fn new(color_a: DVec3, color_b: DVec3) -> GradientPattern {
        GradientPattern {
            color_a,
            color_b,
            transform: DMat4::identity(),
        }
    }
}

impl Pattern for GradientPattern {
    fn pattern_at(&self, point: &DVec4) -> DVec3 {
        let distance = self.color_b - self.color_a;
        let fraction = point.x - point.x.floor();
        self.color_a + distance * fraction
    }

    fn colors(&self) -> Vec<&DVec3> {
        vec![&self.color_a, &self.color_b]
    }
}

impl Transform for GradientPattern {
    fn transform(&self) -> &DMat4 {
        &self.transform
    }

    fn transform_mut(&mut self) -> &mut DMat4 {
        &mut self.transform
    }
}

#[derive(Debug, Clone)]
pub struct RingPattern {
    pattern_a: Box<dyn Pattern>,
    pattern_b: Box<dyn Pattern>,
    transform: DMat4,
}

impl RingPattern {
    pub fn new(pattern_a: Box<dyn Pattern>, pattern_b: Box<dyn Pattern>) -> RingPattern {
        RingPattern {
            pattern_a,
            pattern_b,
            transform: DMat4::identity(),
        }
    }
}

impl Pattern for RingPattern {
    fn pattern_at(&self, point: &DVec4) -> DVec3 {
        if (point.x.powi(2) + point.z.powi(2)).sqrt().floor() as i32 % 2 == 0 {
            self.pattern_a.pattern_at_pattern(self, point)
        } else {
            self.pattern_b.pattern_at_pattern(self, point)
        }
    }

    fn colors(&self) -> Vec<&DVec3> {
        let mut colors = self.pattern_a.colors();
        colors.append(&mut self.pattern_b.colors());
        colors
    }
}

impl Transform for RingPattern {
    fn transform(&self) -> &DMat4 {
        &self.transform
    }

    fn transform_mut(&mut self) -> &mut DMat4 {
        &mut self.transform
    }
}

#[derive(Debug, Clone)]
pub struct CheckeredPattern {
    pattern_a: Box<dyn Pattern>,
    pattern_b: Box<dyn Pattern>,
    transform: DMat4,
}

impl CheckeredPattern {
    pub fn new(pattern_a: Box<dyn Pattern>, pattern_b: Box<dyn Pattern>) -> CheckeredPattern {
        CheckeredPattern {
            pattern_a,
            pattern_b,
            transform: DMat4::identity(),
        }
    }
}

impl Pattern for CheckeredPattern {
    fn pattern_at(&self, point: &DVec4) -> DVec3 {
        if (point.x.floor() + point.y.floor() + point.z.floor()) as i32 % 2 == 0 {
            self.pattern_a.pattern_at_pattern(self, point)
        } else {
            self.pattern_b.pattern_at_pattern(self, point)
        }
    }

    fn colors(&self) -> Vec<&DVec3> {
        let mut colors = self.pattern_a.colors();
        colors.append(&mut self.pattern_b.colors());
        colors
    }
}

impl Transform for CheckeredPattern {
    fn transform(&self) -> &DMat4 {
        &self.transform
    }

    fn transform_mut(&mut self) -> &mut DMat4 {
        &mut self.transform
    }
}

#[derive(Debug, Clone)]
pub struct BlendedPattern {
    pattern_a: Box<dyn Pattern>,
    pattern_b: Box<dyn Pattern>,
    transform: DMat4,
}

impl BlendedPattern {
    pub fn new(pattern_a: Box<dyn Pattern>, pattern_b: Box<dyn Pattern>) -> BlendedPattern {
        BlendedPattern {
            pattern_a,
            pattern_b,
            transform: DMat4::identity(),
        }
    }
}

impl Pattern for BlendedPattern {
    fn pattern_at(&self, point: &DVec4) -> DVec3 {
        let color_a = self.pattern_a.pattern_at_pattern(self, point);
        let color_b = self.pattern_b.pattern_at_pattern(self, point);
        (color_a + color_b) / 2.
    }

    fn colors(&self) -> Vec<&DVec3> {
        let mut colors = self.pattern_a.colors();
        colors.append(&mut self.pattern_b.colors());
        colors
    }
}

impl Transform for BlendedPattern {
    fn transform(&self) -> &DMat4 {
        &self.transform
    }

    fn transform_mut(&mut self) -> &mut DMat4 {
        &mut self.transform
    }
}

#[derive(Debug, Clone)]
pub struct PerturbedPattern {
    pattern: Box<dyn Pattern>,
    transform: DMat4,
}

impl PerturbedPattern {
    pub fn new(pattern: Box<dyn Pattern>) -> PerturbedPattern {
        PerturbedPattern {
            pattern,
            transform: DMat4::identity(),
        }
    }
}

impl Pattern for PerturbedPattern {
    fn pattern_at(&self, point: &DVec4) -> DVec3 {
        let mut noise = FastNoise::new();
        noise.set_noise_type(NoiseType::SimplexFractal);
        noise.set_fractal_octaves(3);
        noise.set_fractal_gain(0.8);
        let scale = 0.2;
        let p_scale = 1.0;

        let x = point.x as f32 * p_scale;
        let y = point.y as f32 * p_scale;
        let z = point.z as f32 * p_scale;

        let noise_x = noise.get_noise3d(x, y, z) as f64 * scale;
        let noise_y = noise.get_noise3d(x, y, z + 1.) as f64 * scale;
        let noise_z = noise.get_noise3d(x, y, z + 2.) as f64 * scale;

        let new_x = point.x + noise_x;
        let new_y = point.y + noise_y;
        let new_z = point.z + noise_z;

        let new_point = vec4(new_x, new_y, new_z, 1.);

        self.pattern.pattern_at_pattern(self, &new_point)
    }

    fn colors(&self) -> Vec<&DVec3> {
        self.pattern.colors()
    }
}

impl Transform for PerturbedPattern {
    fn transform(&self) -> &DMat4 {
        &self.transform
    }

    fn transform_mut(&mut self) -> &mut DMat4 {
        &mut self.transform
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use crate::Sphere;
    use nalgebra_glm::{scale, translate, vec3};

    #[derive(Debug, Clone)]
    pub struct TestPattern {
        transform: DMat4,
        color: DVec3,
    }

    impl Default for TestPattern {
        fn default() -> Self {
            TestPattern {
                transform: DMat4::identity(),
                color: vec3(0., 0., 0.),
            }
        }
    }

    impl Pattern for TestPattern {
        fn pattern_at(&self, point: &DVec4) -> DVec3 {
            nalgebra_glm::vec4_to_vec3(point)
        }

        fn colors(&self) -> Vec<&DVec3> {
            vec![&self.color]
        }
    }

    impl Transform for TestPattern {
        fn transform(&self) -> &DMat4 {
            &self.transform
        }

        fn transform_mut(&mut self) -> &mut DMat4 {
            &mut self.transform
        }
    }

    #[test]
    fn creating_stripe_pattern() {
        let white = vec3(1.0, 1.0, 1.0);
        let black = vec3(0.0, 0.0, 0.0);
        let pattern = StripePattern::new(
            Box::new(SolidPattern::new(white)),
            Box::new(SolidPattern::new(black)),
        );
        assert_eq!(pattern.pattern_a.colors(), vec![&white]);
        assert_eq!(pattern.pattern_b.colors(), vec![&black]);
    }

    #[test]
    fn stripe_pattern_constant_in_y() {
        let white = vec3(1.0, 1.0, 1.0);
        let black = vec3(0.0, 0.0, 0.0);
        let pattern = StripePattern::new(
            Box::new(SolidPattern::new(white)),
            Box::new(SolidPattern::new(black)),
        );
        assert_eq!(pattern.pattern_at(&vec4(0.0, 0.0, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&vec4(0.0, 1.0, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&vec4(0.0, 2.0, 0.0, 1.)), white);
    }

    #[test]
    fn stripe_pattern_constant_in_z() {
        let white = vec3(1.0, 1.0, 1.0);
        let black = vec3(0.0, 0.0, 0.0);
        let pattern = StripePattern::new(
            Box::new(SolidPattern::new(white)),
            Box::new(SolidPattern::new(black)),
        );
        assert_eq!(pattern.pattern_at(&vec4(0.0, 0.0, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&vec4(0.0, 0.0, 1.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&vec4(0.0, 0.0, 2.0, 1.)), white);
    }

    #[test]
    fn stripe_pattern_alternates_in_x() {
        let white = vec3(1.0, 1.0, 1.0);
        let black = vec3(0.0, 0.0, 0.0);
        let pattern = StripePattern::new(
            Box::new(SolidPattern::new(white)),
            Box::new(SolidPattern::new(black)),
        );
        assert_eq!(pattern.pattern_at(&vec4(0.0, 0.0, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&vec4(0.9, 0.0, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&vec4(1.0, 0.0, 0.0, 1.)), black);
        assert_eq!(pattern.pattern_at(&vec4(-0.1, 0.0, 0.0, 1.)), black);
        assert_eq!(pattern.pattern_at(&vec4(-1.0, 0.0, 0.0, 1.)), black);
        assert_eq!(pattern.pattern_at(&vec4(-1.1, 0.0, 0.0, 1.)), white);
    }

    #[test]
    fn stripe_pattern_with_object_transformation() {
        let white = vec3(1.0, 1.0, 1.0);
        let black = vec3(0.0, 0.0, 0.0);
        let mut shape = Sphere::new(0);
        *shape.transform_mut() = scale(&DMat4::identity(), &vec3(2.0, 2.0, 2.0));
        let pattern = StripePattern::new(
            Box::new(SolidPattern::new(white)),
            Box::new(SolidPattern::new(black)),
        );
        let color = pattern.pattern_at_object(&shape, &vec4(1.5, 0.0, 0.0, 1.));
        assert_eq!(color, white);
    }

    #[test]
    fn stripe_pattern_with_pattern_transformation() {
        let white = vec3(1.0, 1.0, 1.0);
        let black = vec3(0.0, 0.0, 0.0);
        let shape = Sphere::new(0);
        let mut pattern = StripePattern::new(
            Box::new(SolidPattern::new(white)),
            Box::new(SolidPattern::new(black)),
        );
        *pattern.transform_mut() = scale(&DMat4::identity(), &vec3(2.0, 2.0, 2.0));
        let color = pattern.pattern_at_object(&shape, &vec4(1.5, 0.0, 0.0, 1.));
        assert_eq!(color, white);
    }

    #[test]
    fn stripe_pattern_with_both_object_and_pattern_transformation() {
        let white = vec3(1.0, 1.0, 1.0);
        let black = vec3(0.0, 0.0, 0.0);
        let mut shape = Sphere::new(0);
        *shape.transform_mut() = scale(&DMat4::identity(), &vec3(2.0, 2.0, 2.0));
        let mut pattern = StripePattern::new(
            Box::new(SolidPattern::new(white)),
            Box::new(SolidPattern::new(black)),
        );
        *pattern.transform_mut() = translate(&DMat4::identity(), &vec3(0.5, 0.0, 0.0));
        let color = pattern.pattern_at_object(&shape, &vec4(2.5, 0.0, 0.0, 1.));
        assert_eq!(color, white);
    }

    #[test]
    fn default_pattern_transformation() {
        let pattern = TestPattern::default();
        assert_eq!(*pattern.transform(), DMat4::identity());
    }

    #[test]
    fn assigning_transformation() {
        let mut pattern = TestPattern::default();
        *pattern.transform_mut() = translate(&DMat4::identity(), &vec3(1.0, 2.0, 3.0));
        assert_eq!(
            *pattern.transform(),
            translate(&DMat4::identity(), &vec3(1.0, 2.0, 3.0))
        );
    }

    #[test]
    fn pattern_with_object_transformation() {
        let mut shape = Sphere::new(0);
        *shape.transform_mut() = scale(&DMat4::identity(), &vec3(2.0, 2.0, 2.0));
        let pattern = TestPattern::default();
        let color = pattern.pattern_at_object(&shape, &vec4(2.0, 3.0, 4.0, 1.));
        assert_eq!(color, vec3(1.0, 1.5, 2.0));
    }

    #[test]
    fn pattern_with_pattern_transformation() {
        let shape = Sphere::new(0);
        let mut pattern = TestPattern::default();
        *pattern.transform_mut() = scale(&DMat4::identity(), &vec3(2.0, 2.0, 2.0));
        let color = pattern.pattern_at_object(&shape, &vec4(2.0, 3.0, 4.0, 1.));
        assert_eq!(color, vec3(1.0, 1.5, 2.0));
    }

    #[test]
    fn pattern_with_both_object_and_pattern_transformation() {
        let mut shape = Sphere::new(0);
        *shape.transform_mut() = scale(&DMat4::identity(), &vec3(2.0, 2.0, 2.0));
        let mut pattern = TestPattern::default();
        *pattern.transform_mut() = translate(&DMat4::identity(), &vec3(0.5, 1.0, 1.5));
        let color = pattern.pattern_at_object(&shape, &vec4(2.5, 3.0, 3.5, 1.));
        assert_eq!(color, vec3(0.75, 0.5, 0.25));
    }

    #[test]
    fn gradient_linearly_interpolates_between_colors() {
        let white = vec3(1.0, 1.0, 1.0);
        let black = vec3(0.0, 0.0, 0.0);
        let pattern = GradientPattern::new(white, black);
        assert_eq!(pattern.pattern_at(&vec4(0.0, 0.0, 0.0, 1.)), white);
        assert_eq!(
            pattern.pattern_at(&vec4(0.25, 0.0, 0.0, 1.)),
            vec3(0.75, 0.75, 0.75)
        );
        assert_eq!(
            pattern.pattern_at(&vec4(0.5, 0.0, 0.0, 1.)),
            vec3(0.5, 0.5, 0.5)
        );
        assert_eq!(
            pattern.pattern_at(&vec4(0.75, 0.0, 0.0, 1.)),
            vec3(0.25, 0.25, 0.25)
        );
    }

    #[test]
    fn ring_should_extend_in_both_x_and_z() {
        let white = vec3(1.0, 1.0, 1.0);
        let black = vec3(0.0, 0.0, 0.0);
        let pattern = RingPattern::new(
            Box::new(SolidPattern::new(white)),
            Box::new(SolidPattern::new(black)),
        );
        assert_eq!(pattern.pattern_at(&vec4(0.0, 0.0, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&vec4(1.0, 0.0, 0.0, 1.)), black);
        assert_eq!(pattern.pattern_at(&vec4(0.0, 0.0, 1.0, 1.)), black);
        assert_eq!(pattern.pattern_at(&vec4(0.708, 0.0, 0.708, 1.)), black);
    }

    #[test]
    fn checkers_should_repeat_in_x() {
        let white = vec3(1.0, 1.0, 1.0);
        let black = vec3(0.0, 0.0, 0.0);
        let pattern = CheckeredPattern::new(
            Box::new(SolidPattern::new(white)),
            Box::new(SolidPattern::new(black)),
        );
        assert_eq!(pattern.pattern_at(&vec4(0.0, 0.0, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&vec4(0.99, 0.0, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&vec4(1.01, 0.0, 0.0, 1.)), black);
    }

    #[test]
    fn checkers_should_repeat_in_y() {
        let white = vec3(1.0, 1.0, 1.0);
        let black = vec3(0.0, 0.0, 0.0);
        let pattern = CheckeredPattern::new(
            Box::new(SolidPattern::new(white)),
            Box::new(SolidPattern::new(black)),
        );
        assert_eq!(pattern.pattern_at(&vec4(0.0, 0.0, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&vec4(0.0, 0.99, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&vec4(0.0, 1.01, 0.0, 1.)), black);
    }

    #[test]
    fn checkers_should_repeat_in_z() {
        let white = vec3(1.0, 1.0, 1.0);
        let black = vec3(0.0, 0.0, 0.0);
        let pattern = CheckeredPattern::new(
            Box::new(SolidPattern::new(white)),
            Box::new(SolidPattern::new(black)),
        );
        assert_eq!(pattern.pattern_at(&vec4(0.0, 0.0, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&vec4(0.0, 0.0, 0.99, 1.)), white);
        assert_eq!(pattern.pattern_at(&vec4(0.0, 0.0, 1.01, 1.)), black);
    }
}
