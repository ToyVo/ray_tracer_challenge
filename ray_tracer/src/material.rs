use crate::{Light, Pattern, Shape, SolidPattern, Tuple};

#[derive(Debug, Clone)]
pub struct Material {
    pub pattern: Box<dyn Pattern>,
    pub ambient: f64,
    pub diffuse: f64,
    pub specular: f64,
    pub shininess: f64,
    pub reflective: f64,
    pub transparency: f64,
    pub refractive_index: f64,
}

impl PartialEq for Material {
    fn eq(&self, other: &Self) -> bool {
        self.pattern.colors() == other.pattern.colors()
            && self.ambient == other.ambient
            && self.diffuse == other.diffuse
            && self.specular == other.specular
            && self.shininess == other.shininess
    }
}

impl Default for Material {
    fn default() -> Self {
        Material {
            pattern: Box::new(SolidPattern::new(Tuple::color(1.0, 1.0, 1.0))),
            ambient: 0.1,
            diffuse: 0.9,
            specular: 0.9,
            shininess: 200.0,
            reflective: 0.0,
            transparency: 0.0,
            refractive_index: 1.0,
        }
    }
}

impl Material {
    pub fn default_with_pattern(pattern: Box<dyn Pattern>) -> Material {
        Material {
            pattern,
            ambient: 0.1,
            diffuse: 0.9,
            specular: 0.9,
            shininess: 200.0,
            reflective: 0.0,
            transparency: 0.0,
            refractive_index: 1.0,
        }
    }

    pub fn lighting(
        &self,
        object: &dyn Shape,
        light: &Light,
        point: &Tuple,
        eye_vector: &Tuple,
        normal_vector: &Tuple,
        in_shadow: bool,
    ) -> Tuple {
        let effective_color = self.pattern.pattern_at_object(object, point) * &light.intensity;
        let ambient = &effective_color * self.ambient;
        if in_shadow {
            return ambient;
        }
        let light_vector = (&light.position - point).normalize();
        let light_dot_normal = light_vector.dot(normal_vector);
        let (diffuse, specular) = if light_dot_normal < 0.0 {
            (Tuple::color(0.0, 0.0, 0.0), Tuple::color(0.0, 0.0, 0.0))
        } else {
            let diffuse = &effective_color * self.diffuse * light_dot_normal;
            let reflect_vector = -light_vector.reflect(normal_vector);
            let reflect_dot_eye = reflect_vector.dot(eye_vector);
            let specular = if reflect_dot_eye <= 0.0 {
                Tuple::color(0.0, 0.0, 0.0)
            } else {
                let factor = reflect_dot_eye.powf(self.shininess);
                &light.intensity * self.specular * factor
            };
            (diffuse, specular)
        };
        ambient + diffuse + specular
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::SQRT_2;

    use super::*;

    use crate::{Sphere, StripePattern};
    use approx::assert_relative_eq;

    #[test]
    fn material_has_default_values() {
        let material = Material::default();
        assert_eq!(material.ambient, 0.1);
        assert_eq!(material.diffuse, 0.9);
        assert_eq!(material.specular, 0.9);
        assert_eq!(material.shininess, 200.0);
    }

    #[test]
    fn lighting_eye_between_light_and_surface() {
        let material = Material::default();
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eye_vector = Tuple::vector(0.0, 0.0, -1.0);
        let normal_vector = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::new(Tuple::point(0.0, 0.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        let shape = Sphere::new(0);
        let result = material.lighting(
            &shape,
            &light,
            &position,
            &eye_vector,
            &normal_vector,
            false,
        );
        assert_eq!(result, Tuple::color(1.9, 1.9, 1.9));
    }

    #[test]
    fn lighting_eye_between_light_and_surface_eye_offset_45_degrees() {
        let material = Material::default();
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eye_vector = Tuple::vector(0.0, SQRT_2 / 2.0, -SQRT_2 / 2.0);
        let normal_vector = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::new(Tuple::point(0.0, 0.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        let shape = Sphere::new(0);
        let result = material.lighting(
            &shape,
            &light,
            &position,
            &eye_vector,
            &normal_vector,
            false,
        );
        assert_eq!(result, Tuple::color(1.0, 1.0, 1.0));
    }

    #[test]
    fn lighting_eye_opposite_surface_light_offset_45_degrees() {
        let material = Material::default();
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eye_vector = Tuple::vector(0.0, 0.0, -1.0);
        let normal_vector = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::new(Tuple::point(0.0, 10.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        let shape = Sphere::new(0);
        let result = material.lighting(
            &shape,
            &light,
            &position,
            &eye_vector,
            &normal_vector,
            false,
        );
        assert_relative_eq!(
            result,
            Tuple::color(
                0.9 * SQRT_2 / 2. + 0.1,
                0.9 * SQRT_2 / 2. + 0.1,
                0.9 * SQRT_2 / 2. + 0.1
            ),
            epsilon = 1e-5f64
        );
    }

    #[test]
    fn lighting_eye_in_path_of_reflection_vector() {
        let material = Material::default();
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eye_vector = Tuple::vector(0.0, -SQRT_2 / 2.0, -SQRT_2 / 2.0);
        let normal_vector = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::new(Tuple::point(0.0, 10.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        let shape = Sphere::new(0);
        let result = material.lighting(
            &shape,
            &light,
            &position,
            &eye_vector,
            &normal_vector,
            false,
        );
        assert_eq!(
            result,
            Tuple::color(
                0.9 * SQRT_2 / 2. + 1.,
                0.9 * SQRT_2 / 2. + 1.,
                0.9 * SQRT_2 / 2. + 1.
            )
        );
    }

    #[test]
    fn lighting_light_behind_surface() {
        let material = Material::default();
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eye_vector = Tuple::vector(0.0, 0.0, -1.0);
        let normal_vector = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::new(Tuple::point(0.0, 0.0, 10.0), Tuple::color(1.0, 1.0, 1.0));
        let shape = Sphere::new(0);
        let result = material.lighting(
            &shape,
            &light,
            &position,
            &eye_vector,
            &normal_vector,
            false,
        );
        assert_eq!(result, Tuple::color(0.1, 0.1, 0.1));
    }

    #[test]
    fn lighting_surface_in_shadow() {
        let material = Material::default();
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eye_vector = Tuple::vector(0.0, 0.0, -1.0);
        let normal_vector = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::new(Tuple::point(0.0, 0.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        let shape = Sphere::new(0);
        let result =
            material.lighting(&shape, &light, &position, &eye_vector, &normal_vector, true);
        assert_eq!(result, Tuple::color(0.1, 0.1, 0.1));
    }

    #[test]
    fn lighting_with_pattern_applied() {
        let mut material = Material::default_with_pattern(Box::new(StripePattern::new(
            Box::new(SolidPattern::new(Tuple::color(1.0, 1.0, 1.0))),
            Box::new(SolidPattern::new(Tuple::color(0.0, 0.0, 0.0))),
        )));
        material.ambient = 1.0;
        material.diffuse = 0.0;
        material.specular = 0.0;
        let eye_vector = Tuple::vector(0.0, 0.0, -1.0);
        let normal_vector = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::new(Tuple::point(0.0, 0.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        let shape = Sphere::new(0);
        let color_a = material.lighting(
            &shape,
            &light,
            &Tuple::point(0.9, 0.0, 0.0),
            &eye_vector,
            &normal_vector,
            false,
        );
        let color_b = material.lighting(
            &shape,
            &light,
            &Tuple::point(1.1, 0.0, 0.0),
            &eye_vector,
            &normal_vector,
            false,
        );
        assert_eq!(color_a, Tuple::color(1.0, 1.0, 1.0));
        assert_eq!(color_b, Tuple::color(0.0, 0.0, 0.0));
    }

    #[test]
    fn reflectivity_for_default_material() {
        let material = Material::default();
        assert_eq!(material.reflective, 0.0);
    }

    #[test]
    fn transparency_and_refractive_index_for_default_material() {
        let material = Material::default();
        assert_eq!(material.transparency, 0.0);
        assert_eq!(material.refractive_index, 1.0);
    }
}
