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

    // function lighting(material, light, point, eyev, normalv)
    //   # combine the surface color with the light's color/intensity
    //   effective_color ← material.color * light.intensity
    //
    //   # find the direction to the light source
    //   lightv ← normalize(light.position - point)
    //
    //   # compute the ambient contribution
    //   ambient ← effective_color * material.ambient
    //
    //   # light_dot_normal represents the cosine of the angle between the
    //   # light vector and the normal vector. A negative number means the
    //   # light is on the other side of the surface.
    //   light_dot_normal ← dot(lightv, normalv)
    //   if light_dot_normal < 0
    //     diffuse ← black
    //     specular ← black
    //
    //   else
    //     # compute the diffuse contribution
    //     diffuse ← effective_color * material.diffuse * light_dot_normal
    //
    //     # reflect_dot_eye represents the cosine of the angle between the
    //     # reflection vector and the eye vector. A negative number means the
    //     # light reflects away from the eye.
    //     reflectv ← reflect(-lightv, normalv)
    //     reflect_dot_eye ← dot(reflectv, eyev)
    //
    //     if reflect_dot_eye <= 0
    //       specular ← black
    //     else
    //       # compute the specular contribution
    //       factor ← pow(reflect_dot_eye, material.shininess)
    //       specular ← light.intensity * material.specular * factor
    //     end if
    //   end if
    //
    //   # Add the three contributions together to get the final shading
    //   return ambient + diffuse + specular
    // end function
    // function lighting(material, light, point, eyev, normalv, in_shadow)
    //   if material has a pattern
    //     color ← stripe_at(material.pattern, point)
    //   else
    //     color ← material.color
    //   end if
    //
    //   # then, compute the lighting as usual, using `color`
    //   # instead of `material.color`
    //
    //   # ...
    // end function
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

// Feature: Materials
// Background:
//   Given m ← material()
//     And position ← point(0, 0, 0)
#[cfg(test)]
mod tests {
    use std::f64::consts::SQRT_2;

    use super::*;

    use crate::{Sphere, StripePattern};
    use approx::assert_relative_eq;

    // Scenario: The default material
    //   Given m ← material()
    //   Then m.color = color(1, 1, 1)
    //     And m.ambient = 0.1
    //     And m.diffuse = 0.9
    //     And m.specular = 0.9
    //     And m.shininess = 200.0
    #[test]
    fn material_has_default_values() {
        let material = Material::default();
        assert_eq!(material.ambient, 0.1);
        assert_eq!(material.diffuse, 0.9);
        assert_eq!(material.specular, 0.9);
        assert_eq!(material.shininess, 200.0);
    }

    // Scenario: Lighting with the eye between the light and the surface
    //   Given eyev ← vector(0, 0, -1)
    //     And normalv ← vector(0, 0, -1)
    //     And light ← point_light(point(0, 0, -10), color(1, 1, 1))
    //   When result ← lighting(m, light, position, eyev, normalv)
    //   Then result = color(1.9, 1.9, 1.9)
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

    // Scenario: Lighting with the eye between light and surface, eye offset 45°
    //   Given eyev ← vector(0, √2/2, -√2/2)
    //     And normalv ← vector(0, 0, -1)
    //     And light ← point_light(point(0, 0, -10), color(1, 1, 1))
    //   When result ← lighting(m, light, position, eyev, normalv)
    //   Then result = color(1.0, 1.0, 1.0)
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

    // Scenario: Lighting with eye opposite surface, light offset 45°
    //   Given eyev ← vector(0, 0, -1)
    //     And normalv ← vector(0, 0, -1)
    //     And light ← point_light(point(0, 10, -10), color(1, 1, 1))
    //   When result ← lighting(m, light, position, eyev, normalv)
    //   Then result = color(0.7364, 0.7364, 0.7364)
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
            )
        );
    }

    // Scenario: Lighting with eye in the path of the reflection vector
    //   Given eyev ← vector(0, -√2/2, -√2/2)
    //     And normalv ← vector(0, 0, -1)
    //     And light ← point_light(point(0, 10, -10), color(1, 1, 1))
    //   When result ← lighting(m, light, position, eyev, normalv)
    //   Then result = color(1.6364, 1.6364, 1.6364)
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

    // Scenario: Lighting with the light behind the surface
    //   Given eyev ← vector(0, 0, -1)
    //     And normalv ← vector(0, 0, -1)
    //     And light ← point_light(point(0, 0, 10), color(1, 1, 1))
    //   When result ← lighting(m, light, position, eyev, normalv)
    //   Then result = color(0.1, 0.1, 0.1)
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

    // Scenario: Lighting with the surface in shadow
    //   Given eyev ← vector(0, 0, -1)
    //     And normalv ← vector(0, 0, -1)
    //     And light ← point_light(point(0, 0, -10), color(1, 1, 1))
    //     And in_shadow ← true
    //   When result ← lighting(m, light, position, eyev, normalv, in_shadow)
    //   Then result = color(0.1, 0.1, 0.1)
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

    // Scenario: Lighting with a pattern applied
    //   Given m.pattern ← stripe_pattern(color(1, 1, 1), color(0, 0, 0))
    //     And m.ambient ← 1
    //     And m.diffuse ← 0
    //     And m.specular ← 0
    //     And eyev ← vector(0, 0, -1)
    //     And normalv ← vector(0, 0, -1)
    //     And light ← point_light(point(0, 0, -10), color(1, 1, 1))
    //   When c1 ← lighting(m, light, point(0.9, 0, 0), eyev, normalv, false)
    //     And c2 ← lighting(m, light, point(1.1, 0, 0), eyev, normalv, false)
    //   Then c1 = color(1, 1, 1)
    //     And c2 = color(0, 0, 0)
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

    // Scenario: Reflectivity for the default material
    //   Given m ← material()
    //   Then m.reflective = 0.0
    #[test]
    fn reflectivity_for_default_material() {
        let material = Material::default();
        assert_eq!(material.reflective, 0.0);
    }

    // Scenario: Transparency and Refractive Index for the default material
    //   Given m ← material()
    //   Then m.transparency = 0.0
    //     And m.refractive_index = 1.0
    #[test]
    fn transparency_and_refractive_index_for_default_material() {
        let material = Material::default();
        assert_eq!(material.transparency, 0.0);
        assert_eq!(material.refractive_index, 1.0);
    }
}
