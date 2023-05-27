use crate::{Light, Tuple};

#[derive(Debug, PartialEq, Clone)]
pub struct Material {
    pub color: Tuple,
    pub ambient: f64,
    pub diffuse: f64,
    pub specular: f64,
    pub shininess: f64,
}

impl Material {
    pub fn new() -> Material {
        Material {
            color: Tuple::color(1.0, 1.0, 1.0),
            ambient: 0.1,
            diffuse: 0.9,
            specular: 0.9,
            shininess: 200.0,
        }
    }

    pub fn lighting(&self, light: &Light, point: &Tuple, eyev: &Tuple, normalv: &Tuple) -> Tuple {
        let effective_color = &self.color * &light.intensity;
        let lightv = (&light.position - point).normalize();
        let ambient = &effective_color * self.ambient;
        let light_dot_normal = lightv.dot(normalv);
        let (diffuse, specular) = if light_dot_normal < 0.0 {
            (Tuple::color(0.0, 0.0, 0.0), Tuple::color(0.0, 0.0, 0.0))
        } else {
            let diffuse = &effective_color * self.diffuse * light_dot_normal;
            let reflectv = -lightv.reflect(normalv);
            let reflect_dot_eye = reflectv.dot(eyev);
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
    use super::*;
    use std::f64::consts::SQRT_2;

    #[test]
    fn material_has_default_values() {
        let m = Material::new();
        assert_eq!(m.color, Tuple::color(1.0, 1.0, 1.0));
        assert_eq!(m.ambient, 0.1);
        assert_eq!(m.diffuse, 0.9);
        assert_eq!(m.specular, 0.9);
        assert_eq!(m.shininess, 200.0);
    }

    #[test]
    fn lighting_eye_between_light_and_surface() {
        let m = Material::new();
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eyev = Tuple::vector(0.0, 0.0, -1.0);
        let normalv = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::new(Tuple::point(0.0, 0.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        let result = m.lighting(&light, &position, &eyev, &normalv);
        assert_eq!(result, Tuple::color(1.9, 1.9, 1.9));
    }

    #[test]
    fn lighting_eye_between_light_and_surface_eye_offset_45_degrees() {
        let m = Material::new();
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eyev = Tuple::vector(0.0, SQRT_2/2.0, -SQRT_2/2.0);
        let normalv = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::new(Tuple::point(0.0, 0.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        let result = m.lighting(&light, &position, &eyev, &normalv);
        assert_eq!(result, Tuple::color(1.0, 1.0, 1.0));
    }

    #[test]
    fn lighting_eye_opposite_surface_light_offset_45_degrees() {
        let m = Material::new();
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eyev = Tuple::vector(0.0, 0.0, -1.0);
        let normalv = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::new(Tuple::point(0.0, 10.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        let result = m.lighting(&light, &position, &eyev, &normalv);
        assert!(result.nearly_equals(&Tuple::color(0.9*SQRT_2/2.+0.1, 0.9*SQRT_2/2.+0.1, 0.9*SQRT_2/2.+0.1), 0.00001));
    }

    #[test]
    fn lighting_eye_in_path_of_reflection_vector() {
        let m = Material::new();
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eyev = Tuple::vector(0.0, -SQRT_2/2.0, -SQRT_2/2.0);
        let normalv = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::new(Tuple::point(0.0, 10.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        let result = m.lighting(&light, &position, &eyev, &normalv);
        assert_eq!(result, Tuple::color(0.9*SQRT_2/2.+1., 0.9*SQRT_2/2.+1., 0.9*SQRT_2/2.+1.));
    }

    #[test]
    fn lighting_light_behind_surface() {
        let m = Material::new();
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eyev = Tuple::vector(0.0, 0.0, -1.0);
        let normalv = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::new(Tuple::point(0.0, 0.0, 10.0), Tuple::color(1.0, 1.0, 1.0));
        let result = m.lighting(&light, &position, &eyev, &normalv);
        assert_eq!(result, Tuple::color(0.1, 0.1, 0.1));
    }
}
