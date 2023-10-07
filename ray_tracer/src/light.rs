use nalgebra_glm::{DVec3, DVec4};

#[derive(Debug, PartialEq, Clone)]
pub struct Light {
    pub intensity: DVec3,
    pub position: DVec4,
}

impl Light {
    pub fn new(position: DVec4, intensity: DVec3) -> Light {
        Light {
            intensity,
            position,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra_glm::{vec3, vec4};

    #[test]
    fn point_light_has_position_and_intensity() {
        let intensity = vec3(1.0, 1.0, 1.0);
        let position = vec4(0.0, 0.0, 0.0, 1.);
        let light = Light {
            intensity,
            position,
        };
        assert_eq!(light.intensity, intensity);
        assert_eq!(light.position, position);
    }
}
