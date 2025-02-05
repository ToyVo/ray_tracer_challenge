use crate::Tuple;

#[derive(Debug, PartialEq, Clone)]
pub struct Light {
    pub intensity: Tuple,
    pub position: Tuple,
}

impl Light {
    pub fn new(position: Tuple, intensity: Tuple) -> Light {
        Light {
            intensity,
            position,
        }
    }
}

// Feature: Lights
#[cfg(test)]
mod tests {
    use super::*;

    // Scenario: A point light has a position and intensity
    //   Given intensity ← color(1, 1, 1)
    //     And position ← point(0, 0, 0)
    //   When light ← point_light(position, intensity)
    //   Then light.position = position
    //     And light.intensity = intensity
    #[test]
    fn point_light_has_position_and_intensity() {
        let intensity = Tuple::color(1.0, 1.0, 1.0);
        let position = Tuple::point(0.0, 0.0, 0.0);
        let light = Light {
            intensity: intensity.clone(),
            position: position.clone(),
        };
        assert_eq!(light.intensity, intensity);
        assert_eq!(light.position, position);
    }
}
