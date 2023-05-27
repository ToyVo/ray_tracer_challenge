use crate::Tuple;

#[derive(Debug, PartialEq, Clone)]
pub struct Light {
    pub intensity: Tuple,
    pub position: Tuple,
}

impl Light {
    pub fn new(position: Tuple, intensity: Tuple) -> Light {
        Light { intensity, position }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_light_has_position_and_intensity() {
        let intensity = Tuple::color(1.0, 1.0, 1.0);
        let position = Tuple::point(0.0, 0.0, 0.0);
        let light = Light { intensity: intensity.clone(), position: position.clone() };
        assert_eq!(light.intensity, intensity);
        assert_eq!(light.position, position);
    }
}
