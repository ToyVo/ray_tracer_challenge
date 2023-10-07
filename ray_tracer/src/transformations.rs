use nalgebra_glm::{mat4, translate, DMat4, DVec3};

pub fn view_transform(from: DVec3, to: DVec3, up: DVec3) -> DMat4 {
    let forward = (to - from).normalize();
    let left = forward.cross(&up.normalize());
    let true_up = left.cross(&forward);
    let orientation = mat4(
        left.x, left.y, left.z, 0., true_up.x, true_up.y, true_up.z, 0., -forward.x, -forward.y,
        -forward.z, 0., 0., 0., 0., 1.,
    );
    translate(&orientation, &-from)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use nalgebra_glm::{scale, vec3, DMat4};

    #[test]
    fn transformation_matrix_default_orientation() {
        let from = vec3(0., 0., 0.);
        let to = vec3(0., 0., -1.);
        let up = vec3(0., 1., 0.);
        let transform = view_transform(from, to, up);
        assert_eq!(transform, DMat4::identity());
    }

    #[test]
    fn transformation_matrix_looking_positive_z() {
        let from = vec3(0., 0., 0.);
        let to = vec3(0., 0., 1.);
        let up = vec3(0., 1., 0.);
        let transform = view_transform(from, to, up);
        assert_eq!(transform, scale(&DMat4::identity(), &vec3(-1., 1., -1.)));
    }

    #[test]
    fn view_transformation_moves_the_world() {
        let from = vec3(0., 0., 8.);
        let to = vec3(0., 0., 0.);
        let up = vec3(0., 1., 0.);
        let transform = view_transform(from, to, up);
        assert_eq!(transform, translate(&DMat4::identity(), &vec3(0., 0., -8.)));
    }

    #[test]
    fn arbitrary_view_transformation() {
        let from = vec3(1., 3., 2.);
        let to = vec3(4., -2., 8.);
        let up = vec3(1., 1., 0.);
        let transform = view_transform(from, to, up);
        assert_relative_eq!(
            transform,
            mat4(
                -0.50709, 0.50709, 0.67612, -2.36643, 0.76772, 0.60609, 0.12122, -2.82843,
                -0.35857, 0.59761, -0.71714, 0., 0., 0., 0., 1.,
            ),
            epsilon = 1e-5f64
        );
    }
}
