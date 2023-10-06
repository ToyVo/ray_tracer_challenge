use crate::{Matrix, Tuple};

pub fn view_transform(from: Tuple, to: Tuple, up: Tuple) -> Matrix {
    let forward = (&to - &from).normalize();
    let left = forward.cross(&up.normalize());
    let true_up = left.cross(&forward);
    let orientation = Matrix::from_vec(
        4,
        4,
        vec![
            left.x(),
            true_up.x(),
            -forward.x(),
            0.,
            left.y(),
            true_up.y(),
            -forward.y(),
            0.,
            left.z(),
            true_up.z(),
            -forward.z(),
            0.,
            0.,
            0.,
            0.,
            1.,
        ],
    );
    orientation * Matrix::translation(-from.x(), -from.y(), -from.z())
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn transformation_matrix_default_orientation() {
        let from = Tuple::point(0., 0., 0.);
        let to = Tuple::point(0., 0., -1.);
        let up = Tuple::vector(0., 1., 0.);
        let transform = view_transform(from, to, up);
        assert_eq!(transform, Matrix::identity(4));
    }

    #[test]
    fn transformation_matrix_looking_positive_z() {
        let from = Tuple::point(0., 0., 0.);
        let to = Tuple::point(0., 0., 1.);
        let up = Tuple::vector(0., 1., 0.);
        let transform = view_transform(from, to, up);
        assert_eq!(transform, Matrix::scaling(-1., 1., -1.));
    }

    #[test]
    fn view_transformation_moves_the_world() {
        let from = Tuple::point(0., 0., 8.);
        let to = Tuple::point(0., 0., 0.);
        let up = Tuple::vector(0., 1., 0.);
        let transform = view_transform(from, to, up);
        assert_eq!(transform, Matrix::translation(0., 0., -8.));
    }

    #[test]
    fn arbitrary_view_transformation() {
        let from = Tuple::point(1., 3., 2.);
        let to = Tuple::point(4., -2., 8.);
        let up = Tuple::vector(1., 1., 0.);
        let transform = view_transform(from, to, up);
        assert_relative_eq!(
            transform,
            Matrix::from_vec(
                4,
                4,
                vec![
                    -0.50709, 0.76772, -0.35857, 0., 0.50709, 0.60609, 0.59761, 0., 0.67612,
                    0.12122, -0.71714, 0., -2.36643, -2.82843, 0., 1.,
                ]
            ),
            epsilon = 1e-5f64
        );
    }
}
