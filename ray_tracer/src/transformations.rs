use crate::{Matrix, Tuple};

// function view_transform(from, to, up)
//   forward ← normalize(to - from)
//   upn ← normalize(up)
//   left ← cross(forward, upn)
//   true_up ← cross(left, forward)
//
//   orientation ← matrix( left.x,     left.y,     left.z,    0,
//                          true_up.x,  true_up.y,  true_up.z, 0,
//                         -forward.x, -forward.y, -forward.z, 0,
//                          0,          0,          0,         1)
//
//   return orientation * translation(-from.x, -from.y, -from.z)
// end function
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
    orientation.translate(-from.x(), -from.y(), -from.z())
}

// Feature: Matrix Transformations
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Matrix, Tuple};
    use approx::assert_relative_eq;

    // Scenario: The transformation matrix for the default orientation
    //   Given from ← point(0, 0, 0)
    //     And to ← point(0, 0, -1)
    //     And up ← vector(0, 1, 0)
    //   When t ← view_transform(from, to, up)
    //   Then t = identity_matrix
    #[test]
    fn transformation_matrix_default_orientation() {
        let from = Tuple::point(0., 0., 0.);
        let to = Tuple::point(0., 0., -1.);
        let up = Tuple::vector(0., 1., 0.);
        let transform = view_transform(from, to, up);
        assert_eq!(transform, Matrix::identity(4));
    }

    // Scenario: A view transformation matrix looking in positive z direction
    //   Given from ← point(0, 0, 0)
    //     And to ← point(0, 0, 1)
    //     And up ← vector(0, 1, 0)
    //   When t ← view_transform(from, to, up)
    //   Then t = scaling(-1, 1, -1)
    #[test]
    fn transformation_matrix_looking_positive_z() {
        let from = Tuple::point(0., 0., 0.);
        let to = Tuple::point(0., 0., 1.);
        let up = Tuple::vector(0., 1., 0.);
        let transform = view_transform(from, to, up);
        assert_eq!(transform, Matrix::scaling(-1., 1., -1.));
    }

    // Scenario: The view transformation moves the world
    //   Given from ← point(0, 0, 8)
    //     And to ← point(0, 0, 0)
    //     And up ← vector(0, 1, 0)
    //   When t ← view_transform(from, to, up)
    //   Then t = translation(0, 0, -8)
    #[test]
    fn view_transformation_moves_the_world() {
        let from = Tuple::point(0., 0., 8.);
        let to = Tuple::point(0., 0., 0.);
        let up = Tuple::vector(0., 1., 0.);
        let transform = view_transform(from, to, up);
        assert_eq!(transform, Matrix::translation(0., 0., -8.));
    }

    // Scenario: An arbitrary view transformation
    //   Given from ← point(1, 3, 2)
    //     And to ← point(4, -2, 8)
    //     And up ← vector(1, 1, 0)
    //   When t ← view_transform(from, to, up)
    //   Then t is the following 4x4 matrix:
    //       | -0.50709 | 0.50709 |  0.67612 | -2.36643 |
    //       |  0.76772 | 0.60609 |  0.12122 | -2.82843 |
    //       | -0.35857 | 0.59761 | -0.71714 |  0.00000 |
    //       |  0.00000 | 0.00000 |  0.00000 |  1.00000 |
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
