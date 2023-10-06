use crate::{Intersection, Material, Matrix, Ray, Shape, Transform, Tuple};

#[derive(Clone, Debug, PartialEq)]
pub struct Plane {
    transform: Matrix,
    material: Material,
    id: u32,
}

impl Plane {
    pub fn new(id: u32) -> Plane {
        Plane {
            transform: Matrix::identity(4),
            material: Material::new(),
            id,
        }
    }
}

impl Transform for Plane {
    fn transform(&self) -> &Matrix {
        &self.transform
    }
    fn transform_mut(&mut self) -> &mut Matrix {
        &mut self.transform
    }
}
impl Shape for Plane {
    fn material(&self) -> &Material {
        &self.material
    }
    fn material_mut(&mut self) -> &mut Material {
        &mut self.material
    }
    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
        if ray.direction.y().abs() < 0.0001 {
            vec![]
        } else {
            let t = -ray.origin.y() / ray.direction.y();
            let shape = Box::new(self.clone());
            vec![Intersection::new(t, shape)]
        }
    }
    fn local_normal_at(&self, _point: &Tuple) -> Tuple {
        Tuple::vector(0.0, 1.0, 0.0)
    }
    fn id(&self) -> u32 {
        self.id
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normal_of_plane_is_constant_everywhere() {
        let plane = Plane::new(0);
        let n1 = plane.local_normal_at(&Tuple::point(0.0, 0.0, 0.0));
        let n2 = plane.local_normal_at(&Tuple::point(10.0, 0.0, -10.0));
        let n3 = plane.local_normal_at(&Tuple::point(-5.0, 0.0, 150.0));
        assert_eq!(n1, Tuple::vector(0.0, 1.0, 0.0));
        assert_eq!(n2, Tuple::vector(0.0, 1.0, 0.0));
        assert_eq!(n3, Tuple::vector(0.0, 1.0, 0.0));
    }

    #[test]
    fn intersect_with_ray_parallel_to_plane() {
        let plane = Plane::new(0);
        let ray = Ray::new(Tuple::point(0.0, 10.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = plane.local_intersect(&ray);
        assert_eq!(xs.len(), 0);
    }

    #[test]
    fn intersect_with_coplanar_ray() {
        let plane = Plane::new(0);
        let ray = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = plane.local_intersect(&ray);
        assert_eq!(xs.len(), 0);
    }

    #[test]
    fn ray_intersecting_plane_from_above() {
        let plane = Plane::new(0);
        let ray = Ray::new(Tuple::point(0.0, 1.0, 0.0), Tuple::vector(0.0, -1.0, 0.0));
        let intersections = plane.local_intersect(&ray);
        assert_eq!(intersections.len(), 1);
        assert_eq!(intersections[0].t, 1.0);
        assert!(<dyn Shape>::eq(&plane, &intersections[0].object));
    }

    #[test]
    fn ray_intersecting_plane_from_below() {
        let plane = Plane::new(0);
        let ray = Ray::new(Tuple::point(0.0, -1.0, 0.0), Tuple::vector(0.0, 1.0, 0.0));
        let intersections = plane.local_intersect(&ray);
        assert_eq!(intersections.len(), 1);
        assert_eq!(intersections[0].t, 1.0);
        assert!(<dyn Shape>::eq(&plane, &intersections[0].object));
    }
}
