use crate::{Computations, Intersection, Light, Matrix, Ray, Sphere, Tuple, Shape, SolidPattern, Transform};

pub struct World {
    pub objects: Vec<Box<dyn Shape>>,
    pub lights: Vec<Light>,
}

impl World {
    pub fn new() -> World {
        World {
            objects: vec![],
            lights: vec![],
        }
    }

    pub fn default() -> World {
        let mut s1 = Sphere::new();
        s1.material_mut().pattern = Box::new(SolidPattern::new(Tuple::color(0.8, 1.0, 0.6)));
        s1.material_mut().diffuse = 0.7;
        s1.material_mut().specular = 0.2;
        let mut s2 = Sphere::new();
        *s2.transform_mut() = Matrix::scaling(0.5, 0.5, 0.5);
        let light = Light::new(Tuple::point(-10.0, 10.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        World {
            objects: vec![Box::new(s1), Box::new(s2)],
            lights: vec![light],
        }
    }

    pub fn intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let mut intersections: Vec<Intersection> = vec![];
        for object in &self.objects {
            let mut object_intersections = object.intersect(ray);
            intersections.append(&mut object_intersections);
        }
        intersections.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
        intersections
    }

    pub fn shade_hit(&self, comps: &Computations) -> Tuple {
        self.lights.iter().fold(Tuple::color(0.0, 0.0, 0.0), |sum, light| {
            let shadowed = self.is_shadowed(light, &comps.over_point);
            sum + comps.object.material().lighting(comps.object.as_ref(), light, &comps.over_point, &comps.eye_vector, &comps.normal_vector, shadowed)
        })
    }

    pub fn color_at(&self, ray: &Ray) -> Tuple {
        let intersections = self.intersect(ray);
        if let Some(hit) = Intersection::hit(&intersections) {
            let comps = hit.prepare_computations(ray);
            self.shade_hit(&comps)
        } else {
            Tuple::color(0., 0., 0.)
        }
    }

    pub fn is_shadowed(&self, light: &Light, point: &Tuple) -> bool {
        let vector = &light.position - point;
        let distance = vector.magnitude();
        let direction = vector.normalize();
        let ray = Ray::new(point.clone(), direction);
        let intersections = self.intersect(&ray);
        if let Some(hit) = Intersection::hit(&intersections) {
            if hit.t < distance {
                return true;
            }
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Material;

    #[test]
    fn creating_world() {
        let world = World::new();
        assert_eq!(world.objects.len(), 0);
        assert_eq!(world.lights.len(), 0);
    }

    #[test]
    fn default_world() {
        let light = Light::new(Tuple::point(-10.0, 10.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        let mut material = Material::new();
        material.pattern = Box::new(SolidPattern::new(Tuple::color(0.8, 1.0, 0.6)));
        material.diffuse = 0.7;
        material.specular = 0.2;
        let transform= Matrix::scaling(0.5, 0.5, 0.5);
        let world = World::default();
        assert_eq!(world.lights.len(), 1);
        assert_eq!(world.lights[0], light);
        assert_eq!(world.objects.len(), 2);
        assert_eq!(world.objects[0].material(), &material);
        assert_eq!(world.objects[1].transform(), &transform);
    }

    #[test]
    fn intersect_world_with_ray() {
        let world = World::default();
        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let intersections = world.intersect(&ray);
        assert_eq!(intersections.len(), 4);
        assert_eq!(intersections[0].t, 4.0);
        assert_eq!(intersections[1].t, 4.5);
        assert_eq!(intersections[2].t, 5.5);
        assert_eq!(intersections[3].t, 6.0);
    }

    #[test]
    fn shading_intersection() {
        let world = World::default();
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let shape = world.objects[0].clone();
        let intersection = Intersection::new(4., shape);
        let comps = intersection.prepare_computations(&ray);
        let color = world.shade_hit(&comps);
        assert!(color.nearly_equals(&Tuple::color(0.38066, 0.47583, 0.2855), 1e-3f64));
    }

    #[test]
    fn shading_intersection_inside() {
        let mut world = World::default();
        world.lights = vec![Light::new(Tuple::point(0., 0.25, 0.), Tuple::color(1., 1., 1.))];
        let ray = Ray::new(Tuple::point(0., 0., 0.), Tuple::vector(0., 0., 1.));
        let shape = world.objects[1].clone();
        let intersection = Intersection::new(0.5, shape);
        let comps = intersection.prepare_computations(&ray);
        let color = world.shade_hit(&comps);
        assert!(color.nearly_equals(&Tuple::color(0.90498, 0.90498, 0.90498), 1e-3f64));
    }

    #[test]
    fn color_when_ray_misses() {
        let world = World::default();
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 1., 0.));
        let color = world.color_at(&ray);
        assert_eq!(color, Tuple::color(0., 0., 0.))
    }

    #[test]
    fn color_when_ray_hits() {
        let world = World::default();
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let color = world.color_at(&ray);
        assert!(color.nearly_equals(&Tuple::color(0.38066, 0.47583, 0.2855), 1e-3f64));
    }

    #[test]
    fn color_intersection_behind_ray() {
        let mut world = World::default();
        Box::get_mut(&mut world.objects[0]).unwrap().material_mut().ambient = 1.;
        Box::get_mut(&mut world.objects[1]).unwrap().material_mut().ambient = 1.;
        let ray = Ray::new(Tuple::point(0., 0., 0.75), Tuple::vector(0., 0., -1.));
        let color = world.color_at(&ray);
        assert_eq!(&color, world.objects[1].material().pattern.colors()[0])
    }

    #[test]
    fn no_shadow_when_nothing_collinear_with_point_and_light() {
        let world = World::default();
        let point = Tuple::point(0., 10., 0.);
        assert!(!world.is_shadowed(&world.lights[0], &point));
    }

    #[test]
    fn shadow_when_object_between_point_and_light() {
        let world = World::default();
        let point = Tuple::point(10., -10., 10.);
        assert!(world.is_shadowed(&world.lights[0], &point));
    }

    #[test]
    fn no_shadow_when_object_behind_light() {
        let world = World::default();
        let point = Tuple::point(-20., 20., -20.);
        assert!(!world.is_shadowed(&world.lights[0], &point));
    }

    #[test]
    fn no_shadow_when_object_behind_point() {
        let world = World::default();
        let point = Tuple::point(-2., 2., -2.);
        assert!(!world.is_shadowed(&world.lights[0], &point));
    }

    #[test]
    fn shade_hit_given_intersection_in_shadow() {
        let mut world = World::default();
        world.lights = vec![Light::new(Tuple::point(0., 0., -10.), Tuple::color(1., 1., 1.))];
        let shape_a = Sphere::new();
        let mut shape_b = Sphere::new();
        shape_b.transform_mut().translate(0., 0., 10.);
        world.objects = vec![Box::new(shape_a), Box::new(shape_b)];
        let ray = Ray::new(Tuple::point(0., 0., 5.), Tuple::vector(0., 0., 1.));
        let intersection = Intersection::new(4., world.objects[1].clone());
        let comps = intersection.prepare_computations(&ray);
        let color = world.shade_hit(&comps);
        assert!(color.nearly_equals(&Tuple::color(0.1, 0.1, 0.1), 1e-3f64));
    }
}
