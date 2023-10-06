use crate::{
    Computations, Intersection, Light, Matrix, Ray, Shape, SolidPattern, Sphere, Transform, Tuple, Counter,
};

pub struct World {
    pub objects: Vec<Box<dyn Shape>>,
    pub lights: Vec<Light>,
    pub counter: Counter,
}

impl World {
    pub fn new() -> World {
        World {
            objects: vec![],
            lights: vec![],
            counter: Counter::new(),
        }
    }

    pub fn default() -> World {
        let mut counter = Counter::new();
        let mut s1 = Sphere::new(counter.next());
        s1.material_mut().pattern = Box::new(SolidPattern::new(Tuple::color(0.8, 1.0, 0.6)));
        s1.material_mut().diffuse = 0.7;
        s1.material_mut().specular = 0.2;
        let mut s2 = Sphere::new(counter.next());
        *s2.transform_mut() = Matrix::scaling(0.5, 0.5, 0.5);
        let light = Light::new(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );
        World {
            objects: vec![Box::new(s1), Box::new(s2)],
            lights: vec![light],
            counter,
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

    pub fn shade_hit(&self, comps: &Computations, remaining: usize) -> Tuple {
        self.lights
            .iter()
            .fold(Tuple::color(0.0, 0.0, 0.0), |sum, light| {
                let shadowed = self.is_shadowed(light, &comps.over_point);
                let surface = comps.object.material().lighting(
                    comps.object.as_ref(),
                    light,
                    &comps.over_point,
                    &comps.eye_vector,
                    &comps.normal_vector,
                    shadowed,
                );
                let reflected = self.reflected_color(comps, remaining);
                let refracted = self.refracted_color(comps, remaining);
                if comps.object.material().reflective > 0.0 && comps.object.material().transparency > 0.0 {
                    let reflectance = comps.schlick();
                    sum + surface + reflected * reflectance + refracted * (1.0 - reflectance)
                } else {
                    sum + surface + reflected + refracted
                }
            })
    }

    pub fn color_at(&self, ray: &Ray, remaining: usize) -> Tuple {
        let intersections = self.intersect(ray);
        if let Some(hit) = Intersection::hit(&intersections) {
            let comps = hit.prepare_computations(ray, Some(&intersections));
            self.shade_hit(&comps, remaining)
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

    pub fn reflected_color(&self, comps: &Computations, remaining: usize) -> Tuple {
        if comps.object.material().reflective == 0.0 || remaining < 1 {
            return Tuple::color(0., 0., 0.);
        }
        let reflect_ray = Ray::new(comps.over_point.clone(), comps.reflect_vector.clone());
        let color = self.color_at(&reflect_ray, remaining - 1);
        color * comps.object.material().reflective
    }

    pub fn refracted_color(&self, comps: &Computations, remaining: usize) -> Tuple {
        if comps.object.material().transparency == 0.0 || remaining < 1 {
            return Tuple::color(0., 0., 0.);
        }
        let n_ratio = comps.n1 / comps.n2;
        let cos_i = comps.eye_vector.dot(&comps.normal_vector);
        let sin2_t = n_ratio.powi(2) * (1.0 - cos_i.powi(2));
        if sin2_t > 1.0 {
            return Tuple::color(0., 0., 0.);
        }
        let cos_t = (1.0 - sin2_t).sqrt();
        let direction = &comps.normal_vector * (n_ratio * cos_i - cos_t) - &comps.eye_vector * n_ratio;
        let refract_ray = Ray::new(comps.under_point.clone(), direction);
        self.color_at(&refract_ray, remaining - 1) * comps.object.material().transparency
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Material, Plane, TestPattern};
    use std::f64::consts::SQRT_2;
    use approx::assert_relative_eq;

    #[test]
    fn creating_world() {
        let world = World::new();
        assert_eq!(world.objects.len(), 0);
        assert_eq!(world.lights.len(), 0);
    }

    #[test]
    fn default_world() {
        let light = Light::new(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );
        let mut material = Material::new();
        material.pattern = Box::new(SolidPattern::new(Tuple::color(0.8, 1.0, 0.6)));
        material.diffuse = 0.7;
        material.specular = 0.2;
        let transform = Matrix::scaling(0.5, 0.5, 0.5);
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
        let comps = intersection.prepare_computations(&ray, None);
        let color = world.shade_hit(&comps, 5);
        assert_relative_eq!(color, Tuple::color(0.38066, 0.47583, 0.2855), epsilon = 1e-5f64);
    }

    #[test]
    fn shading_intersection_inside() {
        let mut world = World::default();
        world.lights = vec![Light::new(
            Tuple::point(0., 0.25, 0.),
            Tuple::color(1., 1., 1.),
        )];
        let ray = Ray::new(Tuple::point(0., 0., 0.), Tuple::vector(0., 0., 1.));
        let shape = world.objects[1].clone();
        let intersection = Intersection::new(0.5, shape);
        let comps = intersection.prepare_computations(&ray, None);
        let color = world.shade_hit(&comps, 5);
        assert_relative_eq!(color, Tuple::color(0.90498, 0.90498, 0.90498), epsilon = 1e-5f64);
    }

    #[test]
    fn color_when_ray_misses() {
        let world = World::default();
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 1., 0.));
        let color = world.color_at(&ray, 5);
        assert_eq!(color, Tuple::color(0., 0., 0.))
    }

    #[test]
    fn color_when_ray_hits() {
        let world = World::default();
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let color = world.color_at(&ray, 5);
        assert_relative_eq!(color, Tuple::color(0.38066, 0.47583, 0.2855), epsilon = 1e-5f64);
    }

    #[test]
    fn color_intersection_behind_ray() {
        let mut world = World::default();
        world.objects[0].material_mut().ambient = 1.;
        world.objects[1].material_mut().ambient = 1.;
        let ray = Ray::new(Tuple::point(0., 0., 0.75), Tuple::vector(0., 0., -1.));
        let color = world.color_at(&ray, 5);
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
        world.lights = vec![Light::new(
            Tuple::point(0., 0., -10.),
            Tuple::color(1., 1., 1.),
        )];
        let shape_a = Sphere::new(world.counter.next());
        let mut shape_b = Sphere::new(world.counter.next());
        shape_b.transform_mut().translate(0., 0., 10.);
        world.objects = vec![Box::new(shape_a), Box::new(shape_b)];
        let ray = Ray::new(Tuple::point(0., 0., 5.), Tuple::vector(0., 0., 1.));
        let intersection = Intersection::new(4., world.objects[1].clone());
        let comps = intersection.prepare_computations(&ray, None);
        let color = world.shade_hit(&comps, 5);
        assert_relative_eq!(color, Tuple::color(0.1, 0.1, 0.1), epsilon = 1e-5f64);
    }

    #[test]
    fn reflected_color_for_nonreflective_material() {
        let mut world = World::default();
        let ray = Ray::new(Tuple::point(0., 0., 0.), Tuple::vector(0., 0., 1.));
        world.objects[1].material_mut().ambient = 1.;
        let intersection = Intersection::new(1., world.objects[1].clone());
        let comps = intersection.prepare_computations(&ray, None);
        let color = world.reflected_color(&comps, 5);
        assert_eq!(color, Tuple::color(0., 0., 0.));
    }

    #[test]
    fn reflected_color_for_reflective_material() {
        let mut world = World::default();
        let mut plane = Plane::new(world.counter.next());
        plane.material_mut().reflective = 0.5;
        *plane.transform_mut() = Matrix::translation(0., -1., 0.);
        let shape = Box::new(plane);
        world.objects.push(shape.clone());
        let ray = Ray::new(
            Tuple::point(0., 0., -3.),
            Tuple::vector(0., -SQRT_2 / 2., SQRT_2 / 2.),
        );
        let intersection = Intersection::new(SQRT_2, shape);
        let comps = intersection.prepare_computations(&ray, None);
        let color = world.reflected_color(&comps, 5);
        assert_relative_eq!(color, Tuple::color(0.19032, 0.2379, 0.14274), epsilon = 1e-4f64);
    }

    #[test]
    fn shade_hit_with_reflective_material() {
        let mut world = World::default();
        let mut plane = Plane::new(world.counter.next());
        plane.material_mut().reflective = 0.5;
        *plane.transform_mut() = Matrix::translation(0., -1., 0.);
        let shape = Box::new(plane);
        world.objects.push(shape.clone());
        let ray = Ray::new(
            Tuple::point(0., 0., -3.),
            Tuple::vector(0., -SQRT_2 / 2., SQRT_2 / 2.),
        );
        let intersection = Intersection::new(SQRT_2, shape);
        let comps = intersection.prepare_computations(&ray, None);
        let color = world.shade_hit(&comps, 5);
        assert_relative_eq!(color, Tuple::color(0.87677, 0.92436, 0.82918), epsilon = 1e-4f64);
    }

    #[test]
    fn color_at_with_mutually_reflective_surfaces() {
        let mut world = World::default();
        world.lights[0].position = Tuple::point(0., 0., 0.);
        let mut lower = Plane::new(world.counter.next());
        lower.material_mut().reflective = 1.;
        *lower.transform_mut() = Matrix::translation(0., -1., 0.);
        world.objects.push(Box::new(lower));
        let mut upper = Plane::new(world.counter.next());
        upper.material_mut().reflective = 1.;
        *upper.transform_mut() = Matrix::translation(0., 1., 0.);
        world.objects.push(Box::new(upper));
        let ray = Ray::new(Tuple::point(0., 0., 0.), Tuple::vector(0., 1., 0.));
        let color = world.color_at(&ray, 5);
        assert_relative_eq!(color, Tuple::color(1.9, 1.9, 1.9), epsilon = 1e-5f64);
    }

    #[test]
    fn reflected_color_at_maximum_recursive_depth() {
        let mut world = World::default();
        let mut plane = Plane::new(world.counter.next());
        plane.material_mut().reflective = 0.5;
        *plane.transform_mut() = Matrix::translation(0., -1., 0.);
        let shape = Box::new(plane);
        world.objects.push(shape.clone());
        let ray = Ray::new(
            Tuple::point(0., 0., -3.),
            Tuple::vector(0., -SQRT_2 / 2., SQRT_2 / 2.),
        );
        let intersection = Intersection::new(SQRT_2, shape);
        let comps = intersection.prepare_computations(&ray, None);
        let color = world.reflected_color(&comps, 0);
        assert_eq!(color, Tuple::color(0., 0., 0.));
    }

    #[test]
    fn refracted_color_with_opaque_surface() {
        let world = World::default();
        let shape = &world.objects[0];
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let intersections = vec![
            Intersection::new(4., shape.clone()),
            Intersection::new(6., shape.clone()),
        ];
        let comps = intersections[0].prepare_computations(&ray, Some(&intersections));
        let c = world.refracted_color(&comps, 5);
        assert_eq!(c, Tuple::color(0., 0., 0.));
    }

    #[test]
    fn refracted_color_at_maximum_recursive_depth() {
        let mut world = World::default();
        world.objects[0].material_mut().transparency = 1.;
        world.objects[0].material_mut().refractive_index = 1.5;
        let shape = &world.objects[0];
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let intersections = vec![
            Intersection::new(4., shape.clone()),
            Intersection::new(6., shape.clone()),
        ];
        let comps = intersections[0].prepare_computations(&ray, Some(&intersections));
        let c = world.refracted_color(&comps, 0);
        assert_eq!(c, Tuple::color(0., 0., 0.));
    }

    #[test]
    fn refracted_color_under_total_internal_reflection() {
        let mut world = World::default();
        world.objects[0].material_mut().transparency = 1.;
        world.objects[0].material_mut().refractive_index = 1.5;
        let shape = &world.objects[0];
        let ray = Ray::new(Tuple::point(0., 0., SQRT_2 / 2.), Tuple::vector(0., 1., 0.));
        let intersections = vec![
            Intersection::new(-SQRT_2 / 2., shape.clone()),
            Intersection::new(SQRT_2 / 2., shape.clone()),
        ];
        let comps = intersections[1].prepare_computations(&ray, Some(&intersections));
        let c = world.refracted_color(&comps, 5);
        assert_eq!(c, Tuple::color(0., 0., 0.));
    }

    #[test]
    fn refracted_color_with_refracted_ray() {
        let mut world = World::default();
        world.objects[0].material_mut().ambient = 1.;
        world.objects[0].material_mut().pattern = Box::new(TestPattern::new());
        world.objects[1].material_mut().transparency = 1.;
        world.objects[1].material_mut().refractive_index = 1.5;
        let shape_a = &world.objects[0];
        let shape_b = &world.objects[1];
        let ray = Ray::new(Tuple::point(0., 0., 0.1), Tuple::vector(0., 1., 0.));
        let intersections = vec![
            Intersection::new(-0.9899, shape_a.clone()),
            Intersection::new(-0.4899, shape_b.clone()),
            Intersection::new(0.4899, shape_b.clone()),
            Intersection::new(0.9899, shape_a.clone()),
        ];
        let comps = intersections[2].prepare_computations(&ray, Some(&intersections));
        let c = world.refracted_color(&comps, 5);
        assert_relative_eq!(c, Tuple::color(0., 0.99888, 0.04725), epsilon = 1e-4f64);
    }

    #[test]
    fn shade_hit_with_transparent_material() {
        let mut world = World::default();
        let mut floor = Box::new(Plane::new(world.counter.next()));
        floor.material_mut().transparency = 0.5;
        floor.material_mut().refractive_index = 1.5;
        *floor.transform_mut() = Matrix::translation(0., -1., 0.);
        world.objects.push(floor.clone());
        let mut ball = Box::new(Sphere::new(world.counter.next()));
        ball.material_mut().pattern = Box::new(SolidPattern::new(Tuple::color(1., 0., 0.)));
        ball.material_mut().ambient = 0.5;
        *ball.transform_mut() = Matrix::translation(0., -3.5, -0.5);
        world.objects.push(ball);
        let ray = Ray::new(Tuple::point(0., 0., -3.), Tuple::vector(0., -SQRT_2 / 2., SQRT_2 / 2.));
        let intersections = vec![
            Intersection::new(SQRT_2, floor),
        ];
        let comps = intersections[0].prepare_computations(&ray, Some(&intersections));
        let color = world.shade_hit(&comps, 5);
        assert_relative_eq!(color, Tuple::color(0.93642, 0.68642, 0.68642), epsilon = 1e-5f64);
    }

    #[test]
    fn shade_hit_with_reflective_transparent_material() {
        let mut world = World::default();
        let mut floor = Box::new(Plane::new(world.counter.next()));
        floor.material_mut().reflective = 0.5;
        floor.material_mut().transparency = 0.5;
        floor.material_mut().refractive_index = 1.5;
        *floor.transform_mut() = Matrix::translation(0., -1., 0.);
        world.objects.push(floor.clone());
        let mut ball = Box::new(Sphere::new(world.counter.next()));
        ball.material_mut().pattern = Box::new(SolidPattern::new(Tuple::color(1., 0., 0.)));
        ball.material_mut().ambient = 0.5;
        *ball.transform_mut() = Matrix::translation(0., -3.5, -0.5);
        world.objects.push(ball);
        let ray = Ray::new(Tuple::point(0., 0., -3.), Tuple::vector(0., -SQRT_2 / 2., SQRT_2 / 2.));
        let intersections = vec![
            Intersection::new(SQRT_2, floor),
        ];
        let comps = intersections[0].prepare_computations(&ray, Some(&intersections));
        let color = world.shade_hit(&comps, 5);
        assert_relative_eq!(color, Tuple::color(0.93391, 0.69643, 0.69243), epsilon = 1e-5f64);
    }
}
