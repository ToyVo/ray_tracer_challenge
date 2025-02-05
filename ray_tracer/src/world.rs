use crate::{
    Computations, Counter, Intersection, Light, Matrix, Ray, Shape, SolidPattern, Sphere,
    Transform, Tuple,
};

pub struct World {
    pub objects: Vec<Box<dyn Shape>>,
    pub lights: Vec<Light>,
    pub counter: Counter,
}

impl Default for World {
    fn default() -> World {
        let mut counter = Counter::default();
        let mut s1 = Sphere::new(counter.increment());
        s1.material_mut().pattern = Box::new(SolidPattern::new(Tuple::color(0.8, 1.0, 0.6)));
        s1.material_mut().diffuse = 0.7;
        s1.material_mut().specular = 0.2;
        let mut s2 = Sphere::new(counter.increment());
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
}

impl World {
    pub fn new() -> World {
        World {
            objects: vec![],
            lights: vec![],
            counter: Counter::default(),
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

    // function shade_hit(world, comps, remaining)
    //   shadowed ← is_shadowed(world, comps.over_point)
    //
    //   surface ← lighting(comps.object.material,
    //                      comps.object,
    //                      world.light,
    //                      comps.over_point, comps.eyev, comps.normalv,
    //                      shadowed)
    //
    //   reflected ← reflected_color(world, comps, remaining)
    //   refracted ← refracted_color(world, comps, remaining)
    //
    //   material ← comps.object.material
    //   if material.reflective > 0 && material.transparency > 0
    //     reflectance ← schlick(comps)
    //     return surface + reflected * reflectance +
    //                      refracted * (1 - reflectance)
    //   else
    //     return surface + reflected + refracted
    //   end
    // end function
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
                if comps.object.material().reflective > 0.0
                    && comps.object.material().transparency > 0.0
                {
                    let reflectance = comps.schlick();
                    sum + surface + reflected * reflectance + refracted * (1.0 - reflectance)
                } else {
                    sum + surface + reflected + refracted
                }
            })
    }

    // function color_at(world, ray, remaining)
    //   # ...
    //   color ← shade_hit(world, comps, remaining)
    //   # ...
    // end function
    pub fn color_at(&self, ray: &Ray, remaining: usize) -> Tuple {
        let intersections = self.intersect(ray);
        if let Some(hit) = Intersection::hit(&intersections) {
            let comps = hit.prepare_computations(ray, Some(&intersections));
            self.shade_hit(&comps, remaining)
        } else {
            Tuple::color(0., 0., 0.)
        }
    }

    // function is_shadowed(world, point)
    //   v ← world.light.position - point
    //   distance ← magnitude(v)
    //   direction ← normalize(v)
    //
    //   r ← ray(point, direction)
    //   intersections ← intersect_world(world, r)
    //
    //   h ← hit(intersections)
    //   if h is present and h.t < distance
    //     return true
    //   else
    //     return false
    //   end if
    // end function
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

    // function reflected_color(world, comps, remaining)
    //   if remaining <= 0
    //     return color(0, 0, 0)
    //   end if
    //   if comps.object.material.reflective = 0
    //     return color(0, 0, 0)
    //   end if
    //
    //   reflect_ray ← ray(comps.over_point, comps.reflectv)
    //   color ← color_at(world, reflect_ray, remaining - 1)
    //
    //   return color * comps.object.material.reflective
    // end function
    pub fn reflected_color(&self, comps: &Computations, remaining: usize) -> Tuple {
        if comps.object.material().reflective == 0.0 || remaining < 1 {
            return Tuple::color(0., 0., 0.);
        }
        let reflect_ray = Ray::new(comps.over_point.clone(), comps.reflect_vector.clone());
        let color = self.color_at(&reflect_ray, remaining - 1);
        color * comps.object.material().reflective
    }

    // function refracted_color(world, comps, remaining)
    //   if comps.object.material.transparency = 0
    //     return color(0, 0, 0)
    //   end if
    //
    //   return color(1, 1, 1)
    // end function
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
        let direction =
            &comps.normal_vector * (n_ratio * cos_i - cos_t) - &comps.eye_vector * n_ratio;
        let refract_ray = Ray::new(comps.under_point.clone(), direction);
        self.color_at(&refract_ray, remaining - 1) * comps.object.material().transparency
    }
}

// Feature: World
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Material, Plane, TestPattern};
    use approx::assert_relative_eq;
    use std::f64::consts::SQRT_2;

    // Scenario: Creating a world
    //   Given w ← world()
    //   Then w contains no objects
    //     And w has no light source
    #[test]
    fn creating_world() {
        let world = World::new();
        assert_eq!(world.objects.len(), 0);
        assert_eq!(world.lights.len(), 0);
    }

    // Scenario: The default world
    //   Given light ← point_light(point(-10, 10, -10), color(1, 1, 1))
    //     And s1 ← sphere() with:
    //       | material.color     | (0.8, 1.0, 0.6)        |
    //       | material.diffuse   | 0.7                    |
    //       | material.specular  | 0.2                    |
    //     And s2 ← sphere() with:
    //       | transform | scaling(0.5, 0.5, 0.5) |
    //   When w ← default_world()
    //   Then w.light = light
    //     And w contains s1
    //     And w contains s2
    #[test]
    fn default_world() {
        let light = Light::new(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );
        let mut material = Material::default_with_pattern(Box::new(SolidPattern::new(
            Tuple::color(0.8, 1.0, 0.6),
        )));
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

    // Scenario: Intersect a world with a ray
    //   Given w ← default_world()
    //     And r ← ray(point(0, 0, -5), vector(0, 0, 1))
    //   When xs ← intersect_world(w, r)
    //   Then xs.count = 4
    //     And xs[0].t = 4
    //     And xs[1].t = 4.5
    //     And xs[2].t = 5.5
    //     And xs[3].t = 6
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

    // Scenario: Shading an intersection
    //   Given w ← default_world()
    //     And r ← ray(point(0, 0, -5), vector(0, 0, 1))
    //     And shape ← the first object in w
    //     And i ← intersection(4, shape)
    //   When comps ← prepare_computations(i, r)
    //     And c ← shade_hit(w, comps)
    //   Then c = color(0.38066, 0.47583, 0.2855)
    #[test]
    fn shading_intersection() {
        let world = World::default();
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let shape = world.objects[0].clone();
        let intersection = Intersection::new(4., shape);
        let comps = intersection.prepare_computations(&ray, None);
        let color = world.shade_hit(&comps, 5);
        assert_relative_eq!(
            color,
            Tuple::color(0.38066, 0.47583, 0.2855)
        );
    }

    // Scenario: Shading an intersection from the inside
    //   Given w ← default_world()
    //     And w.light ← point_light(point(0, 0.25, 0), color(1, 1, 1))
    //     And r ← ray(point(0, 0, 0), vector(0, 0, 1))
    //     And shape ← the second object in w
    //     And i ← intersection(0.5, shape)
    //   When comps ← prepare_computations(i, r)
    //     And c ← shade_hit(w, comps)
    //   Then c = color(0.90498, 0.90498, 0.90498)
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
        assert_relative_eq!(
            color,
            Tuple::color(0.90498, 0.90498, 0.90498)
        );
    }

    // Scenario: The color when a ray misses
    //   Given w ← default_world()
    //     And r ← ray(point(0, 0, -5), vector(0, 1, 0))
    //   When c ← color_at(w, r)
    //   Then c = color(0, 0, 0)
    #[test]
    fn color_when_ray_misses() {
        let world = World::default();
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 1., 0.));
        let color = world.color_at(&ray, 5);
        assert_eq!(color, Tuple::color(0., 0., 0.))
    }

    // Scenario: The color when a ray hits
    //   Given w ← default_world()
    //     And r ← ray(point(0, 0, -5), vector(0, 0, 1))
    //   When c ← color_at(w, r)
    //   Then c = color(0.38066, 0.47583, 0.2855)
    #[test]
    fn color_when_ray_hits() {
        let world = World::default();
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let color = world.color_at(&ray, 5);
        assert_relative_eq!(
            color,
            Tuple::color(0.38066, 0.47583, 0.2855)
        );
    }

    // Scenario: The color with an intersection behind the ray
    //   Given w ← default_world()
    //     And outer ← the first object in w
    //     And outer.material.ambient ← 1
    //     And inner ← the second object in w
    //     And inner.material.ambient ← 1
    //     And r ← ray(point(0, 0, 0.75), vector(0, 0, -1))
    //   When c ← color_at(w, r)
    //   Then c = inner.material.color
    #[test]
    fn color_intersection_behind_ray() {
        let mut world = World::default();
        world.objects[0].material_mut().ambient = 1.;
        world.objects[1].material_mut().ambient = 1.;
        let ray = Ray::new(Tuple::point(0., 0., 0.75), Tuple::vector(0., 0., -1.));
        let color = world.color_at(&ray, 5);
        assert_eq!(&color, world.objects[1].material().pattern.colors()[0])
    }

    // Scenario: There is no shadow when nothing is collinear with point and light
    //   Given w ← default_world()
    //     And p ← point(0, 10, 0)
    //    Then is_shadowed(w, p) is false
    #[test]
    fn no_shadow_when_nothing_collinear_with_point_and_light() {
        let world = World::default();
        let point = Tuple::point(0., 10., 0.);
        assert!(!world.is_shadowed(&world.lights[0], &point));
    }

    // Scenario: The shadow when an object is between the point and the light
    //   Given w ← default_world()
    //     And p ← point(10, -10, 10)
    //    Then is_shadowed(w, p) is true
    #[test]
    fn shadow_when_object_between_point_and_light() {
        let world = World::default();
        let point = Tuple::point(10., -10., 10.);
        assert!(world.is_shadowed(&world.lights[0], &point));
    }

    // Scenario: There is no shadow when an object is behind the light
    //   Given w ← default_world()
    //     And p ← point(-20, 20, -20)
    //    Then is_shadowed(w, p) is false
    #[test]
    fn no_shadow_when_object_behind_light() {
        let world = World::default();
        let point = Tuple::point(-20., 20., -20.);
        assert!(!world.is_shadowed(&world.lights[0], &point));
    }

    // Scenario: There is no shadow when an object is behind the point
    //   Given w ← default_world()
    //     And p ← point(-2, 2, -2)
    //    Then is_shadowed(w, p) is false
    #[test]
    fn no_shadow_when_object_behind_point() {
        let world = World::default();
        let point = Tuple::point(-2., 2., -2.);
        assert!(!world.is_shadowed(&world.lights[0], &point));
    }

    // Scenario: shade_hit() is given an intersection in shadow
    //   Given w ← world()
    //     And w.light ← point_light(point(0, 0, -10), color(1, 1, 1))
    //     And s1 ← sphere()
    //     And s1 is added to w
    //     And s2 ← sphere() with:
    //       | transform | translation(0, 0, 10) |
    //     And s2 is added to w
    //     And r ← ray(point(0, 0, 5), vector(0, 0, 1))
    //     And i ← intersection(4, s2)
    //   When comps ← prepare_computations(i, r)
    //     And c ← shade_hit(w, comps)
    //   Then c = color(0.1, 0.1, 0.1)
    #[test]
    fn shade_hit_given_intersection_in_shadow() {
        let mut world = World::default();
        world.lights = vec![Light::new(
            Tuple::point(0., 0., -10.),
            Tuple::color(1., 1., 1.),
        )];
        let shape_a = Sphere::new(world.counter.increment());
        let mut shape_b = Sphere::new(world.counter.increment());
        shape_b.transform_mut().translate(0., 0., 10.);
        world.objects = vec![Box::new(shape_a), Box::new(shape_b)];
        let ray = Ray::new(Tuple::point(0., 0., 5.), Tuple::vector(0., 0., 1.));
        let intersection = Intersection::new(4., world.objects[1].clone());
        let comps = intersection.prepare_computations(&ray, None);
        let color = world.shade_hit(&comps, 5);
        assert_relative_eq!(color, Tuple::color(0.1, 0.1, 0.1));
    }

    // Scenario: The reflected color for a nonreflective material
    //   Given w ← default_world()
    //     And r ← ray(point(0, 0, 0), vector(0, 0, 1))
    //     And shape ← the second object in w
    //     And shape.material.ambient ← 1
    //     And i ← intersection(1, shape)
    //   When comps ← prepare_computations(i, r)
    //     And color ← reflected_color(w, comps)
    //   Then color = color(0, 0, 0)
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

    // Scenario: The reflected color for a reflective material
    //   Given w ← default_world()
    //     And shape ← plane() with:
    //       | material.reflective | 0.5                   |
    //       | transform           | translation(0, -1, 0) |
    //     And shape is added to w
    //     And r ← ray(point(0, 0, -3), vector(0, -√2/2, √2/2))
    //     And i ← intersection(√2, shape)
    //   When comps ← prepare_computations(i, r)
    //     And color ← reflected_color(w, comps)
    //   Then color = color(0.19032, 0.2379, 0.14274)
    #[test]
    fn reflected_color_for_reflective_material() {
        let mut world = World::default();
        let mut plane = Plane::new(world.counter.increment());
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
        assert_relative_eq!(
            color,
            Tuple::color(0.19032, 0.2379, 0.14274),
            epsilon = 1e-4f64
        );
    }

    // Scenario: shade_hit() with a reflective material
    //   Given w ← default_world()
    //     And shape ← plane() with:
    //       | material.reflective | 0.5                   |
    //       | transform           | translation(0, -1, 0) |
    //     And shape is added to w
    //     And r ← ray(point(0, 0, -3), vector(0, -√2/2, √2/2))
    //     And i ← intersection(√2, shape)
    //   When comps ← prepare_computations(i, r)
    //     And color ← shade_hit(w, comps)
    //   Then color = color(0.87677, 0.92436, 0.82918)
    #[test]
    fn shade_hit_with_reflective_material() {
        let mut world = World::default();
        let mut plane = Plane::new(world.counter.increment());
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
        assert_relative_eq!(
            color,
            Tuple::color(0.87677, 0.92436, 0.82918),
            epsilon = 1e-4f64
        );
    }

    // Scenario: color_at() with mutually reflective surfaces
    //   Given w ← world()
    //     And w.light ← point_light(point(0, 0, 0), color(1, 1, 1))
    //     And lower ← plane() with:
    //       | material.reflective | 1                     |
    //       | transform           | translation(0, -1, 0) |
    //     And lower is added to w
    //     And upper ← plane() with:
    //       | material.reflective | 1                    |
    //       | transform           | translation(0, 1, 0) |
    //     And upper is added to w
    //     And r ← ray(point(0, 0, 0), vector(0, 1, 0))
    //   Then color_at(w, r) should terminate successfully
    #[test]
    fn color_at_with_mutually_reflective_surfaces() {
        let mut world = World::default();
        world.lights[0].position = Tuple::point(0., 0., 0.);
        let mut lower = Plane::new(world.counter.increment());
        lower.material_mut().reflective = 1.;
        *lower.transform_mut() = Matrix::translation(0., -1., 0.);
        world.objects.push(Box::new(lower));
        let mut upper = Plane::new(world.counter.increment());
        upper.material_mut().reflective = 1.;
        *upper.transform_mut() = Matrix::translation(0., 1., 0.);
        world.objects.push(Box::new(upper));
        let ray = Ray::new(Tuple::point(0., 0., 0.), Tuple::vector(0., 1., 0.));
        let color = world.color_at(&ray, 5);
        assert_relative_eq!(color, Tuple::color(1.9, 1.9, 1.9));
    }

    // Scenario: The reflected color at the maximum recursive depth
    //   Given w ← default_world()
    //     And shape ← plane() with:
    //       | material.reflective | 0.5                   |
    //       | transform           | translation(0, -1, 0) |
    //     And shape is added to w
    //     And r ← ray(point(0, 0, -3), vector(0, -√2/2, √2/2))
    //     And i ← intersection(√2, shape)
    //   When comps ← prepare_computations(i, r)
    //     And color ← reflected_color(w, comps, 0)
    //   Then color = color(0, 0, 0)
    #[test]
    fn reflected_color_at_maximum_recursive_depth() {
        let mut world = World::default();
        let mut plane = Plane::new(world.counter.increment());
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

    // Scenario: The refracted color with an opaque surface
    //   Given w ← default_world()
    //     And shape ← the first object in w
    //     And r ← ray(point(0, 0, -5), vector(0, 0, 1))
    //     And xs ← intersections(4:shape, 6:shape)
    //   When comps ← prepare_computations(xs[0], r, xs)
    //     And c ← refracted_color(w, comps, 5)
    //   Then c = color(0, 0, 0)
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

    // Scenario: The refracted color at the maximum recursive depth
    //   Given w ← default_world()
    //     And shape ← the first object in w
    //     And shape has:
    //       | material.transparency     | 1.0 |
    //       | material.refractive_index | 1.5 |
    //     And r ← ray(point(0, 0, -5), vector(0, 0, 1))
    //     And xs ← intersections(4:shape, 6:shape)
    //   When comps ← prepare_computations(xs[0], r, xs)
    //     And c ← refracted_color(w, comps, 0)
    //   Then c = color(0, 0, 0)
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

    // Scenario: The refracted color under total internal reflection
    //   Given w ← default_world()
    //     And shape ← the first object in w
    //     And shape has:
    //       | material.transparency     | 1.0 |
    //       | material.refractive_index | 1.5 |
    //     And r ← ray(point(0, 0, √2/2), vector(0, 1, 0))
    //     And xs ← intersections(-√2/2:shape, √2/2:shape)
    //   # NOTE: this time you're inside the sphere, so you need
    //   # to look at the second intersection, xs[1], not xs[0]
    //   When comps ← prepare_computations(xs[1], r, xs)
    //     And c ← refracted_color(w, comps, 5)
    //   Then c = color(0, 0, 0)
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

    // Scenario: The refracted color with a refracted ray
    //   Given w ← default_world()
    //     And A ← the first object in w
    //     And A has:
    //       | material.ambient | 1.0            |
    //       | material.pattern | test_pattern() |
    //     And B ← the second object in w
    //     And B has:
    //       | material.transparency     | 1.0 |
    //       | material.refractive_index | 1.5 |
    //     And r ← ray(point(0, 0, 0.1), vector(0, 1, 0))
    //     And xs ← intersections(-0.9899:A, -0.4899:B, 0.4899:B, 0.9899:A)
    //   When comps ← prepare_computations(xs[2], r, xs)
    //     And c ← refracted_color(w, comps, 5)
    //   Then c = color(0, 0.99888, 0.04725)
    #[test]
    fn refracted_color_with_refracted_ray() {
        let mut world = World::default();
        world.objects[0].material_mut().ambient = 1.;
        world.objects[0].material_mut().pattern = Box::<TestPattern>::default();
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

    // Scenario: shade_hit() with a transparent material
    //   Given w ← default_world()
    //     And floor ← plane() with:
    //       | transform                 | translation(0, -1, 0) |
    //       | material.transparency     | 0.5                   |
    //       | material.refractive_index | 1.5                   |
    //     And floor is added to w
    //     And ball ← sphere() with:
    //       | material.color     | (1, 0, 0)                  |
    //       | material.ambient   | 0.5                        |
    //       | transform          | translation(0, -3.5, -0.5) |
    //     And ball is added to w
    //     And r ← ray(point(0, 0, -3), vector(0, -√2/2, √2/2))
    //     And xs ← intersections(√2:floor)
    //   When comps ← prepare_computations(xs[0], r, xs)
    //     And color ← shade_hit(w, comps, 5)
    //   Then color = color(0.93642, 0.68642, 0.68642)
    #[test]
    fn shade_hit_with_transparent_material() {
        let mut world = World::default();
        let mut floor = Box::new(Plane::new(world.counter.increment()));
        floor.material_mut().transparency = 0.5;
        floor.material_mut().refractive_index = 1.5;
        *floor.transform_mut() = Matrix::translation(0., -1., 0.);
        world.objects.push(floor.clone());
        let mut ball = Box::new(Sphere::new(world.counter.increment()));
        ball.material_mut().pattern = Box::new(SolidPattern::new(Tuple::color(1., 0., 0.)));
        ball.material_mut().ambient = 0.5;
        *ball.transform_mut() = Matrix::translation(0., -3.5, -0.5);
        world.objects.push(ball);
        let ray = Ray::new(
            Tuple::point(0., 0., -3.),
            Tuple::vector(0., -SQRT_2 / 2., SQRT_2 / 2.),
        );
        let intersections = vec![Intersection::new(SQRT_2, floor)];
        let comps = intersections[0].prepare_computations(&ray, Some(&intersections));
        let color = world.shade_hit(&comps, 5);
        assert_relative_eq!(
            color,
            Tuple::color(0.93642, 0.68642, 0.68642)
        );
    }

    // Scenario: shade_hit() with a reflective, transparent material
    //   Given w ← default_world()
    //     And r ← ray(point(0, 0, -3), vector(0, -√2/2, √2/2))
    //     And floor ← plane() with:
    //       | transform                 | translation(0, -1, 0) |
    //       | material.reflective       | 0.5                   |
    //       | material.transparency     | 0.5                   |
    //       | material.refractive_index | 1.5                   |
    //     And floor is added to w
    //     And ball ← sphere() with:
    //       | material.color     | (1, 0, 0)                  |
    //       | material.ambient   | 0.5                        |
    //       | transform          | translation(0, -3.5, -0.5) |
    //     And ball is added to w
    //     And xs ← intersections(√2:floor)
    //   When comps ← prepare_computations(xs[0], r, xs)
    //     And color ← shade_hit(w, comps, 5)
    //   Then color = color(0.93391, 0.69643, 0.69243)
    #[test]
    fn shade_hit_with_reflective_transparent_material() {
        let mut world = World::default();
        let mut floor = Box::new(Plane::new(world.counter.increment()));
        floor.material_mut().reflective = 0.5;
        floor.material_mut().transparency = 0.5;
        floor.material_mut().refractive_index = 1.5;
        *floor.transform_mut() = Matrix::translation(0., -1., 0.);
        world.objects.push(floor.clone());
        let mut ball = Box::new(Sphere::new(world.counter.increment()));
        ball.material_mut().pattern = Box::new(SolidPattern::new(Tuple::color(1., 0., 0.)));
        ball.material_mut().ambient = 0.5;
        *ball.transform_mut() = Matrix::translation(0., -3.5, -0.5);
        world.objects.push(ball);
        let ray = Ray::new(
            Tuple::point(0., 0., -3.),
            Tuple::vector(0., -SQRT_2 / 2., SQRT_2 / 2.),
        );
        let intersections = vec![Intersection::new(SQRT_2, floor)];
        let comps = intersections[0].prepare_computations(&ray, Some(&intersections));
        let color = world.shade_hit(&comps, 5);
        assert_relative_eq!(
            color,
            Tuple::color(0.93391, 0.69643, 0.69243)
        );
    }
}
