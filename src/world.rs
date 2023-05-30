use crate::{Light, Matrix, Sphere, Tuple, Ray, Intersection, Computations};

pub struct World {
    pub objects: Vec<Sphere>,
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
        s1.material.color = Tuple::color(0.8, 1.0, 0.6);
        s1.material.diffuse = 0.7;
        s1.material.specular = 0.2;
        let mut s2 = Sphere::new();
        s2.transform = Matrix::scaling(0.5, 0.5, 0.5);
        let light = Light::new(Tuple::point(-10.0, 10.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        World {
            objects: vec![s1, s2],
            lights: vec![light],
        }
    }

    pub fn intersects(&self, ray: &Ray) -> Vec<Intersection> {
        let mut intersections: Vec<Intersection> = vec![];
        for object in &self.objects {
            let mut object_intersections = object.intersects(ray);
            intersections.append(&mut object_intersections);
        }
        intersections.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
        intersections
    }

    pub fn shade_hit(&self, comps: &Computations) -> Tuple {
        self.lights.iter().fold(Tuple::color(0.,0.,0.), |acc, light| {
            acc + comps.object.material.lighting(light, &comps.point, &comps.eyev,&comps.normalv)
        })
    }

    pub fn color_at(&self, ray: &Ray) -> Tuple {
        let intersections = self.intersects(ray);
        if intersections.is_empty() {
            Tuple::color(0.,0.,0.)
        } else {
            let hit = Intersection::hit(&intersections).unwrap();
            let comps = hit.prepare_computations(ray);
            self.shade_hit(&comps)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn creating_world() {
        let world = World::new();
        assert_eq!(world.objects.len(), 0);
        assert_eq!(world.lights.len(), 0);
    }

    #[test]
    fn default_world() {
        let light = Light::new(Tuple::point(-10.0, 10.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        let mut s1 = Sphere::new();
        s1.material.color = Tuple::color(0.8, 1.0, 0.6);
        s1.material.diffuse = 0.7;
        s1.material.specular = 0.2;
        let mut s2 = Sphere::new();
        s2.transform = Matrix::scaling(0.5, 0.5, 0.5);
        let world = World::default();
        assert_eq!(world.lights.len(), 1);
        assert_eq!(world.lights[0], light);
        assert_eq!(world.objects.len(), 2);
        assert_eq!(world.objects[0], s1);
        assert_eq!(world.objects[1], s2);
    }

    #[test]
    fn intersect_world_with_ray() {
        let world = World::default();
        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let intersections = world.intersects(&ray);
        assert_eq!(intersections.len(), 4);
        assert_eq!(intersections[0].t, 4.0);
        assert_eq!(intersections[1].t, 4.5);
        assert_eq!(intersections[2].t, 5.5);
        assert_eq!(intersections[3].t, 6.0);
    }

    #[test]
    fn shading_intersection() {
        let world = World::default();
        let ray = Ray::new(Tuple::point(0.,0.,-5.),Tuple::vector(0.,0.,1.));
        let shape = &world.objects[0];
        let intersection = Intersection::new(4., shape);
        let comps = intersection.prepare_computations(&ray);
        let color = world.shade_hit(&comps);
        assert!(color.nearly_equals(&Tuple::color(0.38066,0.47583,0.2855), 0.00001));
    }

    #[test]
    fn shading_intersection_inside() {
        let mut world = World::default();
        world.lights = vec![Light::new(Tuple::point(0.,0.25,0.), Tuple::color(1.,1.,1.))];
        let ray = Ray::new(Tuple::point(0.,0.,0.),Tuple::vector(0.,0.,1.));
        let shape = &world.objects[1];
        let intersection = Intersection::new(0.5, shape);
        let comps = intersection.prepare_computations(&ray);
        let color = world.shade_hit(&comps);
        assert!(color.nearly_equals(&Tuple::color(0.90498,0.90498,0.90498), 0.00001));
    }

    #[test]
    fn color_when_ray_misses() {
        let world = World::default();
        let ray = Ray::new(Tuple::point(0.,0.,-5.), Tuple::vector(0.,1.,0.));
        let color = world.color_at(&ray);
        assert_eq!(color, Tuple::color(0.,0.,0.))
    }

    #[test]
    fn color_when_ray_hits() {
        let world = World::default();
        let ray = Ray::new(Tuple::point(0.,0.,-5.), Tuple::vector(0.,0.,1.));
        let color = world.color_at(&ray);
        assert!(color.nearly_equals(&Tuple::color(0.38066,0.47583,0.2855), 0.00001));
    }

    #[test]
    fn color_intersection_behind_ray() {
        let mut world = World::default();
        world.objects[0].material.ambient = 1.;
        world.objects[1].material.ambient = 1.;
        let ray = Ray::new(Tuple::point(0.,0.,0.75), Tuple::vector(0.,0.,-1.));
        let color = world.color_at(&ray);
        assert_eq!(color, world.objects[1].material.color)
    }
}