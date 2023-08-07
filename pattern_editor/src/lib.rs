use bracket_noise::prelude::*;
use dyn_clone::DynClone;
use nalgebra::{Transform3, Vector3, Vector4};
use wasm_bindgen::prelude::*;
use web_sys::console;
use std::f64::consts::PI;

// When the `wee_alloc` feature is enabled, this uses `wee_alloc` as the global
// allocator.
//
// If you don't want to use `wee_alloc`, you can safely delete this.
#[cfg(feature = "wee_alloc")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

// This is like the `main` function, except for JavaScript.
#[wasm_bindgen(start)]
pub fn main_js() -> Result<(), JsValue> {
    // This provides better error messages in debug mode.
    // It's disabled in release mode so it doesn't bloat up the file size.
    #[cfg(debug_assertions)]
    console_error_panic_hook::set_once();

    // Your code goes here!
    console::log_1(&JsValue::from_str("Hello world!"));

    Ok(())
}

const CELL_SIZE: u32 = 5;
const DEAD_COLOR: u32 = 0xFFFFFFFF;
const ALIVE_COLOR: u32 = 0xFF000000;
const GRID_COLOR: u32 = 0xFFCCCCCC;

fn set_cell(data: &mut Vec<u32>, row: u32, col: u32, canvas_height: u32, color: u32) {
    let row = row * CELL_SIZE + row + 1;
    let col = col * CELL_SIZE + col + 1;
    for x in 0..CELL_SIZE {
        for y in 0..CELL_SIZE {
            let index = ((row + y) * canvas_height + col + x) as usize;
            data[index] = color;
        }
    }
}

#[wasm_bindgen]
pub struct Universe {
    width: u32,
    height: u32,
    data: Vec<u32>,
}

impl Universe {
    fn get_index(&self, row: u32, col: u32) -> usize {
        let row = row * CELL_SIZE + row + 1;
        let col = col * CELL_SIZE + col + 1;
        (row * self.canvas_height() + col) as usize
    }

    fn live_neighbor_count(&self, row: u32, column: u32) -> u8 {
        let mut count = 0;
        for delta_row in [self.height - 1, 0, 1].iter().cloned() {
            for delta_col in [self.width - 1, 0, 1].iter().cloned() {
                if delta_row == 0 && delta_col == 0 {
                    continue;
                }

                let neighbor_row = (row + delta_row) % self.height;
                let neighbor_col = (column + delta_col) % self.width;
                let idx = self.get_index(neighbor_row, neighbor_col);
                count += (self.data[idx] == ALIVE_COLOR) as u8;
            }
        }
        count
    }

    fn set_cell(&mut self, row: u32, col: u32, color: u32) {
        let row = row * CELL_SIZE + row + 1;
        let col = col * CELL_SIZE + col + 1;
        for x in 0..CELL_SIZE {
            for y in 0..CELL_SIZE {
                let index = ((row + y) * self.canvas_height() + col + x) as usize;
                self.data[index] = color;
            }
        }
    }
}

/// Public methods, exported to JavaScript.
#[wasm_bindgen]
impl Universe {
    pub fn tick(&mut self) {
        let mut next = self.data.clone();

        for row in 0..self.height {
            for col in 0..self.width {
                let idx = self.get_index(row, col);
                let color = self.data[idx];
                let live_neighbors = self.live_neighbor_count(row, col);

                let next_color = match (color, live_neighbors) {
                    // Rule 1: Any live color with fewer than two live neighbours
                    // dies, as if caused by underpopulation.
                    (ALIVE_COLOR, x) if x < 2 => DEAD_COLOR,
                    // Rule 2: Any live color with two or three live neighbours
                    // lives on to the next generation.
                    (ALIVE_COLOR, 2) | (ALIVE_COLOR, 3) => ALIVE_COLOR,
                    // Rule 3: Any live color with more than three live
                    // neighbours dies, as if by overpopulation.
                    (ALIVE_COLOR, x) if x > 3 => DEAD_COLOR,
                    // Rule 4: Any dead color with exactly three live neighbours
                    // becomes a live color, as if by reproduction.
                    (DEAD_COLOR, 3) => ALIVE_COLOR,
                    // All other colors remain in the same state.
                    (otherwise, _) => otherwise,
                };

                set_cell(&mut next, row, col, self.canvas_height(), next_color);
            }
        }

        self.data = next;
    }

    pub fn new() -> Universe {
        let width = 64;
        let height = 64;
        let canvas_height = height * CELL_SIZE + height + 1;
        let canvas_width = width * CELL_SIZE + width + 1;
        let size = (canvas_width * canvas_height) as usize;
        let mut data = vec![DEAD_COLOR; size];

        for row in 0..canvas_height {
            for col in 0..canvas_width {
                if row % (CELL_SIZE + 1) == 0 || col % (CELL_SIZE + 1) == 0 {
                    let index = (row * canvas_height + col) as usize;
                    data[index] = GRID_COLOR;
                }
            }
        }

        let mut universe = Universe {
            width,
            height,
            data,
        };

        universe.set_cell(0, 1, ALIVE_COLOR);
        universe.set_cell(1, 2, ALIVE_COLOR);
        universe.set_cell(2, 0, ALIVE_COLOR);
        universe.set_cell(2, 1, ALIVE_COLOR);
        universe.set_cell(2, 2, ALIVE_COLOR);

        universe
    }

    pub fn render(&self) -> String {
        let mut result = String::new();
        for row in 0..self.canvas_height() {
            for col in 0..self.canvas_width() {
                let idx = (row * self.canvas_width() + col) as usize;
                let color = self.data[idx];
                let symbol = match color {
                    ALIVE_COLOR => "◼",
                    DEAD_COLOR => "◻",
                    GRID_COLOR => "#",
                    _ => " ",
                };
                result.push_str(symbol);
            }
            result.push_str("\n");
        }
        result
    }

    pub fn canvas_width(&self) -> u32 {
        self.width * CELL_SIZE + (self.width + 1)
    }

    pub fn canvas_height(&self) -> u32 {
        self.height * CELL_SIZE + (self.height + 1)
    }

    pub fn data(&self) -> *const u32 {
        self.data.as_ptr()
    }

    pub fn size(&self) -> usize {
        self.data.len()
    }
}

#[wasm_bindgen]
pub struct PatternImage {
    width: u32,
    height: u32,
    data: Vec<u32>,
}

impl PatternImage {
    pub fn set_pixel(&mut self, x: u32, y: u32, color: Vector3<f32>) {
        let index = (y * self.width + x) as usize;
        self.data[index] = 0xFF000000;
    }
}

#[wasm_bindgen]
impl PatternImage {
    pub fn new(width: u32, height: u32) -> PatternImage {
        let size = (width * height) as usize;
        let data = vec![0xFF000000; size];

        let mut s1 = StripePattern::new(
            Box::new(SolidPattern::new(Vector3::new(1., 1., 1.))),
            Box::new(SolidPattern::new(Vector3::new(0.75, 0.75, 0.75))),
        );
        // *s1.transform_mut() = Matrix::scaling(0.1, 0.1, 0.1);
        let mut s2 = StripePattern::new(
            Box::new(SolidPattern::new(Vector3::new(0.5, 0.5, 0.5))),
            Box::new(SolidPattern::new(Vector3::new(0.25, 0.25, 0.25))),
        );
        // *s2.transform_mut() = Matrix::scaling(0.1, 0.1, 0.1) * Matrix::rotation_y(PI / 2.);
        let pattern =
            PerturbedPattern::new(Box::new(BlendedPattern::new(Box::new(s1), Box::new(s2))));

        let mut image = PatternImage {
            width,
            height,
            data,
        };

        for row in 0..height {
            for col in 0..width {
                let point = Vector4::new(col as f64, 0., row as f64, 1.);
                let color = pattern.pattern_at(&point);
                image.set_pixel(col, row, color);
            }
        }

        image
    }
    pub fn width(&self) -> u32 {
        self.width
    }
    pub fn height(&self) -> u32 {
        self.height
    }
    pub fn data(&self) -> *const u32 {
        self.data.as_ptr()
    }
}

pub trait Pattern: std::fmt::Debug + DynClone {
    fn pattern_at(&self, point: &Vector4<f64>) -> Vector3<f32>;
    fn colors(&self) -> Vec<&Vector3<f32>>;
    fn transform(&self) -> &Transform3<f64>;
    fn transform_mut(&mut self) -> &mut Transform3<f64>;
    fn pattern_at_pattern(&self, object: &dyn Pattern, world_point: &Vector4<f64>) -> Vector3<f32> {
        let object_point = object.transform().try_inverse().unwrap().matrix() * world_point;
        let pattern_point = self.transform().try_inverse().unwrap().matrix() * &object_point;
        self.pattern_at(&pattern_point)
    }
}

dyn_clone::clone_trait_object!(Pattern);

impl PartialEq for dyn Pattern {
    fn eq(&self, other: &Self) -> bool {
        self.colors() == other.colors()
    }
}

#[derive(Debug, Clone)]
pub struct StripePattern {
    pattern_a: Box<dyn Pattern>,
    pattern_b: Box<dyn Pattern>,
    transform: Transform3<f64>,
}

impl StripePattern {
    pub fn new(pattern_a: Box<dyn Pattern>, pattern_b: Box<dyn Pattern>) -> StripePattern {
        StripePattern {
            pattern_a,
            pattern_b,
            transform: Transform3::identity(),
        }
    }
}

impl Pattern for StripePattern {
    fn pattern_at(&self, point: &Vector4<f64>) -> Vector3<f32> {
        if point.x.floor() as i32 % 2 == 0 {
            self.pattern_a.pattern_at_pattern(self, point)
        } else {
            self.pattern_b.pattern_at_pattern(self, point)
        }
    }

    fn colors(&self) -> Vec<&Vector3<f32>> {
        let mut colors = self.pattern_a.colors();
        colors.append(&mut self.pattern_b.colors());
        colors
    }

    fn transform(&self) -> &Transform3<f64> {
        &self.transform
    }

    fn transform_mut(&mut self) -> &mut Transform3<f64> {
        &mut self.transform
    }
}

#[derive(Debug, Clone)]
pub struct SolidPattern {
    color: Vector3<f32>,
    transform: Transform3<f64>,
}

impl SolidPattern {
    pub fn new(color: Vector3<f32>) -> SolidPattern {
        SolidPattern {
            color,
            transform: Transform3::identity(),
        }
    }
}

impl Pattern for SolidPattern {
    fn pattern_at(&self, _point: &Vector4<f64>) -> Vector3<f32> {
        self.color.clone()
    }

    fn colors(&self) -> Vec<&Vector3<f32>> {
        vec![&self.color]
    }

    fn transform(&self) -> &Transform3<f64> {
        &self.transform
    }

    fn transform_mut(&mut self) -> &mut Transform3<f64> {
        &mut self.transform
    }
}

#[derive(Debug, Clone)]
pub struct GradientPattern {
    color_a: Vector3<f32>,
    color_b: Vector3<f32>,
    transform: Transform3<f64>,
}

impl GradientPattern {
    pub fn new(color_a: Vector3<f32>, color_b: Vector3<f32>) -> GradientPattern {
        GradientPattern {
            color_a,
            color_b,
            transform: Transform3::identity(),
        }
    }
}

impl Pattern for GradientPattern {
    fn pattern_at(&self, point: &Vector4<f64>) -> Vector3<f32> {
        let distance = &self.color_b - &self.color_a;
        let fraction = point.x as f32 - (point.x as f32).floor();
        &self.color_a + &distance * fraction
    }

    fn colors(&self) -> Vec<&Vector3<f32>> {
        vec![&self.color_a, &self.color_b]
    }

    fn transform(&self) -> &Transform3<f64> {
        &self.transform
    }

    fn transform_mut(&mut self) -> &mut Transform3<f64> {
        &mut self.transform
    }
}

#[derive(Debug, Clone)]
pub struct RingPattern {
    pattern_a: Box<dyn Pattern>,
    pattern_b: Box<dyn Pattern>,
    transform: Transform3<f64>,
}

impl RingPattern {
    pub fn new(pattern_a: Box<dyn Pattern>, pattern_b: Box<dyn Pattern>) -> RingPattern {
        RingPattern {
            pattern_a,
            pattern_b,
            transform: Transform3::identity(),
        }
    }
}

impl Pattern for RingPattern {
    fn pattern_at(&self, point: &Vector4<f64>) -> Vector3<f32> {
        if (point.x.powi(2) + point.z.powi(2)).sqrt().floor() as i32 % 2 == 0 {
            self.pattern_a.pattern_at_pattern(self, point)
        } else {
            self.pattern_b.pattern_at_pattern(self, point)
        }
    }

    fn colors(&self) -> Vec<&Vector3<f32>> {
        let mut colors = self.pattern_a.colors();
        colors.append(&mut self.pattern_b.colors());
        colors
    }

    fn transform(&self) -> &Transform3<f64> {
        &self.transform
    }

    fn transform_mut(&mut self) -> &mut Transform3<f64> {
        &mut self.transform
    }
}

#[derive(Debug, Clone)]
pub struct CheckeredPattern {
    pattern_a: Box<dyn Pattern>,
    pattern_b: Box<dyn Pattern>,
    transform: Transform3<f64>,
}

impl CheckeredPattern {
    pub fn new(pattern_a: Box<dyn Pattern>, pattern_b: Box<dyn Pattern>) -> CheckeredPattern {
        CheckeredPattern {
            pattern_a,
            pattern_b,
            transform: Transform3::identity(),
        }
    }
}

impl Pattern for CheckeredPattern {
    fn pattern_at(&self, point: &Vector4<f64>) -> Vector3<f32> {
        if (point.x.floor() + point.y.floor() + point.z.floor()) as i32 % 2 == 0 {
            self.pattern_a.pattern_at_pattern(self, point)
        } else {
            self.pattern_b.pattern_at_pattern(self, point)
        }
    }

    fn colors(&self) -> Vec<&Vector3<f32>> {
        let mut colors = self.pattern_a.colors();
        colors.append(&mut self.pattern_b.colors());
        colors
    }

    fn transform(&self) -> &Transform3<f64> {
        &self.transform
    }

    fn transform_mut(&mut self) -> &mut Transform3<f64> {
        &mut self.transform
    }
}

#[derive(Debug, Clone)]
pub struct BlendedPattern {
    pattern_a: Box<dyn Pattern>,
    pattern_b: Box<dyn Pattern>,
    transform: Transform3<f64>,
}

impl BlendedPattern {
    pub fn new(pattern_a: Box<dyn Pattern>, pattern_b: Box<dyn Pattern>) -> BlendedPattern {
        BlendedPattern {
            pattern_a,
            pattern_b,
            transform: Transform3::identity(),
        }
    }
}

impl Pattern for BlendedPattern {
    fn pattern_at(&self, point: &Vector4<f64>) -> Vector3<f32> {
        let color_a = self.pattern_a.pattern_at_pattern(self, point);
        let color_b = self.pattern_b.pattern_at_pattern(self, point);
        (&color_a + &color_b) / 2.
    }

    fn colors(&self) -> Vec<&Vector3<f32>> {
        let mut colors = self.pattern_a.colors();
        colors.append(&mut self.pattern_b.colors());
        colors
    }

    fn transform(&self) -> &Transform3<f64> {
        &self.transform
    }

    fn transform_mut(&mut self) -> &mut Transform3<f64> {
        &mut self.transform
    }
}

#[derive(Debug, Clone)]
pub struct PerturbedPattern {
    pattern: Box<dyn Pattern>,
    transform: Transform3<f64>,
}

impl PerturbedPattern {
    pub fn new(pattern: Box<dyn Pattern>) -> PerturbedPattern {
        PerturbedPattern {
            pattern,
            transform: Transform3::identity(),
        }
    }
}

impl Pattern for PerturbedPattern {
    fn pattern_at(&self, point: &Vector4<f64>) -> Vector3<f32> {
        let mut noise = FastNoise::new();
        noise.set_noise_type(NoiseType::SimplexFractal);
        self.pattern.pattern_at_pattern(
            self,
            &Vector4::new(
                point.x + noise.get_noise3d(point.x as f32, point.y as f32, point.z as f32) as f64,
                point.y
                    + noise.get_noise3d(point.x as f32, point.y as f32, point.z as f32 + 1.) as f64,
                point.z
                    + noise.get_noise3d(point.x as f32, point.y as f32, point.z as f32 + 2.) as f64,
                1.,
            ),
        )
    }

    fn colors(&self) -> Vec<&Vector3<f32>> {
        self.pattern.colors()
    }

    fn transform(&self) -> &Transform3<f64> {
        &self.transform
    }

    fn transform_mut(&mut self) -> &mut Transform3<f64> {
        &mut self.transform
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn universe_after_1_tick() {
        let mut universe = Universe::new();
        let c = universe.live_neighbor_count(1, 0);
        universe.tick();
    }

    #[test]
    fn creating_stripe_pattern() {
        let white = Vector3::new(1.0, 1.0, 1.0);
        let black = Vector3::new(0.0, 0.0, 0.0);
        let pattern = StripePattern::new(
            Box::new(SolidPattern::new(white.clone())),
            Box::new(SolidPattern::new(black.clone())),
        );
        assert_eq!(pattern.pattern_a.colors(), vec![&white]);
        assert_eq!(pattern.pattern_b.colors(), vec![&black]);
    }

    #[test]
    fn stripe_pattern_constant_in_y() {
        let white = Vector3::new(1.0, 1.0, 1.0);
        let black = Vector3::new(0.0, 0.0, 0.0);
        let pattern = StripePattern::new(
            Box::new(SolidPattern::new(white.clone())),
            Box::new(SolidPattern::new(black.clone())),
        );
        assert_eq!(pattern.pattern_at(&Vector4::new(0.0, 0.0, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&Vector4::new(0.0, 1.0, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&Vector4::new(0.0, 2.0, 0.0, 1.)), white);
    }

    #[test]
    fn stripe_pattern_constant_in_z() {
        let white = Vector3::new(1.0, 1.0, 1.0);
        let black = Vector3::new(0.0, 0.0, 0.0);
        let pattern = StripePattern::new(
            Box::new(SolidPattern::new(white.clone())),
            Box::new(SolidPattern::new(black.clone())),
        );
        assert_eq!(pattern.pattern_at(&Vector4::new(0.0, 0.0, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&Vector4::new(0.0, 0.0, 1.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&Vector4::new(0.0, 0.0, 2.0, 1.)), white);
    }

    #[test]
    fn stripe_pattern_alternates_in_x() {
        let white = Vector3::new(1.0, 1.0, 1.0);
        let black = Vector3::new(0.0, 0.0, 0.0);
        let pattern = StripePattern::new(
            Box::new(SolidPattern::new(white.clone())),
            Box::new(SolidPattern::new(black.clone())),
        );
        assert_eq!(pattern.pattern_at(&Vector4::new(0.0, 0.0, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&Vector4::new(0.9, 0.0, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&Vector4::new(1.0, 0.0, 0.0, 1.)), black);
        assert_eq!(pattern.pattern_at(&Vector4::new(-0.1, 0.0, 0.0, 1.)), black);
        assert_eq!(pattern.pattern_at(&Vector4::new(-1.0, 0.0, 0.0, 1.)), black);
        assert_eq!(pattern.pattern_at(&Vector4::new(-1.1, 0.0, 0.0, 1.)), white);
    }

    #[test]
    fn gradient_linearly_interpolates_between_colors() {
        let white = Vector3::new(1.0, 1.0, 1.0);
        let black = Vector3::new(0.0, 0.0, 0.0);
        let pattern = GradientPattern::new(white.clone(), black.clone());
        assert_eq!(pattern.pattern_at(&Vector4::new(0.0, 0.0, 0.0, 1.)), white);
        assert_eq!(
            pattern.pattern_at(&Vector4::new(0.25, 0.0, 0.0, 1.)),
            Vector3::new(0.75, 0.75, 0.75)
        );
        assert_eq!(
            pattern.pattern_at(&Vector4::new(0.5, 0.0, 0.0, 1.)),
            Vector3::new(0.5, 0.5, 0.5)
        );
        assert_eq!(
            pattern.pattern_at(&Vector4::new(0.75, 0.0, 0.0, 1.)),
            Vector3::new(0.25, 0.25, 0.25)
        );
    }

    #[test]
    fn ring_should_extend_in_both_x_and_z() {
        let white = Vector3::new(1.0, 1.0, 1.0);
        let black = Vector3::new(0.0, 0.0, 0.0);
        let pattern = RingPattern::new(
            Box::new(SolidPattern::new(white.clone())),
            Box::new(SolidPattern::new(black.clone())),
        );
        assert_eq!(pattern.pattern_at(&Vector4::new(0.0, 0.0, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&Vector4::new(1.0, 0.0, 0.0, 1.)), black);
        assert_eq!(pattern.pattern_at(&Vector4::new(0.0, 0.0, 1.0, 1.)), black);
        assert_eq!(
            pattern.pattern_at(&Vector4::new(0.708, 0.0, 0.708, 1.)),
            black
        );
    }

    #[test]
    fn checkers_should_repeat_in_x() {
        let white = Vector3::new(1.0, 1.0, 1.0);
        let black = Vector3::new(0.0, 0.0, 0.0);
        let pattern = CheckeredPattern::new(
            Box::new(SolidPattern::new(white.clone())),
            Box::new(SolidPattern::new(black.clone())),
        );
        assert_eq!(pattern.pattern_at(&Vector4::new(0.0, 0.0, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&Vector4::new(0.99, 0.0, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&Vector4::new(1.01, 0.0, 0.0, 1.)), black);
    }

    #[test]
    fn checkers_should_repeat_in_y() {
        let white = Vector3::new(1.0, 1.0, 1.0);
        let black = Vector3::new(0.0, 0.0, 0.0);
        let pattern = CheckeredPattern::new(
            Box::new(SolidPattern::new(white.clone())),
            Box::new(SolidPattern::new(black.clone())),
        );
        assert_eq!(pattern.pattern_at(&Vector4::new(0.0, 0.0, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&Vector4::new(0.0, 0.99, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&Vector4::new(0.0, 1.01, 0.0, 1.)), black);
    }

    #[test]
    fn checkers_should_repeat_in_z() {
        let white = Vector3::new(1.0, 1.0, 1.0);
        let black = Vector3::new(0.0, 0.0, 0.0);
        let pattern = CheckeredPattern::new(
            Box::new(SolidPattern::new(white.clone())),
            Box::new(SolidPattern::new(black.clone())),
        );
        assert_eq!(pattern.pattern_at(&Vector4::new(0.0, 0.0, 0.0, 1.)), white);
        assert_eq!(pattern.pattern_at(&Vector4::new(0.0, 0.0, 0.99, 1.)), white);
        assert_eq!(pattern.pattern_at(&Vector4::new(0.0, 0.0, 1.01, 1.)), black);
    }
}
