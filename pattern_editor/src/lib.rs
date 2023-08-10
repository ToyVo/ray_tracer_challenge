use wasm_bindgen::prelude::*;

// When the `wee_alloc` feature is enabled, this uses `wee_alloc` as the global
// allocator.
//
// If you don't want to use `wee_alloc`, you can safely delete this.
#[cfg(feature = "wee_alloc")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

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
    // pub fn set_pixel(&mut self, x: u32, y: u32, color: Vector3<f32>) {
    //     let index = (y * self.width + x) as usize;
    //     self.data[index] = 0xFF000000;
    // }
}

// #[wasm_bindgen]
// impl PatternImage {
//     pub fn new(width: u32, height: u32) -> PatternImage {
//         let size = (width * height) as usize;
//         let data = vec![0xFF000000; size];
//
//         let mut s1 = StripePattern::new(
//             Box::new(SolidPattern::new(Vector3::new(1., 1., 1.))),
//             Box::new(SolidPattern::new(Vector3::new(0.75, 0.75, 0.75))),
//         );
//         // *s1.transform_mut() = Matrix::scaling(0.1, 0.1, 0.1);
//         let mut s2 = StripePattern::new(
//             Box::new(SolidPattern::new(Vector3::new(0.5, 0.5, 0.5))),
//             Box::new(SolidPattern::new(Vector3::new(0.25, 0.25, 0.25))),
//         );
//         // *s2.transform_mut() = Matrix::scaling(0.1, 0.1, 0.1) * Matrix::rotation_y(PI / 2.);
//         let pattern =
//             PerturbedPattern::new(Box::new(BlendedPattern::new(Box::new(s1), Box::new(s2))));
//
//         let mut image = PatternImage {
//             width,
//             height,
//             data,
//         };
//
//         for row in 0..height {
//             for col in 0..width {
//                 let point = Vector4::new(col as f64, 0., row as f64, 1.);
//                 let color = pattern.pattern_at(&point);
//                 image.set_pixel(col, row, color);
//             }
//         }
//
//         image
//     }
//     pub fn width(&self) -> u32 {
//         self.width
//     }
//     pub fn height(&self) -> u32 {
//         self.height
//     }
//     pub fn data(&self) -> *const u32 {
//         self.data.as_ptr()
//     }
// }
