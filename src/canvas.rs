use crate::Tuple;

pub struct Canvas {
    data: Vec<Tuple>,
    rows: usize,
    cols: usize,
}

impl Canvas {
    pub fn new(rows: usize, cols: usize, init: Tuple) -> Canvas {
        let data = vec![init; rows * cols];
        Canvas { data, rows, cols }
    }

    pub fn from_vec(rows: usize, cols: usize, data: Vec<Tuple>) -> Canvas {
        Canvas { data, rows, cols }
    }

    pub fn get(&self, row: usize, col: usize) -> &Tuple {
        assert!(row < self.rows && col < self.cols);
        &self.data[row * self.cols + col]
    }

    pub fn set(&mut self, row: usize, col: usize, value: Tuple) {
        assert!(row < self.rows && col < self.cols);
        self.data[row * self.cols + col] = value;
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn width(&self) -> usize {
        self.cols
    }

    pub fn height(&self) -> usize {
        self.rows
    }
}

impl std::fmt::Display for Canvas {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut ppm = String::with_capacity(self.cols() * self.rows() * 12);
        ppm.push_str("P3\n");
        ppm.push_str(&format!("{} {}\n", self.cols(), self.rows()));
        ppm.push_str("255\n");

        for row in 0..self.rows() {
            let mut row_colors: Vec<u8> = Vec::with_capacity(self.cols() * 3);
            for col in 0..self.cols() {
                let color = self.get(row, col);
                row_colors.push((color.r() * 255.0).round() as u8);
                row_colors.push((color.g() * 255.0).round() as u8);
                row_colors.push((color.b() * 255.0).round() as u8);
            }
            let mut row = String::new();
            let mut line_length = 0;
            for (_i, color) in row_colors.iter().enumerate() {
                let color_str = format!("{} ", color);
                let color_len = color_str.len();
                if line_length + color_len > 70 {
                    row = row.trim_end().to_string();
                    row.push('\n');
                    line_length = 0;
                }
                row.push_str(&color_str);
                line_length += color_len;
            }
            ppm.push_str(row.trim_end());
            ppm.push('\n');
        }
        write!(f, "{}", ppm)
    }
}

#[cfg(test)]
mod testing {
    use super::*;

    #[test]
    fn creating_canvas() {
        let c = Canvas::new(20, 10, Tuple::color(0.0, 0.0, 0.0));
        assert_eq!(c.width(), 10);
        assert_eq!(c.height(), 20);
    }

    #[test]
    fn writing_pixels_to_canvas() {
        let mut c = Canvas::new(20, 10, Tuple::color(0.0, 0.0, 0.0));
        let red = Tuple::color(1.0, 0.0, 0.0);
        c.set(2, 3, red.clone());
        assert_eq!(c.get(2, 3), &red);
    }

    #[test]
    fn constructing_ppm_header() {
        let c = Canvas::new(3, 5, Tuple::color(0.0, 0.0, 0.0));
        let ppm = format!("{}", c);
        let lines: Vec<&str> = ppm.split('\n').collect();
        assert_eq!(lines[0], "P3");
        assert_eq!(lines[1], "5 3");
        assert_eq!(lines[2], "255");
    }

    #[test]
    fn constructing_ppm_pixel_data() {
        let mut c = Canvas::new(3, 5, Tuple::color(0.0, 0.0, 0.0));
        let c1 = Tuple::color(1.5, 0.0, 0.0);
        let c2 = Tuple::color(0.0, 0.5, 0.0);
        let c3 = Tuple::color(-0.5, 0.0, 1.0);
        c.set(0, 0, c1);
        c.set(1, 2, c2);
        c.set(2, 4, c3);
        let ppm = format!("{}", c);
        let lines: Vec<&str> = ppm.split('\n').collect();
        assert_eq!(lines[3], "255 0 0 0 0 0 0 0 0 0 0 0 0 0 0");
        assert_eq!(lines[4], "0 0 0 0 0 0 0 128 0 0 0 0 0 0 0");
        assert_eq!(lines[5], "0 0 0 0 0 0 0 0 0 0 0 0 0 0 255");
    }

    #[test]
    fn splitting_long_lines_in_ppm_files() {
        let c = Canvas::new(2, 10, Tuple::color(1.0, 0.8, 0.6));
        let ppm = format!("{}", c);
        let lines: Vec<&str> = ppm.split('\n').collect();
        assert_eq!(
            lines[3],
            "255 204 153 255 204 153 255 204 153 255 204 153 255 204 153 255 204"
        );
        assert_eq!(
            lines[4],
            "153 255 204 153 255 204 153 255 204 153 255 204 153"
        );
        assert_eq!(
            lines[5],
            "255 204 153 255 204 153 255 204 153 255 204 153 255 204 153 255 204"
        );
        assert_eq!(
            lines[6],
            "153 255 204 153 255 204 153 255 204 153 255 204 153"
        );
    }

    #[test]
    fn ppm_files_are_terminated_by_a_newline_character() {
        let c = Canvas::new(3, 5, Tuple::color(0.0, 0.0, 0.0));
        let ppm = format!("{}", c);
        let lines: Vec<&str> = ppm.split('\n').collect();
        assert_eq!(lines[lines.len() - 1], "");
    }
}
