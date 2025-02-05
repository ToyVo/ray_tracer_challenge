use crate::Tuple;

pub struct Canvas {
    data: Vec<Tuple>,
    cols: usize,
    rows: usize,
}

impl Canvas {
    pub fn new(cols: usize, rows: usize, init: Tuple) -> Canvas {
        let data = vec![init; cols * rows];
        Canvas { data, cols, rows }
    }

    pub fn from_vec(cols: usize, rows: usize, data: Vec<Tuple>) -> Canvas {
        Canvas { data, cols, rows }
    }

    pub fn pixel_at(&self, x: usize, y: usize) -> &Tuple {
        assert!(x < self.cols && y < self.rows);
        &self.data[y * self.cols + x]
    }

    pub fn write_pixel(&mut self, x: usize, y: usize, value: Tuple) {
        assert!(x < self.cols && y < self.rows);
        self.data[y * self.cols + x] = value;
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn width(&self) -> usize {
        self.cols
    }

    pub fn height(&self) -> usize {
        self.rows
    }

    pub fn write_ppm(&self, filename: &str) -> std::io::Result<()> {
        use std::fs::File;
        use std::io::Write;
        let mut file = File::create(filename)?;
        file.write_all(self.to_string().as_bytes())?;
        Ok(())
    }

    pub fn to_buffer(&self) -> Vec<u8> {
        let mut buffer = Vec::with_capacity(self.cols * self.rows * 4);
        for col in 0..self.cols {
            for row in 0..self.rows {
                let color = self.pixel_at(col, row);
                buffer.push((color.r() * 255.0).round() as u8);
                buffer.push((color.g() * 255.0).round() as u8);
                buffer.push((color.b() * 255.0).round() as u8);
                buffer.push(255);
            }
        }
        buffer
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
                let color = self.pixel_at(col, row);
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

// Feature: Canvas
#[cfg(test)]
mod testing {
    use super::*;

    // Scenario: Creating a canvas
    //   Given c ← canvas(10, 20)
    //   Then c.width = 10
    //     And c.height = 20
    //     And every pixel of c is color(0, 0, 0)
    #[test]
    fn creating_canvas() {
        let canvas = Canvas::new(10, 20, Tuple::color(0.0, 0.0, 0.0));
        assert_eq!(canvas.width(), 10);
        assert_eq!(canvas.height(), 20);
    }

    // Scenario: Writing pixels to a canvas
    //   Given c ← canvas(10, 20)
    //     And red ← color(1, 0, 0)
    //   When write_pixel(c, 2, 3, red)
    //   Then pixel_at(c, 2, 3) = red
    #[test]
    fn writing_pixels_to_canvas() {
        let mut canvas = Canvas::new(10, 20, Tuple::color(0.0, 0.0, 0.0));
        let red = Tuple::color(1.0, 0.0, 0.0);
        canvas.write_pixel(3, 2, red.clone());
        assert_eq!(canvas.pixel_at(3, 2), &red);
    }

    // Scenario: Constructing the PPM header
    //   Given c ← canvas(5, 3)
    //   When ppm ← canvas_to_ppm(c)
    //   Then lines 1-3 of ppm are
    //     """
    //     P3
    //     5 3
    //     255
    //     """
    #[test]
    fn constructing_ppm_header() {
        let canvas = Canvas::new(5, 3, Tuple::color(0.0, 0.0, 0.0));
        let ppm = format!("{}", canvas);
        let lines: Vec<&str> = ppm.split('\n').collect();
        assert_eq!(lines[0], "P3");
        assert_eq!(lines[1], "5 3");
        assert_eq!(lines[2], "255");
    }

    // Scenario: Constructing the PPM pixel data
    //   Given c ← canvas(5, 3)
    //     And c1 ← color(1.5, 0, 0)
    //     And c2 ← color(0, 0.5, 0)
    //     And c3 ← color(-0.5, 0, 1)
    //   When write_pixel(c, 0, 0, c1)
    //     And write_pixel(c, 2, 1, c2)
    //     And write_pixel(c, 4, 2, c3)
    //     And ppm ← canvas_to_ppm(c)
    //   Then lines 4-6 of ppm are
    //     """
    //     255 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    //     0 0 0 0 0 0 0 128 0 0 0 0 0 0 0
    //     0 0 0 0 0 0 0 0 0 0 0 0 0 0 255
    //     """
    #[test]
    fn constructing_ppm_pixel_data() {
        let mut canvas = Canvas::new(5, 3, Tuple::color(0.0, 0.0, 0.0));
        let color_a = Tuple::color(1.5, 0.0, 0.0);
        let color_b = Tuple::color(0.0, 0.5, 0.0);
        let color_c = Tuple::color(-0.5, 0.0, 1.0);
        canvas.write_pixel(0, 0, color_a);
        canvas.write_pixel(2, 1, color_b);
        canvas.write_pixel(4, 2, color_c);
        let ppm = format!("{}", canvas);
        let lines: Vec<&str> = ppm.split('\n').collect();
        assert_eq!(lines[3], "255 0 0 0 0 0 0 0 0 0 0 0 0 0 0");
        assert_eq!(lines[4], "0 0 0 0 0 0 0 128 0 0 0 0 0 0 0");
        assert_eq!(lines[5], "0 0 0 0 0 0 0 0 0 0 0 0 0 0 255");
    }

    // Scenario: Splitting long lines in PPM files
    //   Given c ← canvas(10, 2)
    //   When every pixel of c is set to color(1, 0.8, 0.6)
    //     And ppm ← canvas_to_ppm(c)
    //   Then lines 4-7 of ppm are
    //     """
    //     255 204 153 255 204 153 255 204 153 255 204 153 255 204 153 255 204
    //     153 255 204 153 255 204 153 255 204 153 255 204 153
    //     255 204 153 255 204 153 255 204 153 255 204 153 255 204 153 255 204
    //     153 255 204 153 255 204 153 255 204 153 255 204 153
    //     """
    #[test]
    fn splitting_long_lines_in_ppm_files() {
        let canvas = Canvas::new(10, 2, Tuple::color(1.0, 0.8, 0.6));
        let ppm = format!("{}", canvas);
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

    // Scenario: PPM files are terminated by a newline character
    //   Given c ← canvas(5, 3)
    //   When ppm ← canvas_to_ppm(c)
    //   Then ppm ends with a newline character
    #[test]
    fn ppm_files_are_terminated_by_a_newline_character() {
        let canvas = Canvas::new(5, 3, Tuple::color(0.0, 0.0, 0.0));
        let ppm = format!("{}", canvas);
        let lines: Vec<&str> = ppm.split('\n').collect();
        assert_eq!(lines[lines.len() - 1], "");
    }
}
