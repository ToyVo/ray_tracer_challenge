use nalgebra_glm::DVec3;

pub struct Canvas {
    data: Vec<DVec3>,
    cols: usize,
    rows: usize,
}

impl Canvas {
    pub fn new(cols: usize, rows: usize, init: DVec3) -> Canvas {
        let data = vec![init; cols * rows];
        Canvas { data, cols, rows }
    }

    pub fn from_vec(cols: usize, rows: usize, data: Vec<DVec3>) -> Canvas {
        Canvas { data, cols, rows }
    }

    pub fn pixel_at(&self, x: usize, y: usize) -> &DVec3 {
        assert!(x < self.cols && y < self.rows);
        &self.data[y * self.cols + x]
    }

    pub fn write_pixel(&mut self, x: usize, y: usize, value: DVec3) {
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
                buffer.push((color.x * 255.0).round() as u8);
                buffer.push((color.y * 255.0).round() as u8);
                buffer.push((color.z * 255.0).round() as u8);
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
                row_colors.push((color.x * 255.0).round() as u8);
                row_colors.push((color.y * 255.0).round() as u8);
                row_colors.push((color.z * 255.0).round() as u8);
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
    use nalgebra_glm::vec3;

    #[test]
    fn creating_canvas() {
        let canvas = Canvas::new(10, 20, vec3(0.0, 0.0, 0.0));
        assert_eq!(canvas.width(), 10);
        assert_eq!(canvas.height(), 20);
    }

    #[test]
    fn writing_pixels_to_canvas() {
        let mut canvas = Canvas::new(10, 20, vec3(0.0, 0.0, 0.0));
        let red = vec3(1.0, 0.0, 0.0);
        canvas.write_pixel(3, 2, red);
        assert_eq!(canvas.pixel_at(3, 2), &red);
    }

    #[test]
    fn constructing_ppm_header() {
        let canvas = Canvas::new(5, 3, vec3(0.0, 0.0, 0.0));
        let ppm = format!("{}", canvas);
        let lines: Vec<&str> = ppm.split('\n').collect();
        assert_eq!(lines[0], "P3");
        assert_eq!(lines[1], "5 3");
        assert_eq!(lines[2], "255");
    }

    #[test]
    fn constructing_ppm_pixel_data() {
        let mut canvas = Canvas::new(5, 3, vec3(0.0, 0.0, 0.0));
        let color_a = vec3(1.5, 0.0, 0.0);
        let color_b = vec3(0.0, 0.5, 0.0);
        let color_c = vec3(-0.5, 0.0, 1.0);
        canvas.write_pixel(0, 0, color_a);
        canvas.write_pixel(2, 1, color_b);
        canvas.write_pixel(4, 2, color_c);
        let ppm = format!("{}", canvas);
        let lines: Vec<&str> = ppm.split('\n').collect();
        assert_eq!(lines[3], "255 0 0 0 0 0 0 0 0 0 0 0 0 0 0");
        assert_eq!(lines[4], "0 0 0 0 0 0 0 128 0 0 0 0 0 0 0");
        assert_eq!(lines[5], "0 0 0 0 0 0 0 0 0 0 0 0 0 0 255");
    }

    #[test]
    fn splitting_long_lines_in_ppm_files() {
        let canvas = Canvas::new(10, 2, vec3(1.0, 0.8, 0.6));
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

    #[test]
    fn ppm_files_are_terminated_by_a_newline_character() {
        let canvas = Canvas::new(5, 3, vec3(0.0, 0.0, 0.0));
        let ppm = format!("{}", canvas);
        let lines: Vec<&str> = ppm.split('\n').collect();
        assert_eq!(lines[lines.len() - 1], "");
    }
}
