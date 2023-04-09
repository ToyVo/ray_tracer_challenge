use crate::color::Color;

pub struct Canvas {
    pub width: usize,
    pub height: usize,
    pub pixels: Vec<Color>,
}

impl Canvas {
    pub fn new(width: usize, height: usize) -> Canvas {
        Canvas {
            width,
            height,
            pixels: vec![Color::new(0.0, 0.0, 0.0); (width * height) as usize],
        }
    }

    pub fn write_pixel(&mut self, x: usize, y: usize, color: Color) {
        if x < self.width && y < self.height {
            let index = (y * self.width + x) as usize;
            self.pixels[index] = color;
        }
    }

    pub fn pixel_at(&self, x: usize, y: usize) -> Color {
        let index = (y * self.width + x) as usize;
        self.pixels[index]
    }
}

impl std::fmt::Display for Canvas {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut ppm = String::with_capacity((self.width * self.height * 12) as usize);
        ppm.push_str("P3\n");
        ppm.push_str(&format!("{} {}\n", self.width, self.height));
        ppm.push_str("255\n");

        for y in 0..self.height {
            let mut row_colors: Vec<u8> = Vec::with_capacity((self.width * 3) as usize);
            for x in 0..self.width {
                let color = self.pixel_at(x, y);
                row_colors.push((color.red * 255.0).round() as u8);
                row_colors.push((color.green * 255.0).round() as u8);
                row_colors.push((color.blue * 255.0).round() as u8);
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
            ppm.push_str(&row.trim_end());
            ppm.push_str("\n");
        }
        write!(f, "{}", ppm)
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn creating_canvas() {
        let c = Canvas::new(10, 20);
        assert_eq!(c.width, 10);
        assert_eq!(c.height, 20);
        assert_eq!(c.pixels.len(), 200);
        for pixel in c.pixels {
            assert_eq!(pixel, Color::new(0.0, 0.0, 0.0));
        }
    }

    #[test]
    fn writing_pixels_to_canvas() {
        let mut c = Canvas::new(10, 20);
        let red = Color::new(1.0, 0.0, 0.0);
        c.write_pixel(2, 3, red);
        assert_eq!(c.pixel_at(2, 3), red);
    }

    #[test]
    fn constructing_ppm_header() {
        let c = Canvas::new(5, 3);
        let ppm = format!("{}", c);
        let lines: Vec<&str> = ppm.split("\n").collect();
        assert_eq!(lines[0], "P3");
        assert_eq!(lines[1], "5 3");
        assert_eq!(lines[2], "255");
    }

    #[test]
    fn constructing_ppm_pixel_data() {
        let mut c = Canvas::new(5, 3);
        let c1 = Color::new(1.5, 0.0, 0.0);
        let c2 = Color::new(0.0, 0.5, 0.0);
        let c3 = Color::new(-0.5, 0.0, 1.0);
        c.write_pixel(0, 0, c1);
        c.write_pixel(2, 1, c2);
        c.write_pixel(4, 2, c3);
        let ppm = format!("{}", c);
        let lines: Vec<&str> = ppm.split("\n").collect();
        assert_eq!(lines[3], "255 0 0 0 0 0 0 0 0 0 0 0 0 0 0");
        assert_eq!(lines[4], "0 0 0 0 0 0 0 128 0 0 0 0 0 0 0");
        assert_eq!(lines[5], "0 0 0 0 0 0 0 0 0 0 0 0 0 0 255");
    }

    #[test]
        fn splitting_long_lines_in_ppm_files() {
        let mut c = Canvas::new(10, 2);
        let color = Color::new(1.0, 0.8, 0.6);
        for y in 0..2 {
            for x in 0..10 {
                c.write_pixel(x, y, color);
            }
        }
        let ppm = format!("{}", c);
        let lines: Vec<&str> = ppm.split("\n").collect();
        assert_eq!(lines[3], "255 204 153 255 204 153 255 204 153 255 204 153 255 204 153 255 204");
        assert_eq!(lines[4], "153 255 204 153 255 204 153 255 204 153 255 204 153");
        assert_eq!(lines[5], "255 204 153 255 204 153 255 204 153 255 204 153 255 204 153 255 204");
        assert_eq!(lines[6], "153 255 204 153 255 204 153 255 204 153 255 204 153");
    }

    #[test]
        fn ppm_files_are_terminated_by_a_newline_character() {
        let c = Canvas::new(5, 3);
        let ppm = format!("{}", c);
        let lines: Vec<&str> = ppm.split("\n").collect();
        assert_eq!(lines[lines.len() - 1], "");
    }
}
