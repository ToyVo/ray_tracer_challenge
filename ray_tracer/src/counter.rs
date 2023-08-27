pub struct Counter {
    pub count: u32,
}

impl Counter {
    pub fn new() -> Self {
        Self { count: 0 }
    }

    pub fn next(&mut self) -> u32 {
        let id = self.count;
        self.count += 1;
        id
    }
}
