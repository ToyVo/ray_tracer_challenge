#[derive(Default)]
pub struct Counter {
    pub count: u32,
}

impl Counter {
    pub fn increment(&mut self) -> u32 {
        let id = self.count;
        self.count += 1;
        id
    }
}
