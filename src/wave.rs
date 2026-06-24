pub struct WaveField {
    pub width: usize,
    pub height: usize,
    pub previous: Vec<f32>,
    pub current: Vec<f32>,
    pub next: Vec<f32>,
}

impl WaveField {
    pub fn new(width: usize, height: usize) -> Self {
        let size = width * height;

        Self {
            width,
            height,
            previous: vec![0.0; size],
            current: vec![0.0; size],
            next: vec![0.0; size],
        }
    }

    #[inline]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.width + x
    }

    pub fn inject(&mut self, x: usize, y: usize, value: f32) {
        let idx = self.idx(x, y);
        self.current[idx] += value;
    }

    pub fn swap_buffers(&mut self) {
        std::mem::swap(&mut self.previous, &mut self.current);
        std::mem::swap(&mut self.current, &mut self.next);
        self.next.fill(0.0);
    }
}
