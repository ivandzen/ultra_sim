/// Staggered acoustic field storage.
///
/// Pressure lives at cell centers:
///   p[i, j] -> width * height
///
/// Velocity components live on cell faces:
///   vx[i + 1/2, j] -> (width + 1) * height
///   vy[i, j + 1/2] -> width * (height + 1)
pub struct WaveField {
    pub width: usize,
    pub height: usize,
    pub pressure: Vec<f32>,
    pub pressure_next: Vec<f32>,
    pub vx: Vec<f32>,
    pub vy: Vec<f32>,
}

impl WaveField {
    pub fn new(width: usize, height: usize) -> Self {
        let pressure_size = width * height;
        let vx_size = (width + 1) * height;
        let vy_size = width * (height + 1);

        Self {
            width,
            height,
            pressure: vec![0.0; pressure_size],
            pressure_next: vec![0.0; pressure_size],
            vx: vec![0.0; vx_size],
            vy: vec![0.0; vy_size],
        }
    }

    #[inline]
    pub fn pressure_idx(&self, x: usize, y: usize) -> usize {
        y * self.width + x
    }

    #[inline]
    pub fn vx_idx(&self, x_face: usize, y: usize) -> usize {
        y * (self.width + 1) + x_face
    }

    #[inline]
    pub fn vy_idx(&self, x: usize, y_face: usize) -> usize {
        y_face * self.width + x
    }

    #[allow(dead_code)]
    #[inline]
    pub fn inject_pressure(&mut self, x: usize, y: usize, value: f32) {
        let idx = self.pressure_idx(x, y);
        self.pressure[idx] += value;
    }

    #[inline]
    pub fn pressure_at(&self, x: usize, y: usize) -> f32 {
        self.pressure[self.pressure_idx(x, y)]
    }

    pub fn swap_pressure_buffers(&mut self) {
        std::mem::swap(&mut self.pressure, &mut self.pressure_next);
        self.pressure_next.fill(0.0);
    }
}
