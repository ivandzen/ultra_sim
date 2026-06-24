pub struct Source {
    pub x: usize,
    pub y: usize,
    pub amplitude: f32,
    pub frequency: f32,
    pub delay: f32,
    pub width: f32,
}

impl Source {
    pub fn new(x: usize, y: usize) -> Self {
        Self {
            x,
            y,
            amplitude: 1.0,
            frequency: 0.04,
            delay: 40.0,
            width: 15.0,
        }
    }

    pub fn sample(&self, t: f32) -> f32 {
        let envelope = -((t - self.delay).powi(2)) / (2.0 * self.width * self.width);
        let carrier = (2.0 * std::f32::consts::PI * self.frequency * t).sin();

        self.amplitude * carrier * envelope.exp()
    }
}