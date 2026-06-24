use crate::wave::WaveField;

pub struct Sensor {
    pub x: usize,
    pub y: usize,
    pub samples: Vec<f32>,
}

impl Sensor {
    pub fn new(x: usize, y: usize) -> Self {
        Self {
            x,
            y,
            samples: Vec::new(),
        }
    }

    pub fn record(&mut self, field: &WaveField) {
        let idx = field.idx(self.x, self.y);
        self.samples.push(field.current[idx]);
    }
}
