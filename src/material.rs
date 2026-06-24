#[derive(Debug, Clone, Copy)]
pub struct Material {
    pub speed: f32, // relative speed of sound, water = 1.0
}

impl Material {
    pub fn water() -> Self {
        Self { speed: 1.0 }
    }

    pub fn fat() -> Self {
        Self { speed: 0.95 }
    }

    pub fn muscle() -> Self {
        Self { speed: 1.05 }
    }

    pub fn bone() -> Self {
        Self { speed: 2.0 }
    }
}