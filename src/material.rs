#[derive(Debug, Clone, Copy)]
pub struct Material {
    pub c: f32, // relative speed of sound, water = 1.0
}

impl Material {
    pub fn water() -> Self {
        Self { c: 1.0 }
    }

    pub fn fat() -> Self {
        Self { c: 0.95 }
    }

    pub fn muscle() -> Self {
        Self { c: 1.05 }
    }

    pub fn bone() -> Self {
        Self { c: 2.0 }
    }
}