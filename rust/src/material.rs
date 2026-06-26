#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub enum MaterialKind {
    Water,
    Fat,
    Muscle,
    Bone,
    Custom,
}

#[derive(Debug, Clone, Copy)]
pub struct Material {
    pub kind: MaterialKind,
    pub density: f32,
    pub speed: f32,
    pub attenuation: f32,
}

#[allow(dead_code)]
impl Material {
    pub fn bulk_modulus(&self) -> f32 {
        self.density * self.speed * self.speed
    }

    pub fn impedance(&self) -> f32 {
        self.density * self.speed
    }

    pub fn water() -> Self {
        Self {
            kind: MaterialKind::Water,
            density: 1.0,
            speed: 1.0,
            attenuation: 0.0,
        }
    }

    pub fn fat() -> Self {
        Self {
            kind: MaterialKind::Fat,
            density: 0.92,
            speed: 0.95,
            attenuation: 0.01,
        }
    }

    pub fn muscle() -> Self {
        Self {
            kind: MaterialKind::Muscle,
            density: 1.06,
            speed: 1.05,
            attenuation: 0.02,
        }
    }

    pub fn bone() -> Self {
        Self {
            kind: MaterialKind::Bone,
            density: 1.85,
            speed: 2.4,
            attenuation: 0.2,
        }
    }

    pub fn custom(density: f32, speed: f32, attenuation: f32) -> Self {
        Self {
            kind: MaterialKind::Custom,
            density,
            speed,
            attenuation,
        }
    }
}
