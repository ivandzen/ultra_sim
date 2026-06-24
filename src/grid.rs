use crate::material::Material;

pub struct Grid {
    pub width: usize,
    pub height: usize,
    pub cells: Vec<Material>,
}

impl Grid {
    pub fn new(width: usize, height: usize, material: Material) -> Self {
        Self {
            width,
            height,
            cells: vec![material; width * height],
        }
    }

    #[inline]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.width + x
    }

    #[inline]
    pub fn get(&self, x: usize, y: usize) -> Material {
        self.cells[self.idx(x, y)]
    }

    pub fn set(&mut self, x: usize, y: usize, material: Material) {
        let idx = self.idx(x, y);
        self.cells[idx] = material;
    }

    pub fn add_circle(&mut self, cx: usize, cy: usize, radius: usize, material: Material) {
        let r2 = (radius * radius) as isize;

        for y in 0..self.height {
            for x in 0..self.width {
                let dx = x as isize - cx as isize;
                let dy = y as isize - cy as isize;

                if dx * dx + dy * dy <= r2 {
                    self.set(x, y, material);
                }
            }
        }
    }
}