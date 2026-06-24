use crate::grid::Grid;
use crate::source::Source;
use crate::wave::WaveField;

pub struct Solver {
    pub dt: f32,
    pub dx: f32,
    pub damping_border: usize,
}

impl Solver {
    pub fn new(dt: f32, dx: f32) -> Self {
        Self {
            dt,
            dx,
            damping_border: 80,
        }
    }

    pub fn step(&self, grid: &Grid, field: &mut WaveField, source: &Source, step_index: usize) {
        let w = grid.width;
        let h = grid.height;

        for y in 1..h - 1 {
            for x in 1..w - 1 {
                let idx = field.idx(x, y);

                let center = field.current[idx];

                let laplacian = field.current[field.idx(x - 1, y)]
                    + field.current[field.idx(x + 1, y)]
                    + field.current[field.idx(x, y - 1)]
                    + field.current[field.idx(x, y + 1)]
                    - 4.0 * center;

                let c = grid.get(x, y).speed;
                let coeff = (c * self.dt / self.dx).powi(2);

                field.next[idx] =
                    2.0 * field.current[idx] - field.previous[idx] + coeff * laplacian;
            }
        }

        self.apply_simple_absorbing_boundary(field);

        let value = source.sample(step_index as f32);
        field.inject(source.x, source.y, value);

        field.swap_buffers();
    }

    fn apply_simple_absorbing_boundary(&self, field: &mut WaveField) {
        let w = field.width;
        let h = field.height;
        let border = self.damping_border;
        let absorb_strength = 0.025;

        for y in 0..h {
            for x in 0..w {
                let d = x.min(w - 1 - x).min(y).min(h - 1 - y);

                if d < border {
                    let k = d as f32 / border as f32;
                    let edge_factor = 1.0 - k;
                    let damping = (-absorb_strength * edge_factor * edge_factor).exp();

                    let idx = field.idx(x, y);
                    field.next[idx] *= damping;
                    field.current[idx] *= damping;
                    field.previous[idx] *= damping;
                }
            }
        }
    }
}
