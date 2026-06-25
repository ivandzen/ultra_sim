use crate::grid::Grid;
use crate::material::Material;
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
        debug_assert!(
            self.cfl_number(grid) <= 1.0 / 2.0_f32.sqrt(),
            "CFL condition violated: c_max * dt / dx must be <= 1/sqrt(2)"
        );

        self.update_velocities(grid, field);
        self.update_pressure(grid, field);

        let source_idx = field.pressure_idx(source.x, source.y);
        field.pressure_next[source_idx] += source.sample(step_index as f32);

        self.apply_simple_absorbing_boundary(field);
        field.swap_pressure_buffers();
    }

    fn update_velocities(&self, grid: &Grid, field: &mut WaveField) {
        let w = grid.width;
        let h = grid.height;

        let dt_over_dx = self.dt / self.dx;

        if w >= 2 && h > 0 {
            for y in 0..h {
                for x_face in 1..w {
                    let left = grid.get(x_face - 1, y);
                    let right = grid.get(x_face, y);
                    let rho = harmonic(left.density, right.density);

                    if rho <= f32::EPSILON {
                        continue;
                    }

                    let p_right = field.pressure[field.pressure_idx(x_face, y)];
                    let p_left = field.pressure[field.pressure_idx(x_face - 1, y)];
                    let idx = field.vx_idx(x_face, y);

                    field.vx[idx] -= dt_over_dx / rho * (p_right - p_left);
                }
            }
        }

        if h >= 2 && w > 0 {
            for y_face in 1..h {
                for x in 0..w {
                    let upper = grid.get(x, y_face - 1);
                    let lower = grid.get(x, y_face);
                    let rho = harmonic(upper.density, lower.density);

                    if rho <= f32::EPSILON {
                        continue;
                    }

                    let p_down = field.pressure[field.pressure_idx(x, y_face)];
                    let p_up = field.pressure[field.pressure_idx(x, y_face - 1)];
                    let idx = field.vy_idx(x, y_face);

                    field.vy[idx] -= dt_over_dx / rho * (p_down - p_up);
                }
            }
        }
    }

    fn update_pressure(&self, grid: &Grid, field: &mut WaveField) {
        let w = grid.width;
        let h = grid.height;

        if w == 0 || h == 0 {
            return;
        }

        for y in 0..h {
            for x in 0..w {
                let material = grid.get(x, y);
                let idx = field.pressure_idx(x, y);

                let vx_right = if x + 1 < field.width + 1 {
                    field.vx[field.vx_idx(x + 1, y)]
                } else {
                    0.0
                };
                let vx_left = if x > 0 {
                    field.vx[field.vx_idx(x, y)]
                } else {
                    0.0
                };
                let vy_down = if y + 1 < field.height + 1 {
                    field.vy[field.vy_idx(x, y + 1)]
                } else {
                    0.0
                };
                let vy_up = if y > 0 {
                    field.vy[field.vy_idx(x, y)]
                } else {
                    0.0
                };

                let divergence = ((vx_right - vx_left) + (vy_down - vy_up)) / self.dx;
                let attenuation = attenuation_factor(&material, self.dt);

                field.pressure_next[idx] = (field.pressure[idx]
                    - material.bulk_modulus() * self.dt * divergence)
                    * attenuation;
            }
        }
    }

    fn cfl_number(&self, grid: &Grid) -> f32 {
        let c_max = grid
            .cells
            .iter()
            .map(|material| material.speed)
            .fold(0.0_f32, f32::max);

        c_max * self.dt / self.dx
    }

    fn apply_simple_absorbing_boundary(&self, field: &mut WaveField) {
        let border = self.damping_border;
        let absorb_strength = 0.025;

        apply_sponge(
            &mut field.pressure_next,
            field.width,
            field.height,
            border,
            absorb_strength,
        );
        apply_sponge(
            &mut field.vx,
            field.width + 1,
            field.height,
            border,
            absorb_strength,
        );
        apply_sponge(
            &mut field.vy,
            field.width,
            field.height + 1,
            border,
            absorb_strength,
        );
    }
}

fn harmonic(a: f32, b: f32) -> f32 {
    let denom = a + b;

    if denom.abs() <= f32::EPSILON {
        0.0
    } else {
        2.0 * a * b / denom
    }
}

fn attenuation_factor(material: &Material, dt: f32) -> f32 {
    (-material.attenuation * dt).exp().clamp(0.0, 1.0)
}

fn apply_sponge(
    buffer: &mut [f32],
    width: usize,
    height: usize,
    border: usize,
    absorb_strength: f32,
) {
    if width == 0 || height == 0 || border == 0 {
        return;
    }

    for y in 0..height {
        for x in 0..width {
            let distance_to_edge = x.min(width - 1 - x).min(y).min(height - 1 - y);

            if distance_to_edge < border {
                let k = distance_to_edge as f32 / border as f32;
                let edge_factor = 1.0 - k;
                let damping = (-absorb_strength * edge_factor * edge_factor).exp();
                let idx = y * width + x;

                buffer[idx] *= damping;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::Material;

    fn make_solver() -> Solver {
        Solver {
            dt: 0.25,
            dx: 1.0,
            damping_border: 4,
        }
    }

    fn centered_source(x: usize, y: usize) -> Source {
        Source {
            x,
            y,
            amplitude: 1.0,
            frequency: 0.12,
            delay: 0.0,
            width: 20.0,
        }
    }

    fn max_abs(values: &[f32]) -> f32 {
        values.iter().map(|v| v.abs()).fold(0.0, f32::max)
    }

    #[test]
    fn homogeneous_water_remains_symmetric_and_stable() {
        let width = 64;
        let height = 64;
        let grid = Grid::new(width, height, Material::water());
        let mut field = WaveField::new(width, height);
        let solver = make_solver();
        let source = centered_source(width / 2, height / 2);

        for step in 0..80 {
            solver.step(&grid, &mut field, &source, step);
        }

        let c = width / 2;
        let r = 8;
        let p_left = field.pressure_at(c - r, c);
        let p_right = field.pressure_at(c + r, c);
        let p_up = field.pressure_at(c, c - r);
        let p_down = field.pressure_at(c, c + r);

        assert!(p_left.is_finite());
        assert!(p_right.is_finite());
        assert!(p_up.is_finite());
        assert!(p_down.is_finite());
        assert!((p_left - p_right).abs() < 5e-2);
        assert!((p_up - p_down).abs() < 5e-2);
        assert!(max_abs(&field.pressure) < 10.0);
    }

    #[test]
    fn constant_density_field_propagates_without_instability() {
        let width = 48;
        let height = 48;
        let mut grid = Grid::new(width, height, Material::water());
        for y in 0..height {
            for x in 0..width {
                grid.set(x, y, Material::custom(1.0, 1.0, 0.0));
            }
        }

        let mut field = WaveField::new(width, height);
        let solver = make_solver();
        let source = centered_source(width / 2, height / 2);

        for step in 0..100 {
            solver.step(&grid, &mut field, &source, step);
        }

        assert!(field.pressure.iter().all(|value| value.is_finite()));
        assert!(field.vx.iter().all(|value| value.is_finite()));
        assert!(field.vy.iter().all(|value| value.is_finite()));
        assert!(max_abs(&field.pressure) > 0.0);
    }

    #[test]
    fn bone_disk_reflects_and_shadows_before_transmission() {
        let width = 80;
        let height = 80;
        let source = centered_source(width / 2, 8);
        let solver = make_solver();

        let water_grid = Grid::new(width, height, Material::water());
        let mut bone_grid = Grid::new(width, height, Material::water());
        bone_grid.add_circle(width / 2, height / 2, 10, Material::bone());

        let mut water_field = WaveField::new(width, height);
        let mut bone_field = WaveField::new(width, height);

        let shadow_point = (width / 2, height - 24);
        let bone_interior = (width / 2, height / 2);
        let mut water_shadow_peak = 0.0_f32;
        let mut bone_shadow_peak = 0.0_f32;

        for step in 0..220 {
            solver.step(&water_grid, &mut water_field, &source, step);
            solver.step(&bone_grid, &mut bone_field, &source, step);

            water_shadow_peak = water_shadow_peak.max(
                water_field
                    .pressure_at(shadow_point.0, shadow_point.1)
                    .abs(),
            );
            bone_shadow_peak =
                bone_shadow_peak.max(bone_field.pressure_at(shadow_point.0, shadow_point.1).abs());

            if step < 60 {
                assert!(
                    bone_field
                        .pressure_at(bone_interior.0, bone_interior.1)
                        .abs()
                        < 1e-4,
                    "pressure reached the bone interior before the wavefront"
                );
            }
        }

        let bone_center = bone_field
            .pressure_at(bone_interior.0, bone_interior.1)
            .abs();

        assert!(bone_shadow_peak < water_shadow_peak * 0.7);
        assert!(bone_center > 0.0);
    }

    #[test]
    #[should_panic(expected = "CFL condition violated")]
    fn cfl_violation_triggers_debug_assertion() {
        let grid = Grid::new(8, 8, Material::custom(1.0, 2.0, 0.0));
        let mut field = WaveField::new(8, 8);
        let solver = Solver::new(1.0, 1.0);
        let source = centered_source(4, 4);

        solver.step(&grid, &mut field, &source, 0);
    }
}
