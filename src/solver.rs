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
        // For visual stability with higher-order stencils near heterogeneous boundaries,
        // keep c_max * dt / dx <= 0.45 when possible.
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

        if w >= 2 && h > 0 {
            for y in 0..h {
                for x_face in 1..w {
                    let left = grid.get(x_face - 1, y);
                    let right = grid.get(x_face, y);
                    let rho = harmonic(left.density, right.density);

                    if rho <= f32::EPSILON {
                        continue;
                    }

                    let dpdx = pressure_dx_at_vx_face(field, x_face, y, self.dx);
                    let idx = field.vx_idx(x_face, y);

                    field.vx[idx] -= self.dt / rho * dpdx;
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

                    let dpdy = pressure_dy_at_vy_face(field, x, y_face, self.dx);
                    let idx = field.vy_idx(x, y_face);

                    field.vy[idx] -= self.dt / rho * dpdy;
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

                let dvx_dx = vx_dx_at_pressure_cell(field, x, y, self.dx);
                let dvy_dy = vy_dy_at_pressure_cell(field, x, y, self.dx);
                let attenuation = attenuation_factor(&material, self.dt);

                field.pressure_next[idx] = (field.pressure[idx]
                    - material.bulk_modulus() * self.dt * (dvx_dx + dvy_dy))
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

fn staggered_diff_2(right: f32, left: f32, dx: f32) -> f32 {
    (right - left) / dx
}

fn staggered_diff_4(f_plus_1: f32, f_0: f32, f_minus_1: f32, f_minus_2: f32, dx: f32) -> f32 {
    ((9.0 / 8.0) * (f_0 - f_minus_1) - (1.0 / 24.0) * (f_plus_1 - f_minus_2)) / dx
}

fn pressure_dx_at_vx_face(field: &WaveField, x_face: usize, y: usize, dx: f32) -> f32 {
    let right = field.pressure[field.pressure_idx(x_face, y)];
    let left = field.pressure[field.pressure_idx(x_face - 1, y)];

    if x_face >= 2 && x_face + 1 < field.width {
        let plus_1 = field.pressure[field.pressure_idx(x_face + 1, y)];
        let minus_2 = field.pressure[field.pressure_idx(x_face - 2, y)];
        staggered_diff_4(plus_1, right, left, minus_2, dx)
    } else {
        staggered_diff_2(right, left, dx)
    }
}

fn pressure_dy_at_vy_face(field: &WaveField, x: usize, y_face: usize, dx: f32) -> f32 {
    let down = field.pressure[field.pressure_idx(x, y_face)];
    let up = field.pressure[field.pressure_idx(x, y_face - 1)];

    if y_face >= 2 && y_face + 1 < field.height {
        let plus_1 = field.pressure[field.pressure_idx(x, y_face + 1)];
        let minus_2 = field.pressure[field.pressure_idx(x, y_face - 2)];
        staggered_diff_4(plus_1, down, up, minus_2, dx)
    } else {
        staggered_diff_2(down, up, dx)
    }
}

fn vx_dx_at_pressure_cell(field: &WaveField, x: usize, y: usize, dx: f32) -> f32 {
    let right = field.vx[field.vx_idx(x + 1, y)];
    let left = field.vx[field.vx_idx(x, y)];

    if x >= 1 && x + 2 <= field.width {
        let plus_1 = field.vx[field.vx_idx(x + 2, y)];
        let minus_2 = field.vx[field.vx_idx(x - 1, y)];
        staggered_diff_4(plus_1, right, left, minus_2, dx)
    } else {
        staggered_diff_2(right, left, dx)
    }
}

fn vy_dy_at_pressure_cell(field: &WaveField, x: usize, y: usize, dx: f32) -> f32 {
    let down = field.vy[field.vy_idx(x, y + 1)];
    let up = field.vy[field.vy_idx(x, y)];

    if y >= 1 && y + 2 <= field.height {
        let plus_1 = field.vy[field.vy_idx(x, y + 2)];
        let minus_2 = field.vy[field.vy_idx(x, y - 1)];
        staggered_diff_4(plus_1, down, up, minus_2, dx)
    } else {
        staggered_diff_2(down, up, dx)
    }
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

    fn fill_pressure_pattern(field: &mut WaveField) {
        for y in 0..field.height {
            for x in 0..field.width {
                let idx = field.pressure_idx(x, y);
                field.pressure[idx] = (x as f32).powi(3) + 3.0 * (y as f32).powi(2);
            }
        }
    }

    fn fill_face_velocity_pattern(field: &mut WaveField) {
        for y in 0..field.height {
            for x in 0..(field.width + 1) {
                let idx = field.vx_idx(x, y);
                field.vx[idx] = (x as f32).powi(3) + 2.0 * y as f32;
            }
        }

        for y in 0..(field.height + 1) {
            for x in 0..field.width {
                let idx = field.vy_idx(x, y);
                field.vy[idx] = 2.0 * x as f32 + (y as f32).powi(3);
            }
        }
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
    fn fourth_order_stencils_are_used_in_the_bulk_region() {
        let width = 6;
        let height = 6;
        let grid = Grid::new(width, height, Material::water());
        let solver = Solver {
            dt: 1.0,
            dx: 1.0,
            damping_border: 0,
        };

        let mut velocity_field = WaveField::new(width, height);
        fill_pressure_pattern(&mut velocity_field);
        solver.update_velocities(&grid, &mut velocity_field);

        let x_face = 3;
        let y = 2;
        let expected_dpdx = staggered_diff_4(
            velocity_field.pressure[velocity_field.pressure_idx(x_face + 1, y)],
            velocity_field.pressure[velocity_field.pressure_idx(x_face, y)],
            velocity_field.pressure[velocity_field.pressure_idx(x_face - 1, y)],
            velocity_field.pressure[velocity_field.pressure_idx(x_face - 2, y)],
            1.0,
        );
        let expected_vx = -expected_dpdx;
        let vx_idx = velocity_field.vx_idx(x_face, y);
        assert!((velocity_field.vx[vx_idx] - expected_vx).abs() < 1e-6);

        let mut pressure_field = WaveField::new(width, height);
        fill_face_velocity_pattern(&mut pressure_field);
        solver.update_pressure(&grid, &mut pressure_field);

        let x = 2;
        let y = 2;
        let expected_dvx_dx = staggered_diff_4(
            pressure_field.vx[pressure_field.vx_idx(x + 2, y)],
            pressure_field.vx[pressure_field.vx_idx(x + 1, y)],
            pressure_field.vx[pressure_field.vx_idx(x, y)],
            pressure_field.vx[pressure_field.vx_idx(x - 1, y)],
            1.0,
        );
        let expected_dvy_dy = staggered_diff_4(
            pressure_field.vy[pressure_field.vy_idx(x, y + 2)],
            pressure_field.vy[pressure_field.vy_idx(x, y + 1)],
            pressure_field.vy[pressure_field.vy_idx(x, y)],
            pressure_field.vy[pressure_field.vy_idx(x, y - 1)],
            1.0,
        );
        let expected_pressure = -(expected_dvx_dx + expected_dvy_dy);
        let p_idx = pressure_field.pressure_idx(x, y);
        assert!((pressure_field.pressure_next[p_idx] - expected_pressure).abs() < 1e-6);
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
