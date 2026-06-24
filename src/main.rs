mod grid;
mod material;
mod sensor;
mod solver;
mod source;
mod wave;

use std::fs::File;
use std::io::{BufWriter, Write};

use grid::Grid;
use material::Material;
use sensor::Sensor;
use solver::Solver;
use source::Source;
use wave::WaveField;

fn main() {
    let width = 256;
    let height = 256;

    let mut grid = Grid::new(width, height, Material::water());

    grid.add_circle(
        width / 2,
        height / 2,
        35,
        Material::muscle(),
    );

    let mut field = WaveField::new(width, height);

    let source = Source::new(width / 2, 40);

    let mut sensors = vec![
        Sensor::new(width / 2, height - 40),
        Sensor::new(width / 2 - 40, height - 50),
        Sensor::new(width / 2 + 40, height - 50),
    ];

    let solver = Solver::new(0.4, 1.0);

    let total_steps = 800;

    std::fs::create_dir_all("frames").unwrap();

    for step in 0..total_steps {
        solver.step(&grid, &mut field, &source, step);

        for sensor in &mut sensors {
            sensor.record(&field);
        }

        if step % 10 == 0 {
            let path = format!("frames/frame_{:04}.pgm", step);
            write_pgm(&path, &field.current, width, height).unwrap();
            println!("saved {path}");
        }
    }

    for (i, sensor) in sensors.iter().enumerate() {
        let path = format!("sensor_{i}.csv");
        write_sensor_csv(&path, &sensor.samples).unwrap();
    }
}

fn write_pgm(path: &str, data: &[f32], width: usize, height: usize) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "P2")?;
    writeln!(writer, "{} {}", width, height)?;
    writeln!(writer, "255")?;

    let max_abs = data
        .iter()
        .fold(0.0_f32, |m, v| m.max(v.abs()))
        .max(1e-6);

    for y in 0..height {
        for x in 0..width {
            let v = data[y * width + x];
            let normalized = ((v / max_abs) * 0.5 + 0.5).clamp(0.0, 1.0);
            let pixel = (normalized * 255.0) as u8;
            write!(writer, "{} ", pixel)?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

fn write_sensor_csv(path: &str, samples: &[f32]) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "step,value")?;

    for (i, value) in samples.iter().enumerate() {
        writeln!(writer, "{},{}", i, value)?;
    }

    Ok(())
}