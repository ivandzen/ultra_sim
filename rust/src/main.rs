mod grid;
mod material;
mod sensor;
mod solver;
mod source;
mod vtk_export;
mod wave;

use std::io::{BufWriter, Write};
use std::path::Path;

use grid::Grid;
use material::Material;
use sensor::Sensor;
use solver::Solver;
use source::Source;
use vtk_export::{write_frames_pvd, write_sensors_vtp, write_source_vtp, write_vti_frame};
use wave::WaveField;

fn main() {
    let width = 1024;
    let height = 1024;

    let mut grid = Grid::new(width, height, Material::water());

    grid.add_circle(width / 2, height / 2, 35, Material::bone());

    let mut field = WaveField::new(width, height);

    let source = Source::new(width / 2, 300);

    let mut sensors = vec![
        Sensor::new(width / 2, height - 40),
        Sensor::new(width / 2 - 40, height - 50),
        Sensor::new(width / 2 + 40, height - 50),
    ];

    let solver = Solver::new(0.25, 1.0);

    let total_steps = 1500;
    let run_dir = Path::new("output/run_001");
    std::fs::create_dir_all(run_dir).unwrap();

    let mut saved_frames = Vec::new();

    for step in 0..total_steps {
        solver.step(&grid, &mut field, &source, step);

        for sensor in &mut sensors {
            sensor.record(&field);
        }

        if step % 10 == 0 {
            let file_name = format!("frame_{:04}.vti", step);
            let path = run_dir.join(&file_name);
            write_vti_frame(&path, &grid, &field).unwrap();
            saved_frames.push((step, file_name));
            println!("saved {}", path.display());
        }
    }

    for (i, sensor) in sensors.iter().enumerate() {
        let path = run_dir.join(format!("sensor_{i}.csv"));
        write_sensor_csv(&path, &sensor.samples).unwrap();
    }

    write_frames_pvd(&run_dir.join("frames.pvd"), &saved_frames).unwrap();
    write_sensors_vtp(&run_dir.join("sensors.vtp"), &sensors).unwrap();
    write_source_vtp(&run_dir.join("source.vtp"), &source).unwrap();
}

fn write_sensor_csv(path: &Path, samples: &[f32]) -> std::io::Result<()> {
    let file = std::fs::File::create(path)?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "step,value")?;

    for (i, value) in samples.iter().enumerate() {
        writeln!(writer, "{},{}", i, value)?;
    }

    Ok(())
}
