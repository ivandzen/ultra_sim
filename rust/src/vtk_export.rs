use crate::grid::Grid;
use crate::material::MaterialKind;
use crate::sensor::Sensor;
use crate::source::Source;
use crate::wave::WaveField;

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

fn material_id(kind: MaterialKind) -> i32 {
    match kind {
        MaterialKind::Water => 0,
        MaterialKind::Fat => 1,
        MaterialKind::Muscle => 2,
        MaterialKind::Bone => 3,
        MaterialKind::Custom => 4,
    }
}

fn write_vtk_open_tag(writer: &mut dyn Write, file_type: &str) -> std::io::Result<()> {
    writeln!(
        writer,
        r#"<VTKFile type="{file_type}" version="0.1" byte_order="LittleEndian">"#
    )
}

fn write_vtk_close_tag(writer: &mut dyn Write) -> std::io::Result<()> {
    writeln!(writer, "</VTKFile>")
}

fn write_ascii_data_array<T: std::fmt::Display>(
    writer: &mut dyn Write,
    name: &str,
    data_type: &str,
    num_components: usize,
    values: &[T],
) -> std::io::Result<()> {
    if num_components == 1 {
        writeln!(
            writer,
            r#"<DataArray type="{data_type}" Name="{name}" format="ascii">"#
        )?;
    } else {
        writeln!(
            writer,
            r#"<DataArray type="{data_type}" Name="{name}" NumberOfComponents="{num_components}" format="ascii">"#
        )?;
    }

    for (i, value) in values.iter().enumerate() {
        write!(writer, "{value}")?;
        if i + 1 != values.len() {
            write!(writer, " ")?;
        }
    }
    writeln!(writer)?;
    writeln!(writer, "</DataArray>")
}

pub fn write_vti_frame(path: &Path, grid: &Grid, field: &WaveField) -> std::io::Result<()> {
    if grid.width != field.width || grid.height != field.height {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "grid and field dimensions do not match",
        ));
    }
    if grid.cells.len() != field.pressure.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "grid and field sizes do not match",
        ));
    }

    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    write_vtk_open_tag(&mut writer, "ImageData")?;
    writeln!(
        writer,
        r#"  <ImageData WholeExtent="0 {} 0 {} 0 0" Origin="0 0 0" Spacing="1 1 1">"#,
        grid.width.saturating_sub(1),
        grid.height.saturating_sub(1)
    )?;
    writeln!(
        writer,
        r#"    <Piece Extent="0 {} 0 {} 0 0">"#,
        grid.width.saturating_sub(1),
        grid.height.saturating_sub(1)
    )?;
    writeln!(writer, r#"      <PointData Scalars="pressure">"#)?;

    write_ascii_data_array(&mut writer, "pressure", "Float32", 1, &field.pressure)?;

    let speed: Vec<f32> = grid.cells.iter().map(|cell| cell.speed).collect();
    let density: Vec<f32> = grid.cells.iter().map(|cell| cell.density).collect();
    let attenuation: Vec<f32> = grid.cells.iter().map(|cell| cell.attenuation).collect();
    let material_ids: Vec<i32> = grid
        .cells
        .iter()
        .map(|cell| material_id(cell.kind))
        .collect();

    write_ascii_data_array(&mut writer, "speed", "Float32", 1, &speed)?;
    write_ascii_data_array(&mut writer, "density", "Float32", 1, &density)?;
    write_ascii_data_array(&mut writer, "attenuation", "Float32", 1, &attenuation)?;
    write_ascii_data_array(&mut writer, "material_id", "Int32", 1, &material_ids)?;

    writeln!(writer, "      </PointData>")?;
    writeln!(writer, "      <CellData/>")?;
    writeln!(writer, "    </Piece>")?;
    writeln!(writer, "  </ImageData>")?;
    write_vtk_close_tag(&mut writer)
}

pub fn write_frames_pvd(path: &Path, frames: &[(usize, String)]) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    write_vtk_open_tag(&mut writer, "Collection")?;
    writeln!(writer, "  <Collection>")?;

    for (step, file_name) in frames {
        writeln!(
            writer,
            r#"    <DataSet timestep="{step}" file="{file_name}"/>"#
        )?;
    }

    writeln!(writer, "  </Collection>")?;
    write_vtk_close_tag(&mut writer)
}

pub fn write_sensors_vtp(path: &Path, sensors: &[Sensor]) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    write_vtk_open_tag(&mut writer, "PolyData")?;
    writeln!(writer, "  <PolyData>")?;
    writeln!(
        writer,
        r#"    <Piece NumberOfPoints="{}" NumberOfVerts="{}" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">"#,
        sensors.len(),
        sensors.len()
    )?;
    writeln!(writer, r#"      <PointData Scalars="sensor_id">"#)?;

    let sensor_ids: Vec<i32> = (0..sensors.len()).map(|id| id as i32).collect();
    write_ascii_data_array(&mut writer, "sensor_id", "Int32", 1, &sensor_ids)?;
    writeln!(writer, "      </PointData>")?;

    writeln!(writer, "      <Points>")?;
    let mut coords = Vec::with_capacity(sensors.len() * 3);
    for sensor in sensors {
        coords.push(sensor.x as f32);
        coords.push(sensor.y as f32);
        coords.push(0.0);
    }
    write_ascii_data_array(&mut writer, "Points", "Float32", 3, &coords)?;
    writeln!(writer, "      </Points>")?;

    writeln!(writer, "      <Verts>")?;
    let connectivity: Vec<i32> = (0..sensors.len()).map(|i| i as i32).collect();
    let offsets: Vec<i32> = (1..=sensors.len()).map(|i| i as i32).collect();
    write_ascii_data_array(&mut writer, "connectivity", "Int32", 1, &connectivity)?;
    write_ascii_data_array(&mut writer, "offsets", "Int32", 1, &offsets)?;
    writeln!(writer, "      </Verts>")?;

    writeln!(writer, "    </Piece>")?;
    writeln!(writer, "  </PolyData>")?;
    write_vtk_close_tag(&mut writer)
}

pub fn write_source_vtp(path: &Path, source: &Source) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    write_vtk_open_tag(&mut writer, "PolyData")?;
    writeln!(writer, "  <PolyData>")?;
    writeln!(
        writer,
        r#"    <Piece NumberOfPoints="1" NumberOfVerts="1" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">"#
    )?;
    writeln!(writer, "      <PointData/>")?;
    writeln!(writer, "      <Points>")?;
    let coords = [source.x as f32, source.y as f32, 0.0];
    write_ascii_data_array(&mut writer, "Points", "Float32", 3, &coords)?;
    writeln!(writer, "      </Points>")?;
    writeln!(writer, "      <Verts>")?;
    let connectivity = [0_i32];
    let offsets = [1_i32];
    write_ascii_data_array(&mut writer, "connectivity", "Int32", 1, &connectivity)?;
    write_ascii_data_array(&mut writer, "offsets", "Int32", 1, &offsets)?;
    writeln!(writer, "      </Verts>")?;
    writeln!(writer, "    </Piece>")?;
    writeln!(writer, "  </PolyData>")?;
    write_vtk_close_tag(&mut writer)
}
