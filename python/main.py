import os
from pathlib import Path
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk

from kwave.data import Vector
from kwave.kgrid import kWaveGrid
from kwave.kmedium import kWaveMedium
from kwave.ksource import kSource
from kwave.ksensor import kSensor
from kwave.kspaceFirstOrder import kspaceFirstOrder
from kwave.options.simulation_execution_options import SimulationExecutionOptions
from kwave.options.simulation_options import SimulationOptions
from kwave.utils.mapgen import make_disc
from kwave.utils.signals import tone_burst
from kwave.compat import options_to_kwargs


def export_pressure_timeseries_to_vti(
        pressure_data,
        width,
        height,
        sound_speed,
        density,
        alpha_coeff,
        output_dir,
        frame_step=10,
):
    os.makedirs(output_dir, exist_ok=True)

    pressure_data = np.asarray(pressure_data)
    num_points = width * height

    print(f"pressure_data.shape = {pressure_data.shape}")

    if pressure_data.ndim != 2:
        raise ValueError(f"Expected 2D pressure_data, got shape {pressure_data.shape}")

    if pressure_data.shape[0] == num_points:
        # shape: (points, time)
        num_time_steps = pressure_data.shape[1]

        def get_frame(t):
            return pressure_data[:, t]

    elif pressure_data.shape[1] == num_points:
        # shape: (time, points)
        num_time_steps = pressure_data.shape[0]

        def get_frame(t):
            return pressure_data[t, :]

    else:
        raise ValueError(
            f"Cannot interpret pressure_data shape {pressure_data.shape}; "
            f"expected one axis to be width*height={num_points}"
        )

    saved_frames = []

    for t in range(0, num_time_steps, frame_step):
        pressure = get_frame(t).reshape((height, width))

        path = os.path.join(output_dir, f"frame_{t:04d}.vti")
        write_vti(
            path=path,
            pressure=pressure,
            sound_speed=sound_speed,
            density=density,
            alpha_coeff=alpha_coeff,
            width=width,
            height=height,
        )

        saved_frames.append((t, f"frame_{t:04d}.vti"))

    write_pvd(os.path.join(output_dir, "frames.pvd"), saved_frames)


def write_vti(path, pressure, sound_speed, density, alpha_coeff, width, height):
    image_data = vtk.vtkImageData()
    image_data.SetDimensions(width, height, 1)
    image_data.SetSpacing(1.0, 1.0, 1.0)
    image_data.SetOrigin(0.0, 0.0, 0.0)

    point_data = image_data.GetPointData()

    def add_array(name, array, set_active=False):
        contiguous = np.ascontiguousarray(array, dtype=np.float32)
        vtk_array = numpy_to_vtk(contiguous.ravel(order="C"), deep=True)
        vtk_array.SetName(name)
        point_data.AddArray(vtk_array)
        if set_active:
            point_data.SetScalars(vtk_array)

    # NumPy arrays are stored as (height, width), which matches VTK's x-fastest
    # ordering when flattened in C order with dimensions set to (width, height, 1).
    add_array("pressure", pressure, set_active=True)
    add_array("sound_speed", sound_speed)
    add_array("density", density)
    add_array("alpha_coeff", alpha_coeff)

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(path)
    writer.SetInputData(image_data)
    writer.SetDataModeToAscii()
    writer.SetCompressorTypeToNone()

    write_result = writer.Write()
    has_nan = bool(np.isnan(pressure).any())
    has_inf = bool(np.isinf(pressure).any())
    finite = np.isfinite(pressure)
    if finite.any():
        pressure_min = float(pressure[finite].min())
        pressure_max = float(pressure[finite].max())
    else:
        pressure_min = float("nan")
        pressure_max = float("nan")

    print(
        f"{path} write={write_result} pressure_min={pressure_min} "
        f"pressure_max={pressure_max} has_nan={has_nan} has_inf={has_inf}"
    )

    if write_result != 1:
        raise IOError(f"Failed to write VTI file: {path}")


def write_pvd(path, frames):
    with open(path, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
        f.write("  <Collection>\n")

        for timestep, filename in frames:
            f.write(f'    <DataSet timestep="{timestep}" group="" part="0" file="{filename}"/>\n')

        f.write("  </Collection>\n")
        f.write("</VTKFile>\n")


def main():
    width = 512
    height = 512

    dx = 1e-4  # 0.1 mm
    dy = dx

    kgrid = kWaveGrid([height, width], [dy, dx])

    # Реальные единицы, не relative.
    c_water = 1500.0
    rho_water = 1000.0

    sound_speed = c_water * np.ones((height, width), dtype=np.float32)
    density = rho_water * np.ones((height, width), dtype=np.float32)
    alpha_coeff = np.zeros((height, width), dtype=np.float32)

    # Bone disk
    bone_mask = make_disc(
        Vector([height, width]),
        Vector([height // 2, width // 2]),
        35,
    )

    sound_speed[bone_mask == 1] = 3000.0
    density[bone_mask == 1] = 1850.0
    alpha_coeff[bone_mask == 1] = 20.0

    medium = kWaveMedium(
        sound_speed=sound_speed,
        density=density,
        alpha_coeff=alpha_coeff,
        alpha_power=1.1,
    )

    # Время симуляции
    c_max = float(np.max(sound_speed))
    cfl = 0.25
    dt = cfl * dx / c_max
    total_steps = 3500
    kgrid.setTime(total_steps, dt)

    # Source
    source = kSource()
    source.p_mask = np.zeros((height, width), dtype=bool)
    source_y = 50
    source_x = width // 2
    source.p_mask[source_y, source_x] = True

    # Аналог Source::new + gaussian-windowed tone burst
    source_freq = 1.0e6
    source_cycles = 5
    source.p = tone_burst(1 / dt, source_freq, source_cycles)

    # Sensors
    sensor = kSensor()
    sensor.mask = np.ones((height, width), dtype=bool)
    sensor.record = ["p"]

    simulation_options = SimulationOptions(
        pml_inside=False,
        pml_size=[80, 80],
        data_cast="single",
        save_to_disk=True,
    )

    sensor_data = kspaceFirstOrder(
        kgrid=kgrid,
        source=source,
        sensor=sensor,
        medium=medium,
        **options_to_kwargs(simulation_options),
        backend="cpp",
        device="gpu",
    )

    if isinstance(sensor_data, dict):
        pressure_data = sensor_data["p"]
    else:
        pressure_data = sensor_data.p

    export_pressure_timeseries_to_vti(
        pressure_data=pressure_data,
        width=width,
        height=height,
        sound_speed=sound_speed,
        density=density,
        alpha_coeff=alpha_coeff,
        output_dir="output/run_001",
        frame_step=10,
    )


if __name__ == "__main__":
    main()
