import os.path
from pathlib import Path

import h5py
import numpy as np

from kwave.data import Vector
from kwave.kgrid import kWaveGrid
from kwave.kmedium import kWaveMedium
from kwave.ksource import kSource
from kwave.ksensor import kSensor
from kwave.kspaceFirstOrder import kspaceFirstOrder
from kwave.options.simulation_options import SimulationOptions
from kwave.utils.mapgen import make_disc
from kwave.utils.signals import tone_burst
from kwave.compat import options_to_kwargs

def make_point_source(
        height: float,
        width: float,
        source_x: float,
        source_y: float,
        dt: float,
        source_freq: float,
        source_cycles: int,
):
    source = kSource()
    source.p_mask = np.zeros((height, width), dtype=bool)
    source.p_mask[source_y, source_x] = True
    source.p = tone_burst(1 / dt, source_freq, source_cycles)
    return source


def make_electronically_focused_linear_aperture_source(
        height,
        width,
        dx,
        dt,
        total_steps,
        sound_speed,
        aperture_y,
        aperture_width,
        focus_position,
        source_freq,
        source_cycles,
        source_pressure,
):
    """Build a finite-aperture focused pressure source for a 2D transducer."""
    source = kSource()
    source.p_mask = np.zeros((height, width), dtype=bool)

    aperture_points = max(3, int(round(aperture_width / dx)))
    if aperture_points % 2 == 0:
        aperture_points += 1

    aperture_center_x = width // 2
    half_aperture = aperture_points // 2
    source_x = np.arange(
        aperture_center_x - half_aperture,
        aperture_center_x + half_aperture + 1,
    )
    source_x = source_x[(source_x >= 0) & (source_x < width)]
    source_y = np.full_like(source_x, aperture_y)

    source.p_mask[source_y, source_x] = True

    base_signal = tone_burst(1 / dt, source_freq, source_cycles).reshape(-1)
    focus_y, focus_x = focus_position
    element_y_m = source_y * dx
    element_x_m = source_x * dx
    focus_y_m = focus_y * dx
    focus_x_m = focus_x * dx
    distances = np.sqrt((element_y_m - focus_y_m) ** 2 + (element_x_m - focus_x_m) ** 2)

    # Emit edge elements first and center elements last so all wavefronts meet at the focus.
    delays = (np.max(distances) - distances) / sound_speed
    sample_offsets = np.rint(delays / dt).astype(int)

    apodization = 0.1 + 0.9 * np.hanning(len(source_x))
    if np.max(apodization) > 0:
        apodization /= np.max(apodization)

    source.p = np.zeros((len(source_x), total_steps), dtype=np.float32)
    for element_index, (offset, weight) in enumerate(zip(sample_offsets, apodization)):
        end = min(total_steps, offset + base_signal.size)
        source.p[element_index, offset:end] = (
                source_pressure * weight * base_signal[:end - offset]
        )

    source.p_mode = "additive"
    return source


def write_kwave_xdmf(
        input_h5="./kwave/kwave_output.h5",
        output_xdmf="./kwave/kwave_output.xdmf",
        width=None,
        height=None,
        dx=1e-4,
        dy=1e-4,
        frame_step=10,
):
    input_h5 = Path(input_h5)
    output_xdmf = Path(output_xdmf)

    with h5py.File(input_h5, "r") as f:
        p = f["p"]
        dt = float(f["dt"][0, 0, 0])
        nt = int(p.shape[1])
        num_points = int(p.shape[2])

    if width is None or height is None:
        side = int(num_points ** 0.5)
        if side * side != num_points:
            raise ValueError(
                f"Cannot infer square grid from {num_points} points; "
                f"pass width/height explicitly"
            )
        width = side
        height = side

    if width * height != num_points:
        raise ValueError(
            f"width*height={width*height}, but /p has {num_points} points"
        )

    h5_ref = input_h5.name
    frame_indices = list(range(0, nt, frame_step))

    with open(output_xdmf, "w") as x:
        x.write('<?xml version="1.0" ?>\n')
        x.write('<Xdmf Version="3.0">\n')
        x.write('  <Domain>\n')
        x.write('    <Grid Name="PressureTimeSeries" GridType="Collection" CollectionType="Temporal">\n')

        for t in frame_indices:
            time_value = t * dt

            x.write(f'      <Grid Name="Frame_{t}" GridType="Uniform">\n')
            x.write(f'        <Time Value="{time_value}"/>\n')
            x.write(f'        <Topology TopologyType="2DCoRectMesh" Dimensions="{height} {width}"/>\n')
            x.write('        <Geometry GeometryType="ORIGIN_DXDY">\n')
            x.write('          <DataItem Dimensions="2" Format="XML">0 0</DataItem>\n')
            x.write(f'          <DataItem Dimensions="2" Format="XML">{dy} {dx}</DataItem>\n')
            x.write('        </Geometry>\n')
            x.write('        <Attribute Name="pressure" AttributeType="Scalar" Center="Node">\n')
            x.write(f'          <DataItem ItemType="HyperSlab" Dimensions="{height} {width}" Type="HyperSlab">\n')
            x.write('            <DataItem Dimensions="3 3" Format="XML">\n')
            x.write(f'              0 {t} 0\n')
            x.write('              1 1 1\n')
            x.write(f'              1 1 {num_points}\n')
            x.write('            </DataItem>\n')
            x.write(
                f'            <DataItem Dimensions="1 {nt} {num_points}" '
                f'NumberType="Float" Precision="4" Format="HDF">'
                f'{h5_ref}:/p</DataItem>\n'
            )
            x.write('          </DataItem>\n')
            x.write('        </Attribute>\n')
            x.write('      </Grid>\n')

        x.write('    </Grid>\n')
        x.write('  </Domain>\n')
        x.write('</Xdmf>\n')

    print(f"Wrote {output_xdmf}")


def main(run_name: str):
    save_path = os.path.join("runs", run_name)
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

    # Finite-aperture ultrasonic transducer, focused at the bone center.
    source_freq = 1.0e6
    source_cycles = 5
    source = make_electronically_focused_linear_aperture_source(
        height=height,
        width=width,
        dx=dx,
        dt=dt,
        total_steps=total_steps,
        sound_speed=c_water,
        aperture_y=50,
        aperture_width=16e-3,
        focus_position=(height // 2, width // 2),
        source_freq=source_freq,
        source_cycles=source_cycles,
        source_pressure=1.0e5,
    )

    # Sensors
    sensor = kSensor()
    sensor.mask = np.ones((height, width), dtype=bool)
    sensor.record = ["p"]

    simulation_options = SimulationOptions(
        pml_inside=False,
        pml_size=[80, 80],
        data_cast="single",
        save_to_disk=True,
        data_path=save_path,
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

    write_kwave_xdmf(
        input_h5=os.path.join(save_path, "kwave_output.h5"),
        output_xdmf=os.path.join(save_path, "kwave_output.xdmf"),
        frame_step=10,
    )


if __name__ == "__main__":
    main("003")
