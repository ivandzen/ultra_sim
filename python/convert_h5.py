import h5py
from pathlib import Path


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


write_kwave_xdmf(
    input_h5="./kwave/kwave_output.h5",
    output_xdmf="./kwave/kwave_output.xdmf",
    frame_step=10,
)