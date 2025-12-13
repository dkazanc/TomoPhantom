#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import h5py
import numpy as np
import pathlib
import psutil
import tomophantom


def main(args: argparse.Namespace):
    model = args.model_number
    path = pathlib.Path(tomophantom.__file__).parent
    path_library3D = str(path / "phantomlib" / "Phantom3DLibrary.dat")

    sinogram_dtype = np.float32()
    Horiz_det = args.sinogram_shape[2]
    Vert_det = args.sinogram_shape[0]
    angles_num = args.sinogram_shape[1]
    square_phantom_slice_width = int(Horiz_det / np.sqrt(2))

    angles = np.linspace(*args.angle_range, angles_num, dtype="float32")  # in degrees

    available_memory_bytes = psutil.virtual_memory().available
    slice_memory_bytes = Horiz_det * angles_num * sinogram_dtype.itemsize
    if slice_memory_bytes > available_memory_bytes:
        print(
            f"A signle slice ({slice_memory_bytes} bytes) is would be bigger than available memory ({available_memory_bytes} bytes)."
        )
        return

    chunk_size = available_memory_bytes // slice_memory_bytes
    chunk_count = int(np.ceil(Vert_det / chunk_size))
    print(
        f"Creating phantom in {chunk_count} number of chunks of at max {chunk_size} slices per chunk."
    )

    with h5py.File(args.output_path, "w") as file:
        file.attrs["default"] = "entry"
        entry_group = file.create_group("entry")
        entry_group.attrs["NX_class"] = "NXentry"
        entry_group.attrs["default"] = "data"
        entry_group.attrs["NX_application"] = "NXtomo"

        instrument_group = entry_group.create_group("instrument")
        detector_group = instrument_group.create_group("detector")
        detector_group.create_dataset(
            "image_key", data=np.zeros_like(angles, dtype=np.int8)
        )

        entry_group.create_dataset("definition", data="NXtomo")
        data_group = entry_group.create_group("data")
        sinogram_dataset = data_group.create_dataset(
            "data", (angles_num, Vert_det, Horiz_det), sinogram_dtype
        )
        data_group.create_dataset("rotation_angle", data=angles)

        file.create_dataset(args.output_path.name, data=[0, 0])

        data_dims_group = file.create_group("data_dims")
        data_dims_group.create_dataset("detector_x_y", data=[Horiz_det, Vert_det])

        for i in range(chunk_count):
            chunk_start = i * chunk_size
            chunk_end = min((i + 1) * chunk_size, Vert_det)

            projData3D_analyt = tomophantom.TomoP3D.ModelSinoSub(
                model,
                square_phantom_slice_width,
                Horiz_det,
                Vert_det,
                (chunk_start, chunk_end),
                angles,
                path_library3D,
            )

            swapped_projData3D_analyt = np.swapaxes(projData3D_analyt, 0, 1)
            sinogram_dataset[:angles_num, chunk_start:chunk_end, :Horiz_det] = (
                swapped_projData3D_analyt
            )

    print(
        f"#slices: {Vert_det}, #projections: {angles_num}, detector width: {Horiz_det}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="TomoPhantom phantom .nxs generator. Output can be used with httomo, e.g. 'httomo run phantom.nxs pipeline.yaml output'."
    )
    parser.add_argument(
        "-m",
        "--model-number",
        type=int,
        default=13,  # Shepp-Logan
        help="Model number in TomoPhantom's model library. See https://github.com/dkazanc/TomoPhantom/tree/master/tomophantom/phantomlib.",
    )
    parser.add_argument(
        "-s",
        "--sinogram-shape",
        nargs=3,
        metavar=("detector height", "projection count", "detector width"),
        type=float,
        default=(2368, 256, 4416),
    )
    parser.add_argument(
        "-a",
        "--angle-range",
        nargs=2,
        metavar=("start", "end"),
        type=float,
        default=(0.0, 179.9),
    )
    parser.add_argument(
        "-o", "--output-path", type=pathlib.Path, default="./phantom.nxs"
    )

    args = parser.parse_args()
    main(args)
