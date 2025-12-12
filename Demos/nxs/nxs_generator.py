#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import pathlib
import psutil
import numpy as np
import tomophantom

import h5py
from httomo.methods import save_intermediate_data
from httomo.runner.auxiliary_data import AuxiliaryData
from httomo.runner.dataset import DataSetBlock


def save_nxs(
    file: h5py.File,
    aux_data: AuxiliaryData,
    data,
    chunk_start,
    global_shape,
    detector_x,
    detector_y,
):
    b1 = DataSetBlock(
        data=data,
        aux_data=aux_data,
        slicing_dim=1,
        block_start=0,
        chunk_start=chunk_start,
        chunk_shape=data.shape,
        global_shape=global_shape,
    )

    save_intermediate_data(
        b1.data,
        b1.global_shape,
        b1.global_index,
        b1.slicing_dim,
        file,
        frames_per_chunk=0,
        minimum_block_length=global_shape[0],
        path="/data",
        detector_x=detector_x,
        detector_y=detector_y,
        angles=b1.angles,
    )


def main(args: argparse.Namespace):
    model = args.model_number
    path = pathlib.Path(tomophantom.__file__).parent
    path_library3D = str(path / "phantomlib" / "Phantom3DLibrary.dat")

    Horiz_det = args.sinogram_shape[2]
    Vert_det = args.sinogram_shape[0]
    angles_num = args.sinogram_shape[1]
    global_shape = (angles_num, Vert_det, Horiz_det)
    square_phantom_slice_width = int(Horiz_det / np.sqrt(2))

    angles = np.linspace(*args.angle_range, angles_num, dtype="float32")  # in degrees
    aux_data = AuxiliaryData(angles=angles)

    available_memory_bytes = psutil.virtual_memory().available
    slice_memory_bytes = Horiz_det * angles_num * np.float32().itemsize
    if slice_memory_bytes > available_memory_bytes:
        print(f"A signle slice ({slice_memory_bytes} bytes) is would be bigger than available memory ({available_memory_bytes} bytes).")
        return

    chunk_size = available_memory_bytes // slice_memory_bytes
    chunk_count = int(np.ceil(Vert_det / chunk_size))
    print(f"Creating phantom in {chunk_count} number of chunks of at max {chunk_size} slices per chunk.")

    with h5py.File(args.output_path, "w") as file:
        file.attrs["default"] = "entry"
        entry = file.create_group("entry")
        entry.attrs["NX_class"] = "NXentry"
        entry.attrs["default"] = "data"
        entry.attrs["NX_application"] = "NXtomo"

        entry.create_dataset("definition", data="NXtomo")
        data_group = entry.create_group("data")
        data_group["data"] = h5py.SoftLink("/data")
        data_group["rotation_angle"] = h5py.SoftLink("/angles")

        instrument = entry.create_group("instrument")
        detector = instrument.create_group("detector")
        image_key_data = np.zeros_like(angles, dtype=np.int8)
        detector.create_dataset("image_key", data=image_key_data, dtype=np.int8)

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

        with h5py.File(args.output_path, "a") as file:
            save_nxs(
                file,
                aux_data,
                swapped_projData3D_analyt,
                chunk_start,
                global_shape,
                Horiz_det,
                Vert_det,
            )

    print(
        f"#slices: {Vert_det}, #projections: {angles_num}, detector width: {Horiz_det}"
    )
    print(f"phantom shape in output file: {global_shape}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TomoPhantom phantom .nxs generator. Output can be used with httomo, e.g. 'httomo run phantom.nxs pipeline.yaml output'.")
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
