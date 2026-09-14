#!/usr/bin/env python3
"""Generate TomoPhantom projections with optional realistic detector data/artefacts.

Adapted from TomoPhantom's Demos/nxs/nxs_generator.py (Apache-2.0).
Axes inside TomoPhantom: (detector_y, angle, detector_x).
Axes on disk: (frame, detector_y, detector_x).

Clean generation is streamed in detector-row chunks with ONE global scale.
Realistic synthesis and artefacts require a full volume: synth_flats has global
normalisation/background state, and shifts/PVE need neighbouring detector rows.
The memory estimate is a conservative planning heuristic, not an OOM guarantee.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.metadata
import json
import math
from pathlib import Path
import random
import tempfile

import h5py
import numpy as np
import psutil

ARTEFACT_KEYS = {
    "noise_type",
    "noise_amplitude",
    "noise_seed",
    "noise_prelog",
    "zingers_percentage",
    "zingers_modulus",
    "stripes_percentage",
    "stripes_maxthickness",
    "stripes_intensity",
    "stripes_type",
    "stripes_variability",
    "datashifts_maxamplitude_pixel",
    "datashifts_maxamplitude_subpixel",
    "pve_strength",
    "fresnel_dist_observation",
    "fresnel_scale_factor",
    "fresnel_wavelenght",
    "verbose",
}


def parser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("-m", "--model-number", type=int, default=13)
    p.add_argument(
        "-s",
        "--sinogram-shape",
        type=int,
        nargs=3,
        metavar=("HEIGHT", "ANGLES", "WIDTH"),
        default=(128, 256, 362),
    )
    p.add_argument("-a", "--angle-range", type=float, nargs=2, default=(0, 179.9))
    p.add_argument("-o", "--output-path", type=Path, default=Path("phantom.nxs"))
    p.add_argument(
        "--realistic", action="store_true", help="Use synth_flats; add darks."
    )
    p.add_argument("--flats", type=int, default=20)
    p.add_argument("--darks", type=int, default=10)
    p.add_argument("--source-intensity", type=int, default=50000)
    p.add_argument("--flat-settings", type=Path, help="JSON overrides for synth_flats.")
    p.add_argument(
        "--dark-mean", type=float, default=100.0, help="Dark current in ADU."
    )
    p.add_argument("--read-noise", type=float, default=2.0, help="Gaussian SD in ADU.")
    p.add_argument(
        "--signal-scale",
        type=float,
        default=0.9,
        help="Common gain on synth_flats uint16 signals; reserves headroom.",
    )
    p.add_argument(
        "--artefacts", type=Path, help="JSON keyword arguments to artefacts_mix."
    )
    p.add_argument(
        "--artefacts-stage",
        choices=["raw", "line-integrals"],
        default="raw",
        help="raw: before dark offset, on normalised detector counts.",
    )
    p.add_argument("--seed", type=int, default=1)
    p.add_argument(
        "--chunk-rows",
        type=int,
        default=32,
        help="Rows per analytic-generation/write chunk (not synth_flats).",
    )
    p.add_argument("--memory-fraction", type=float, default=0.5)
    p.add_argument(
        "--max-memory-gb", type=float, help="Additional working-memory cap (GiB)."
    )
    p.add_argument("--compression", choices=["none", "gzip", "lzf"], default="none")
    p.add_argument("--overwrite", action="store_true")
    return p


def read_json(path):
    if path is None:
        return {}
    obj = json.loads(path.read_text())
    if not isinstance(obj, dict):
        raise ValueError(f"{path}: expected a JSON object")
    return obj


def settings(args):
    import inspect
    from tomophantom.flatsgen import synth_flats

    if min(args.sinogram_shape) < 2:
        raise ValueError("Each shape dimension must be at least 2.")
    if args.chunk_rows < 1 or not 0 < args.memory_fraction <= 1:
        raise ValueError(
            "chunk-rows must be positive; memory-fraction must be in (0,1]."
        )
    if (
        not np.isfinite(args.angle_range).all()
        or args.angle_range[1] <= args.angle_range[0]
    ):
        raise ValueError("The angle range must be finite and increasing.")
    if args.max_memory_gb is not None and (
        not math.isfinite(args.max_memory_gb) or args.max_memory_gb <= 0
    ):
        raise ValueError("max-memory-gb must be finite and positive.")
    if not 0 <= args.seed < 2**32:
        raise ValueError("seed must be between 0 and 2**32-1.")
    if not np.isfinite([args.dark_mean, args.read_noise, args.signal_scale]).all():
        raise ValueError("Detector settings must be finite.")
    if args.dark_mean < 0 or args.read_noise < 0 or not 0 < args.signal_scale <= 1:
        raise ValueError("Invalid dark mean, read noise or signal scale.")
    if args.realistic and (min(args.flats, args.darks, args.source_intensity) < 1):
        raise ValueError(
            "Realistic generation requires positive flat/dark counts and intensity."
        )
    flat = read_json(args.flat_settings)
    allowed = set(inspect.signature(synth_flats).parameters) - {
        "projData3D_clean",
        "source_intensity",
        "flatsnum",
    }
    if set(flat) - allowed:
        raise ValueError(f"Unsupported flat settings: {sorted(set(flat) - allowed)}")
    if flat and not args.realistic:
        raise ValueError("--flat-settings requires --realistic.")
    flat = dict(source_intensity=args.source_intensity, flatsnum=args.flats, **flat)
    if not 1 <= flat.get("variations_number", 3) <= 3:
        raise ValueError("synth_flats supports 1 to 3 detector-response variations.")
    art = read_json(args.artefacts)
    if set(art) - ARTEFACT_KEYS:
        raise ValueError(f"Unknown artefact keys: {sorted(set(art) - ARTEFACT_KEYS)}")
    if art and not args.realistic and args.artefacts_stage == "raw":
        raise ValueError("Use --realistic or --artefacts-stage line-integrals.")
    if art.get("noise_prelog"):
        raise ValueError(
            "noise_prelog is unsupported; --realistic provides raw data and flats."
        )
    if art.get("noise_type") not in (None, "Gaussian", "Poisson"):
        raise ValueError("noise_type must be Gaussian or Poisson.")
    if args.artefacts_stage == "raw" and (
        art.get("noise_type") == "Poisson"
        or art.get("fresnel_dist_observation") is not None
    ):
        raise ValueError(
            "Poisson/Fresnel in artefacts_mix expects line-integral input. "
            "Use --artefacts-stage line-integrals. synth_flats already adds noise."
        )
    shifts = [k for k in art if k.startswith("datashifts_") and art[k] is not None]
    if len(shifts) > 1 or any(art[k] <= 0 for k in shifts):
        raise ValueError(
            "Specify at most one shift type, with strictly positive amplitude."
        )
    return flat, art


def checked(data, shape, label):
    if data.shape != shape or not np.isfinite(data).all():
        raise ValueError(
            f"{label}: invalid shape or non-finite data; expected {shape}."
        )
    return data


def apply_artefacts(data, config, seed):
    """Use the API's own order; preserve shifts when it returns [data, shifts]."""
    from tomophantom.artefacts import artefacts_mix

    if not config:
        return data, None
    config = dict(config)
    noise_seed = config.pop("noise_seed", seed)
    # The cited API checks `seed is True`, not an arbitrary integer seed.
    random.seed(seed)
    np.random.seed(seed if noise_seed is None else int(noise_seed))
    config.update(noise_seed=False, noise_prelog=False)
    config.setdefault("verbose", False)
    result = artefacts_mix(np.array(data, dtype=np.float32, copy=True), **config)
    shifts = None
    if isinstance(result, (tuple, list)):
        result, shifts = result
    checked(result, data.shape, "artefacts_mix")
    return result, shifts


def nxgroup(parent, name, nxclass):
    group = parent.create_group(name)
    group.attrs["NX_class"] = nxclass
    return group


def initialise(file, args, angles, metadata):
    height, count, width = args.sinogram_shape
    nd, nf = (args.darks, args.flats) if args.realistic else (0, 0)
    keys = np.concatenate([np.full(nd, 2), np.full(nf, 1), np.zeros(count)]).astype(
        "int8"
    )
    # Calibration frames are taken at the initial stage angle and identified by key.
    frame_angles = np.concatenate([np.full(nd + nf, angles[0]), angles]).astype(
        "float32"
    )
    file.attrs.update(
        default="entry",
        file_name=args.output_path.name,
        file_time=datetime.now(timezone.utc).isoformat(),
    )
    entry = nxgroup(file, "entry", "NXentry")
    entry.attrs["default"] = "data"
    entry.create_dataset("definition", data="NXtomo")
    entry.create_dataset("title", data=f"TomoPhantom model {args.model_number}")
    entry.create_dataset("start_time", data=datetime.now(timezone.utc).isoformat())
    instrument = nxgroup(entry, "instrument", "NXinstrument")
    detector = nxgroup(instrument, "detector", "NXdetector")
    data = detector.create_dataset(
        "data",
        shape=(len(keys), height, width),
        dtype="uint16",
        chunks=(1, min(height, args.chunk_rows), min(width, 512)),
        compression=None if args.compression == "none" else args.compression,
    )
    data.attrs.update(
        units="counts" if args.realistic else "arbitrary",
        interpretation="image",
        axis_order="frame,detector_y,detector_x",
    )
    key = detector.create_dataset("image_key", data=keys)
    key.attrs["meaning"] = "0=projection; 1=flat; 2=dark"
    sample = nxgroup(entry, "sample", "NXsample")
    sample.create_dataset("name", data=f"Synthetic model {args.model_number}")
    rotation = sample.create_dataset("rotation_angle", data=frame_angles)
    rotation.attrs["units"] = "deg"
    nxdata = nxgroup(entry, "data", "NXdata")
    nxdata.attrs["signal"] = "data"
    nxdata.attrs["axes"] = np.asarray(
        ["rotation_angle", ".", "."], dtype=h5py.string_dtype()
    )
    nxdata["data"], nxdata["image_key"], nxdata["rotation_angle"] = data, key, rotation
    process = nxgroup(entry, "generator", "NXprocess")
    process.create_dataset("program", data="nxs_generator.py")
    process.create_dataset("version", data="2.0")
    process.create_dataset("date", data=datetime.now(timezone.utc).isoformat())
    process.create_dataset("sequence_index", data=1)
    process.create_dataset(
        "configuration", data=json.dumps(metadata, default=str, sort_keys=True)
    )
    process.create_dataset(
        "projection_frame_indices", data=np.arange(nd + nf, len(keys))
    )
    return data, process, nd, nf


def main(args):
    import tomophantom
    from tomophantom import TomoP3D
    from tomophantom.flatsgen import synth_flats

    flat_config, art = settings(args)
    args.output_path = Path(args.output_path)
    if args.output_path.exists() and not args.overwrite:
        raise FileExistsError(
            f"{args.output_path} exists; use --overwrite to replace it."
        )
    height, count, width = args.sinogram_shape
    rows = min(height, args.chunk_rows)
    angles = np.linspace(*args.angle_range, count, dtype="float32")
    library = str(
        Path(tomophantom.__file__).parent / "phantomlib" / "Phantom3DLibrary.dat"
    )
    full = args.realistic or bool(art)
    budget = int(psutil.virtual_memory().available * args.memory_fraction)
    if args.max_memory_gb is not None:
        budget = min(budget, int(args.max_memory_gb * 2**30))
    # Includes float work arrays, calibration stacks, maps and projection scratch.
    estimate = (
        height * count * width * 64 + height * width * (8 * args.flats + 256)
        if full
        else rows * count * width * 16
    ) + 64 * 2**20
    if estimate > budget:
        raise MemoryError(
            f"Estimated working memory {estimate / 2**30:.2f} GiB exceeds "
            f"budget {budget / 2**30:.2f} GiB. Reduce shape"
            + (
                "; realistic/artefact modes require a full volume."
                if full
                else " or --chunk-rows."
            )
        )
    print(
        f"Estimated working memory: {estimate / 2**30:.2f} GiB; "
        f'mode: {"full-volume synthesis" if full else "streamed clean projections"}'
    )
    random.seed(args.seed)
    np.random.seed(args.seed)
    rng = np.random.default_rng(args.seed)

    def analytic(start, stop):
        data = TomoP3D.ModelSinoSub(
            args.model_number,
            int(width / np.sqrt(2)),
            width,
            height,
            (start, stop),
            angles,
            library,
        )
        checked(data, (stop - start, count, width), "ModelSinoSub")
        if data.min() < 0:
            raise ValueError(
                "This generator requires non-negative attenuation projections."
            )
        return data

    projections = None
    if full:
        projections = np.empty((height, count, width), dtype="float32")
        for start in range(0, height, rows):
            stop = min(start + rows, height)
            projections[start:stop] = analytic(start, stop)
        clean_max = float(projections.max())
    else:
        # First pass gives one common scaling for all chunks, avoiding seams.
        clean_max = max(
            float(analytic(s, min(s + rows, height)).max())
            for s in range(0, height, rows)
        )
    if clean_max <= 0:
        raise ValueError("The selected phantom has no positive projection data.")
    shifts = None
    if art and args.artefacts_stage == "line-integrals":
        projections, shifts = apply_artefacts(projections, art, args.seed)
        np.maximum(projections, 0, out=projections)
    flats = None
    synthesis_peak = None
    if args.realistic:
        synthesis_peak = float(projections.max())
        if synthesis_peak <= 0:
            raise ValueError(
                "Artefacts removed all positive signal before synth_flats."
            )
        projections, flats, response_map = synth_flats(projections, **flat_config)
        checked(projections, (height, count, width), "synth_flats projections")
        checked(flats, (height, args.flats, width), "synth_flats flats")
        if not np.isfinite(response_map).all():
            raise ValueError("synth_flats produced an invalid detector-response map.")
        del response_map
        if art and args.artefacts_stage == "raw":
            projections, shifts = apply_artefacts(
                projections.astype("float32") / 65535, art, args.seed
            )
            projections *= 65535  # No new max normalisation after corruption.
    clean_scale = float(projections.max()) if full and not args.realistic else clean_max
    clean_scale = max(clean_scale, np.finfo("float32").tiny)
    try:
        version = importlib.metadata.version("tomophantom")
    except importlib.metadata.PackageNotFoundError:
        version = getattr(tomophantom, "__version__", "source checkout / unknown")
    metadata = dict(
        arguments=vars(args),
        flat_settings=flat_config if args.realistic else {},
        artefacts=art,
        tomophantom_version=version,
        clean_projection_max=clean_max,
        clean_storage_scale=None if args.realistic else clean_scale,
        synth_flats_input_max=synthesis_peak,
        data_domain=(
            "raw_detector_counts" if args.realistic else "scaled_line_integrals"
        ),
        frame_order="darks, flats, projections",
        estimated_working_bytes=estimate,
    )
    args.output_path.parent.mkdir(parents=True, exist_ok=True)
    # Write beside destination, then publish only a complete file.
    temp = tempfile.NamedTemporaryFile(
        dir=args.output_path.parent,
        prefix=args.output_path.name + ".",
        suffix=".partial",
        delete=False,
    )
    temporary = Path(temp.name)
    temp.close()
    clipped = 0

    def quantise(block, detector_signal=False):
        nonlocal clipped
        values = np.array(block, dtype="float32", copy=True)
        if detector_signal:
            values *= args.signal_scale
            values += rng.poisson(args.dark_mean, size=values.shape).astype("float32")
            if args.read_noise:
                values += rng.normal(0, args.read_noise, size=values.shape).astype(
                    "float32"
                )
        if not np.isfinite(values).all():
            raise ValueError("Non-finite values encountered before uint16 conversion.")
        clipped += int(np.count_nonzero((values < 0) | (values > 65535)))
        return np.rint(np.clip(values, 0, 65535)).astype("uint16")

    try:
        with h5py.File(temporary, "w") as f:
            output, process, nd, nf = initialise(f, args, angles, metadata)
            # One shared dark-current distribution for projections, flats and darks.
            for frame in range(nd):
                output[frame] = quantise(np.zeros((height, width), "float32"), True)
            for frame in range(nf):
                output[nd + frame] = quantise(flats[:, frame, :], True)
            for start in range(0, height, rows):
                stop = min(start + rows, height)
                block = projections[start:stop] if full else analytic(start, stop)
                if not args.realistic:
                    block = block / clean_scale * 65535
                output[nd + nf :, start:stop] = np.swapaxes(
                    quantise(block, args.realistic), 0, 1
                )
            if shifts is not None:
                ds = process.create_dataset("artefact_shifts", data=shifts)
                ds.attrs.update(
                    units="pixel",
                    description="As returned by TomoPhantom; projection frames only",
                )
            process.create_dataset("clipped_voxels", data=clipped)
            f["entry"].create_dataset(
                "end_time", data=datetime.now(timezone.utc).isoformat()
            )
        if args.overwrite:
            temporary.replace(args.output_path)
        else:
            # Atomic no-clobber publication on the same filesystem.
            import os

            os.link(temporary, args.output_path)
            temporary.unlink()
    finally:
        temporary.unlink(missing_ok=True)
    print(
        f"Wrote {args.output_path}: {nd} darks + {nf} flats + {count} projections; "
        f"frame shape {height} x {width}; clipped voxels: {clipped}"
    )


if __name__ == "__main__":
    p = parser()
    try:
        main(p.parse_args())
    except (ValueError, MemoryError, FileExistsError, OSError) as exc:
        p.exit(2, f"Error: {exc}\n")
