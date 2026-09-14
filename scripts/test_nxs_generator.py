"""Run with pytest in an environment with the real TomoPhantom library installed."""

import json

import h5py
import numpy as np
import pytest

import nxs_generator as gen


def args_for(tmp_path, *extra):
    return gen.parser().parse_args(
        [
            "-s",
            "64",
            "24",
            "64",
            "--flats",
            "4",
            "--darks",
            "3",
            "-o",
            str(tmp_path / "phantom.nxs"),
            *extra,
        ]
    )


@pytest.mark.parametrize(
    "stage,config",
    [
        ("raw", {}),
        (
            "raw",
            {
                "zingers_percentage": 0.02,
                "stripes_percentage": 2,
                "stripes_maxthickness": 1,
                "stripes_intensity": 0.05,
            },
        ),
        (
            "line-integrals",
            {"pve_strength": 1, "noise_type": "Gaussian", "noise_amplitude": 0.01},
        ),
        ("raw", {"datashifts_maxamplitude_pixel": 1}),
    ],
)
def test_realistic_frames_and_links(tmp_path, stage, config):
    config_file = tmp_path / "artefacts.json"
    config_file.write_text(json.dumps(config))
    a = args_for(
        tmp_path,
        "--realistic",
        "--artefacts",
        str(config_file),
        "--artefacts-stage",
        stage,
    )
    gen.main(a)
    with h5py.File(a.output_path) as f:
        detector = f["entry/instrument/detector"]
        data = f["entry/data/data"]
        assert data.shape == (31, 64, 64)
        assert data.dtype == np.dtype("uint16")
        assert data.id == detector["data"].id
        assert f["entry/data/image_key"].id == detector["image_key"].id
        assert f["entry/data/rotation_angle"].id == f["entry/sample/rotation_angle"].id
        np.testing.assert_array_equal(
            detector["image_key"][:], [2] * 3 + [1] * 4 + [0] * 24
        )
        angles = f["entry/sample/rotation_angle"][:]
        np.testing.assert_allclose(
            angles[7:], np.linspace(0, 179.9, 24, dtype="float32")
        )
        assert len(angles) == len(detector["image_key"])
        assert f["entry/sample/rotation_angle"].attrs["units"] == "deg"
        assert data[:3].mean() < data[3:7].mean()
        assert data[3:7].mean() > data[7:].mean()
        assert data.maxshape == data.shape
        if "datashifts_maxamplitude_pixel" in config:
            assert f["entry/generator/artefact_shifts"].shape == (24, 2)
            np.testing.assert_array_equal(
                f["entry/generator/projection_frame_indices"][:], np.arange(7, 31)
            )


def test_clean_scaling_is_independent_of_chunk_rows(tmp_path):
    a = args_for(tmp_path, "--chunk-rows", "7")
    gen.main(a)
    with h5py.File(a.output_path) as f:
        expected = f["entry/data/data"][:]
        assert np.all(f["entry/data/image_key"][:] == 0)
    a.chunk_rows = 64
    a.output_path = tmp_path / "whole.nxs"
    gen.main(a)
    with h5py.File(a.output_path) as f:
        np.testing.assert_array_equal(expected, f["entry/data/data"][:])


def test_same_seed_repeats(tmp_path):
    a = args_for(tmp_path, "--realistic")
    gen.main(a)
    with h5py.File(a.output_path) as f:
        first = f["entry/data/data"][:]
    a.output_path = tmp_path / "repeat.nxs"
    gen.main(a)
    with h5py.File(a.output_path) as f:
        np.testing.assert_array_equal(first, f["entry/data/data"][:])


def test_memory_guard_leaves_no_file(tmp_path):
    a = args_for(tmp_path, "--realistic", "--max-memory-gb", ".001")
    with pytest.raises(MemoryError):
        gen.main(a)
    assert not a.output_path.exists()


def test_existing_file_is_preserved(tmp_path):
    a = args_for(tmp_path)
    a.output_path.write_bytes(b"existing file")
    with pytest.raises(FileExistsError):
        gen.main(a)
    assert a.output_path.read_bytes() == b"existing file"


@pytest.mark.parametrize(
    "config",
    [
        {"stripe_type": "full"},  # Actual API keyword is stripes_type.
        {
            "noise_type": "Poisson"
        },  # Cannot use line-integral Poisson API on raw counts.
        {"noise_prelog": True},
        {"datashifts_maxamplitude_pixel": 0},
    ],
)
def test_invalid_artefact_config(tmp_path, config):
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(config))
    with pytest.raises(ValueError):
        gen.main(args_for(tmp_path, "--realistic", "--artefacts", str(path)))
    assert not (tmp_path / "phantom.nxs").exists()


def test_failure_does_not_replace_existing_output(tmp_path, monkeypatch):
    a = args_for(tmp_path, "--overwrite")
    a.output_path.write_bytes(b"keep this")

    def fail(*args):
        raise RuntimeError("simulated writer failure")

    monkeypatch.setattr(gen, "initialise", fail)
    with pytest.raises(RuntimeError):
        gen.main(a)
    assert a.output_path.read_bytes() == b"keep this"
    assert not list(tmp_path.glob("*.partial"))
