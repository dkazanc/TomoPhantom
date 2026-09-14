# TomoPhantom NXtomo generator: realistic data and optional artefacts

This is a replacement for `Demos/nxs/nxs_generator.py`, with the original
`-m`, `-s`, `-a`, and `-o` arguments retained. The default shape is deliberately
small: `(128, 256, 362)` in **detector height, projection count, detector width**
order. Install TomoPhantom in your environment; this bundle does not contain it.

```bash
conda install -c httomo -c conda-forge tomophantom h5py psutil
```

Some artefacts, particularly subpixel shifts, also require `scikit-image`.
No GPU, ASTRA or reconstruction package is needed to run this generator.

## Generate data

Realistic background, flats and darks:

```bash
python nxs_generator.py --realistic -m 13 -s 128 256 362 \
  --flats 20 --darks 10 --source-intensity 50000 \
  --flat-settings flat_settings.json --seed 1 -o realistic.nxs
```

Add stripes and zingers to the simulated detector signal:

```bash
python nxs_generator.py --realistic -m 13 -s 128 256 362 \
  --artefacts artefacts.json --artefacts-stage raw \
  --flat-settings flat_settings.json --seed 1 -o artefacts.nxs
```

The example JSON files pass the documented keyword arguments to
`synth_flats` and `artefacts_mix`. Unknown keys are rejected; the implementation's
stripe keyword is `stripes_type`, not the singular spelling in its docstring.
The `detectors_miscallibration` and `fresnel_wavelenght` spellings are upstream API names.

For blur, Fresnel propagation, or line-integral noise experiments, use
`--artefacts-stage line-integrals`. For example, save this as `blur.json`:

```json
{"pve_strength": 1, "verbose": false}
```

```bash
python nxs_generator.py --realistic -s 128 256 362 \
  --artefacts blur.json --artefacts-stage line-integrals -o blurred.nxs
```

To simulate motion, use one positive `datashifts_maxamplitude_pixel` or
`datashifts_maxamplitude_subpixel` setting. Returned shifts are saved under
`/entry/generator/artefact_shifts`, in the API's native component order. The
`projection_frame_indices` dataset maps these to the combined frame stack.

Clean line-integral data without calibration frames:

```bash
python nxs_generator.py -m 13 -s 128 256 362 --chunk-rows 16 -o clean.nxs
```

Use `python nxs_generator.py --help` for all options. Existing output is protected
unless `--overwrite` is supplied. Output is published only after writing completes.

## Data layout and image keys

The frame order is **darks, flats, projections**. All are stored in one uint16
dataset; calibration images are selected by their image keys, not by separate files.

| Content | HDF5 path |
| --- | --- |
| Combined frame data | `/entry/data/data` |
| Image keys | `/entry/instrument/detector/image_key` |
| Frame-aligned rotation angles | `/entry/data/rotation_angle` |
| Detector data (same hard-linked dataset) | `/entry/instrument/detector/data` |
| Image keys (same hard-linked dataset) | `/entry/data/image_key` |
| Sample angles (same hard-linked dataset) | `/entry/sample/rotation_angle` |
| Configuration and provenance | `/entry/generator` |

| Image key | Frame type | Count in the default realistic example |
| --- | --- | --- |
| `2` | Dark field | 10 |
| `1` | Flat field | 20 |
| `0` | Sample projection | 256 |

For that example, stored shape is `(286, 128, 362)`: **frame, detector y, detector x**.
The key and angle arrays each have length 286. Calibration frames have the initial
stage angle; use `image_key == 0` to obtain the 256 projection angles.
An inclusive `linspace` uses the supplied angular endpoints, as in the original script.

## Realistic generation model

1. Generate non-negative analytical projections with `TomoP3D.ModelSinoSub`.
2. Optionally call `artefacts_mix` in the line-integral domain.
3. Call `synth_flats` once for the entire detector volume, returning raw projections,
   simulated flat frames and a detector-response map. The map is not stored.
4. Optionally call `artefacts_mix` on the raw signal, scaled to approximately `[0, 1]`.
5. Apply the same signal gain to projections and flats, add a dark-current/read-noise
   model, then round, clip and convert to uint16. Generate dark frames with no beam signal.

For raw projections/flats, before rounding and clipping:

`stored = signal_scale * synth_flats_signal + Poisson(dark_mean) + Gaussian(0, read_noise)`

Dark frames use the same expression with the beam signal set to zero. Defaults are
`signal_scale=0.9`, `dark_mean=100` and `read_noise=2` standard deviation.
This is a simple additive detector-dark model, not measured detector calibration.
The original realistic example uses `darks=None`; `synth_flats` does not return darks.

Both projections and flats use a common added gain, and no new maximum rescaling is
performed after raw artefacts. Out-of-range values are clipped, not wrapped; the
clipped voxel count is reported and stored. `synth_flats` itself normalises projections
and flats internally, so its output is synthetic detector signal rather than absolutely
calibrated photon counts. `source_intensity` controls its noise model, not the stored ADU maximum.

`artefacts_mix` has its own fixed application order (PVE, Fresnel, zingers, stripes,
shifts, noise), regardless of JSON key order. Artefacts affect projections, not
calibration frames. Gaussian `noise_amplitude` follows upstream behaviour: it is the
standard deviation relative to that function's normalised data.

The Poisson branch in `artefacts_mix` exponentiates its input; it must not receive
already-raw detector counts. Raw-stage Poisson/Fresnel and `noise_prelog=True` are
therefore rejected. Realistic synthesis already adds its own noise, so explicitly
enabling line-integral noise adds another noise contribution.

## Memory, scaling and repeatability

- Clean mode without artefacts remains chunked. It uses a first pass to find a
  **global maximum**, then a second pass to write uint16 data. This removes the original
  per-chunk min/max scaling, which could cause seams and inconsistent attenuation.
- Realistic mode or any artefact configuration requires a **full-volume working array**.
  The chunk-rows setting limits analytic generation and writing, not realistic synthesis.
  Independently synthesising chunks would reset flat-field background, normalisation
  and motion/PVE context at each boundary. This version does not claim arbitrary-size
  out-of-core realistic generation; that requires a stateful streaming synthesis API.
- A heuristic working-memory check uses 50% of available RAM by default. Use
  `--memory-fraction` or `--max-memory-gb` to adjust the budget. It is not a hard process
  memory limit; reduce the dimensions if the check fails.
- Clean data are scaled line integrals, not raw transmission: do not apply dark/flat
  correction and minus-log to clean-mode output. The storage scale is recorded.
- Use the same seed, parameters, chunk-rows and library version for reproducibility.
  Different chunk sizes can change the added dark/read-noise realisation. Upstream
  routines also contain their own RNG behaviour; cross-version bitwise reproducibility
  is not guaranteed.

Sources:
- [Original generator](https://github.com/dkazanc/TomoPhantom/blob/master/Demos/nxs/nxs_generator.py)
- [Realistic example](https://github.com/dkazanc/TomoPhantom/blob/afadc7c522919e08cdae5aeef6db41cded0fc213/Demos/3D/ReconASTRA3D_realistic.py)
- [synth_flats](https://github.com/dkazanc/TomoPhantom/blob/afadc7c522919e08cdae5aeef6db41cded0fc213/tomophantom/flatsgen.py)
- [artefacts_mix](https://github.com/dkazanc/TomoPhantom/blob/afadc7c522919e08cdae5aeef6db41cded0fc213/tomophantom/artefacts.py)
- [NXtomo application definition](https://manual.nexusformat.org/classes/applications/NXtomo.html)
