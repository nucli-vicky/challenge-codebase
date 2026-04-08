# PET Reconstruction

You do **not** need to understand reconstruction or attenuation correction to participate, however, having an intuition of the first reconstruction steps can be a **significant advantage** when designing your pseudo-CT algorithm and loss function. 

## Introduction
In PET, the radioactive tracer (FDG, a glucose analogue) emits positrons that immediately annihilate with nearby electrons, releasing two 511 keV photons traveling in exactly opposite directions. The detector rings record both photons simultaneously — a *coincidence event* — telling us which straight line through the patient the emission occurred on, but not where along that line<sup>1</sup>. Tallying all such lines produces a **sinogram**, the raw projection data from which the activity image is reconstructed. However, some photons are absorbed (attenuated) by tissue before reaching the detectors — an effect that is more pronounced for deep structures and dense bone. **Attenuation correction** compensates for this, and it is what your pseudo-CT enables.

The PET reconstruction algorithm is titled Ordinary Poisson Ordered Subsets Expectation Maximization (OP-OSEM). OP-OSEM is simply a maximum likelihood expectation maximization algorithm where the observed data (the sinograms) are partioned and processed in chunks (subsets) to  dramatically accelerate the reconstruction speed. The reconstruction for BIC-MAC performs OP-OPSEM for `5 subsets x 4 iterations = 20 subiterations`. Three sinograms are used in OP-OSEM of which the last two must be attenuation corrected: `prompts_rd85[.s/.hs]`, `add_nac_rd85[.s/.hs]`, and `mult_nac_rd85[.s/.hs]` <sup>2</sup> 

The reconstruction pipeline for BIC-MAC (`src/recon/`) uses [STIR](http://stir.sourceforge.net/) (Software for Tomographic Image Reconstruction) to perform reconstruction. (see [Further Reading](#further-reading) for a primer)

> <sup>1</sup> A simplification. Modern scanners like the Siemens Quadra used for BIC-MAC can detect the tiny time of flight (TOF) difference between the arriving photons and infer their origin (roughly) on the LOR. Consequently, TOF-enabled sinograms have an extra dimension.
> 
> <sup>2</sup> RD85 means maximum ring difference of 85 and defines how oblique the LORs are allowed to be. The Siemens Quadra allows up to RD322, but RD85 was chosen to reduce the sinogram size.

---

## BIC-MAC reconstruction steps
Given a CT (ground-truth or pseudo-CT) and the subject's sinogram data (`recon/`), the the reconstruction pipeline exectutes the following steps to arrive at a reconstructed attenuation-corrected PET NIfTI image:

1. **Superimpose bed pixelated face** - The pseudo-CT face is replaced by a pre-saved pixelated face. Likewise, everything outside a ~1cm rim (pillows, bed, hair, air) is replaced by the ground truth image (see `ct_face_and_bed.nii.gz` and `face_and_bed_mask.nii.gz`). Consequently, the pseudo-CT algorithm will not benefit from trying to predict these areas. The `prediction_mask.nii.gz` under `ct-label` is the *inverse* mask of `face_and_bed_mask.nii.gz` and may be used to restrict training to the relevant body region. You can inspect the face-swapped intermediate file (`intermediates/ct_face_swapped.nii.gz`)

2. **HU → μ-map** — The pseudo-CT is a volume of Hounsfield units (HU), a relative X-ray density scale where air = −1000 and water = 0. PET reconstruction needs instead the linear attenuation coefficient μ (cm⁻¹) at the 511 keV photon energy of PET annihilation events. The conversion uses a bilinear model (Carney et al. 2006): a steep linear segment maps soft tissue (HU ≤ 47, dominated by water) and a flatter segment maps dense materials like bone (HU > 47). This means that errors in HU are not equally costly in soft tissue and in bone when converting to a mumap. You can inspect the mu_map intermediate file (`intermediates/mu_map.nii.gz`).

3. **Smooth μ-map** — A 4 mm FWHM Gaussian blur is applied to the μ-map before any sinogram operations. This is standard clinical practice to reduce the effect of CT noise and slight patient movement. Consequently, very fine structural detail in your pseudo-CT may be blurred away before it ever influences the PET sinogram, and by extension, the PET-based metrics. You can inspect the smoothed mu_map intermediate file (`intermediates/mumap_smoothed.nii.gz`)

4. **Resample to STIR** — Resamples the μ-map onto STIR's z-axis grid (ring spacing 3.29 mm), snapping the origin to the STIR coordinate system. A technical prerequisite for STIR's forward projection. The intermediate files are (`intermediates/mumap_stir[.hv/.ahv/.v]`)

5. **Compute ACF sinogram** — The μ-map is *forward projected* along every line of response (LOR) in the PET scanner geometry, computing the total integrated attenuation each annihilation photon pair experiences along that path. The result is the **attenuation correction factor (ACF)** sinogram, which has the same shape as the other PET sinograms, except from the fact that it lacks a TOF-dimension. The intermediate files are(`intermediates/acf[.hs/.v]`)

6.–7. **Apply ACF to sinograms** — The ACF sinogram is multiplied into both the *multiplicative* sinogram and the *additive* sinogram. This encodes the predicted attenuation into the reconstruction inputs. The intermediate files are (`intermediates/add[.hs/.s]` and `intermediates/mult[.hs/.s]`)

8. **OSEM reconstruction** — OP-OSEM is run for 20 subiterations (5 subsets, 4 iterations). The result is smoothed by a 4 mm FWHM Gaussian post-filter to reduce noise. The intermediate files are (`intermediates/pet_20[.ahv/.hv/.v]`).

9. **Convert to NIfTI** — Writes the reconstructed PET volume as a NIfTI file with the correct origin, accounting for the scanner bed position and gantry offset stored in `recon/offset.json`. The final output is `pet.nii.gz`

Intermediate outputs are written to `output_dir/intermediates/`. The pipeline skips steps whose outputs already exist, so it resumes automatically from a partial run.

---

> [!WARNING]
> Running reconstruction requires **~20 GB of RAM** and takes **20–120 minutes** depending on CPU speed. The `intermediates/` folder uses **~50 GB** of additional disk space — consider deleting it after a successful reconstruction.

## Running the Pipeline

### Option 1: Docker (recommended)

A pre-built image with STIR and all dependencies is available.

```bash
docker pull ghcr.io/bic-mac-challenge/recon:latest

docker run --rm \
  -v /path/to/sub-000/recon:/data/recon \
  -v /path/to/ct.nii.gz:/data/ct/ct.nii.gz \
  -v /path/to/output:/data/output \
  ghcr.io/bic-mac-challenge/recon:latest
```

The reconstructed PET is written to `/data/output/pet.nii.gz`. A full debug log is written to `/data/output/intermediates/recon.log`.

**Environment variables:**

| Variable | Default | Effect |
|----------|---------|--------|
| `OVERWRITE` | `0` | Set to `1` to ignore existing intermediates and rerun from scratch |
| `VERBOSE` | `0` | Set to `1` to stream STIR subprocess output to the terminal |

```bash
docker run --rm \
  -e OVERWRITE=1 \
  -e VERBOSE=1 \
  -v /path/to/sub-000/recon:/data/recon \
  -v /path/to/ct.nii.gz:/data/ct/ct.nii.gz \
  -v /path/to/output:/data/output \
  ghcr.io/bic-mac-challenge/recon:latest
```

### Option 2: Direct Python (requires local STIR)

```bash
python src/recon/main.py \
  --recon_dir <subject_recon_dir> \
  --ct <ct.nii.gz> \
  --output_dir <output_dir> \
  [-w] [-v]
```

| Argument | Description |
|----------|-------------|
| `--recon_dir` | Subject's `recon/` directory (contains sinograms, offset.json, face mask) |
| `--ct` | CT NIfTI file in Hounsfield units |
| `--output_dir` | Directory where `pet.nii.gz` and `intermediates/` will be written |
| `-w` / `--overwrite` | Rerun from scratch, ignoring existing intermediates |
| `-v` / `--verbose` | Stream STIR subprocess output to the terminal |

---

## Expected `recon/` Contents

```
recon/
├── add_nac_rd85.hs / .s       # additive sinogram (scatter + randoms)
├── mult_nac_rd85.hs / .s      # multiplicative sinogram (normalisation, decay)
├── prompts_rd85.hs / .s       # raw prompt coincidences
├── offset.json                # bed position and gantry offset
├── ct_face_and_bed.nii.gz     # ground-truth CT values at face + scanner bed
└── face_and_bed_mask.nii.gz   # binary face + scanner bed mask
```

Only subjects in the full `train/` split and the `val/` split include `recon/` data (see the main README).

---

## Further Reading

- Carney et al. (2006) — *"Method for Transforming CT Images for Attenuation Correction in PET/CT Scanners"*, Medical Physics. The bilinear HU→μ model used in this pipeline.
- Thielemans et al. (2012) — *"STIR: Software for Tomographic Image Reconstruction Release 2"*, Physics in Medicine and Biology. The reconstruction library used here.
