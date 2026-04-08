[WIP]

# PET Imaging Background

This guide is for participants who are familiar with medical imaging (MRI, CT) but have not worked extensively with PET. It covers the concepts you need to understand the challenge task, the data, and the evaluation metrics.

---

## What is PET?

**Positron Emission Tomography (PET)** is a functional imaging modality. Unlike CT or MRI, which image anatomy, PET images *metabolic activity*: how actively tissues are taking up a radiotracer.

The most common radiotracer is **FDG** (fluorodeoxyglucose), a glucose analogue. Tissues with high metabolic activity (tumors, brain, heart) accumulate more FDG. The radiotracer emits positrons, which annihilate with electrons to produce two 511 keV gamma rays traveling in opposite directions. Coincident detection of these pairs is what forms the signal.

The raw detector data — called a **sinogram** — captures projection counts at different angles and offsets, analogous to a CT sinogram. From this, an image is reconstructed (usually with iterative algorithms like **OSEM**).

---

## The Attenuation Problem

511 keV photons passing through tissue are attenuated (absorbed or scattered). Without correcting for this, tissues deep inside the body appear artificially dim — the reconstructed image would be quantitatively wrong.

**Attenuation correction (AC)** compensates for this by estimating how much signal was lost along each line of response. The correction factor for each line depends on the total attenuation integral through the body along that path.

In clinical **PET/CT** scanners, a CT scan acquired immediately before the PET scan provides the attenuation map:

1. Convert CT Hounsfield units → **linear attenuation coefficients at 511 keV** (the μ-map)
2. Forward-project the μ-map to compute the **Attenuation Correction Factor (ACF)** sinogram
3. Apply the ACF to the raw PET sinogram before reconstruction

CT is the primary source of radiation dose in a PET/CT exam. Eliminating it by predicting a pseudo-CT from non-ionizing inputs (NAC-PET, MRI) is particularly important for radiation-sensitive populations: children, pregnant patients, and patients requiring frequent follow-up scans.

---

## HU → μ Conversion

Hounsfield units (HU) encode X-ray attenuation relative to water (water = 0 HU, air = −1000 HU). The conversion to linear attenuation coefficients at 511 keV follows a **bilinear model** (Carney et al. 2006):

| Tissue | HU range | Formula |
|--------|----------|---------|
| Air / soft tissue | HU ≤ 0 | μ = 9.6 × 10⁻⁵ × (HU + 1000) |
| Bone | HU > 0 | μ = 9.6 × 10⁻⁵ × 1000 + bone_slope × HU |

The bone slope depends on X-ray tube voltage (kVp). This challenge uses 120 kVp (bone_slope = 5.10 × 10⁻⁵ cm⁻¹/HU).

The resulting μ-map has units of cm⁻¹ and is used directly in the reconstruction pipeline.

---

## What is a Sinogram?

A sinogram stores the raw measured (or corrected) projection data from the PET detector ring. Each row corresponds to a different angular view; each column to a different radial offset. Together they encode the line-integral of activity along every line of response sampled by the detector.

The challenge dataset provides three sinograms per subject (under `recon/`):

| File | Contents |
|------|---------|
| `prompts_rd85.*` | Raw prompt coincidences (signal + background) |
| `mult_nac_rd85.*` | Multiplicative corrections (normalization, detector efficiency, decay) |
| `add_nac_rd85.*` | Additive background estimate (scatter + randoms) |

You do **not** need to work with sinograms directly — the reconstruction pipeline (`src/recon/`) handles everything. But understanding that your pseudo-CT affects the ACF, which is applied to these sinograms before reconstruction, explains why CT accuracy matters for PET quality.

---

## Standardized Uptake Value (SUV)

Raw PET voxel values are proportional to activity concentration (MBq/mL) but vary with injected dose and patient weight — making cross-patient comparisons difficult. **SUV** normalizes for this:

$$\text{SUV} = \frac{\text{voxel activity concentration [kBq/mL]}}{\text{injected dose [kBq] / body weight [g]}}$$

In this challenge, injected dose metadata is not available. Instead, SUV is estimated using the total PET signal and the body mask volume as a weight proxy — consistent across all comparisons.

A perfect pseudo-CT that introduces no attenuation error would produce the same SUV distribution as the ground-truth CT-AC PET. Errors in the μ-map cause regional SUV biases, which the evaluation metrics quantify.

---

## Why MRI + Topogram?

Predicting a CT-quality attenuation map from PET alone is difficult — the NAC-PET has poor tissue contrast and geometric distortion from attenuation effects. MRI and the topogram add complementary anatomical information:

- **DIXON MRI**: A fat/water separation sequence. The in-phase and out-of-phase images allow segmentation of fat vs. soft tissue, which is critical for accurate μ values in adipose-rich regions. Four bed positions (chunks) cover the whole body.
- **Topogram (scout)**: A 2D projection radiograph (like a low-dose planar X-ray) acquired before the main CT. It shows the patient silhouette and bone structure from one projection angle.

Participants are expected to incorporate all available modalities. The baseline model uses NAC-PET only and serves as a lower bound.

---

## The Scanner Setup

Subjects in this challenge were scanned on:

- **Siemens Biograph Vision Quadra** (PET/CT): A long axial field of view (LAFOV) PET scanner with a 106 cm detector ring — it can image the full body in a single bed position but uses multiple bed positions for optimal sensitivity. Ring spacing: 3.29114 mm.
- **Siemens MAGNETOM Vida** (3T MRI): Standard clinical 3T scanner. The DIXON sequence acquires four echoes per TR, enabling fat/water separation from the phase difference of in-phase and out-of-phase images.

Both scanners produce images in the same physical coordinate frame after registration. All data in the dataset has been resampled to the CT grid for consistency.

---

## Summary: What You Are Predicting

Your model receives (per subject):

- `features/nacpet.nii.gz` — NAC-PET volume (low tissue contrast, correlated with uptake)
- `features/mri_chunk_*_*.nii.gz` — DIXON MRI bed positions (good soft tissue contrast)
- `features/mri_combined_*.nii.gz` — Stitched whole-body DIXON
- `features/topogram.nii.gz` — 2D scout image
- `features/metadata.json` — sex, age, height, weight

And must output:

- `ct.nii.gz` — pseudo-CT in Hounsfield units, same shape and affine as the input NAC-PET

This pseudo-CT is then fed into the reconstruction pipeline, which produces an AC-corrected PET image. Both the CT and PET outputs are evaluated against ground truth.

---

## Further Reading

- Carney et al. (2006) — *"Method for Transforming CT Images for Attenuation Correction in PET/CT Scanners"*, Medical Physics. The bilinear HU→μ model used in this challenge.
- Townsend (2008) — *"Multimodality Imaging of Structure and Function"*, Physics in Medicine and Biology. Good clinical PET/CT overview.
- Thielemans et al. (2012) — *"STIR: Software for Tomographic Image Reconstruction Release 2"*, Physics in Medicine and Biology. The reconstruction library used in this challenge.
