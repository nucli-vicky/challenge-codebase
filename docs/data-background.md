# Data Background

Understanding how the input features were acquired — and crucially, *when* — helps you reason about their alignment with the ground-truth CT and each other. All participants are healthy controls. 

## Acquisition timeline

All data was acquired on two scanners at Rigshospitalet, Copenhagen:

- **PET/CT scanner**: Siemens Biograph Vision Quadra (long axial field-of-view)
- **MRI scanner**: Siemens MAGNETOM Vida 3T (separate room, different scanner)

A typical session ran as follows:

```
  t = -2 min   Topogram (2D scout X-ray, ~2 min before CT, duration ~7s)
  t =  0 min   Low-dose CT acquisition (duration ~7s)
  t =  2 min   18F-FDG injection
  t =  2-72 min  Dynamic PET acquisition (70 minutes continuous)
               └─ NAC-PET and PET reconstructed from t = 50–70 min window
  (same day)   MRI on a separate scanner (duration ~25s per chunk)
```

## Feature modalities and alignment

The three input modalities span a spectrum from *well-aligned but low-quality* to *high-quality but poorly aligned* relative to the ground-truth CT. All images under `features/` have been resampled to the CT grid (512×512×531). The topogram is 2D and therefore only resampled in the frontal plane (512×1×531).

### Topogram (`topogram.nii.gz`)
A 2D sagittal scout X-ray acquired roughly 2 minutes before the CT. It is the most temporally and spatially aligned feature, but has very limited anatomical detail — it is a single 2D projection with coarse resolution. The topogram and CT share the same bed position and scanner table.

### NAC-PET (`nacpet.nii.gz`)
The non-attenuation-corrected PET, reconstructed from the **50–70 minute** window of the dynamic scan. This is the same scanner as the CT, so bed positioning is roughly preserved, but ~50 minutes elapsed between the CT and the PET window. Over this time, gross body position is stable, but peripheral structures (arms, fingers, legs) may have shifted a bit. NAC-PET encodes useful soft-tissue contrast and body shape, and has the same matrix and affine as the ground-truth CT. 

### MRI (`mri_chunk_*` / `mri_combined_*`)
Whole-body DIXON MRI acquired on a **separate scanner** the same day as the PET/CT (for 95/99 participants; 4 were scanned 9–89 days later). The MRI provides excellent soft-tissue contrast and anatomical detail, but it is not inherently aligned to PET/CT space. The DIXON MRI sequence used is designed specifically to capture anatomical informing relevant for attenuation correction. **The whole-body image was acquired in four sequential scans of 25s coressponding to four body chunks**. Prior to each acquisition the participant was asked to take a deep breath and hold still. Since participants were not asked to hold their breath during Topogram, CT or PET acquisition, the lung volume is often larger on the MRI. Participants were placed head first supine with arms down the side in both PET/CT and MRI, but the use of coils, different beds and raisable head rests makes the MRI-to-PET/CT alignment poor. A rigid translation (no rotation) was applied to bring the MRI into approximate PET/CT alignment. Whether to incorporate additional registration as a preprocessing step is left to the participant.



### Summary

| Feature | Alignment to CT | Anatomical quality |
|---|---|---|
| Topogram | Excellent (same scan, ~2 min prior) | Low (2D projection) |
| NAC-PET | Good (same scanner, ~50 min later) | Moderate |
| MRI | Approximate (rigid translation only) | High |

 Despite all images sharing the same matrix, this does not mean they are perfectly registered.
