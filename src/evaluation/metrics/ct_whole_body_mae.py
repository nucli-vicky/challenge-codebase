"""
Whole-body CT MAE metric.
"""

import numpy as np
import nibabel as nib

import numpy as np
import nibabel as nib

def hu_to_mu(ct_path):
    """
    Carney et al. 2006 (Med Phys 33:976-983) 
    bilinear HU to mu at 511 keV for 120 kVp.
    """
    # Carney parameters for KVP 120: (slope 'a', intercept 'b', breakpoint 'bp' in HU+1000)
    a, b, bp = (5.10e-5, 4.71e-2, 1047) # 1047 corresponds to 47 HU

    # Load NIfTI CT image
    ct = nib.load(ct_path)
    hu = ct.get_fdata(dtype=np.float32)

    # Pre-calculate (HU + 1000) to optimize the np.where evaluations
    hu1000 = hu + 1000

    # Apply the Carney bilinear scaling
    mu = np.where(hu1000 < bp,
                  9.6e-5 * hu1000,
                  a * hu1000 + b)

    # Ensure no negative attenuation values from extreme CT noise
    mu = np.clip(mu, 0, None)

    return nib.Nifti1Image(mu, ct.affine, ct.header)


LIVER_LABEL = 5


def compute_whole_body_mu_mae(
    pred_ct_path,
    gt_ct_path,
    body_seg_path,
    organ_seg_path,
    exclusion_cm=4.0,
    save_mask_path=None,
):
    """
    Voxel-wise MAE of linear attenuation coefficient (mu, cm⁻¹) across the body,
    excluding ±4 cm around the superior liver slice to reduce sensitivity to
    respiratory misalignment. CT images are converted from HU → mu using the
    Carney et al. 2006 bilinear model at 511 keV (see hu_to_mu()).

    Parameters
    ----------
    pred_ct_path   : path to predicted CT NIfTI (in HU)
    gt_ct_path     : ct-label/ct.nii.gz
    body_seg_path  : ct-label/body_seg.nii.gz
    organ_seg_path : ct-label/organ_seg.nii.gz
    save_mask_path : optional path to save the evaluation mask as NIfTI
    """

    # Convert CT → mu
    pred = hu_to_mu(pred_ct_path).get_fdata()
    gt = hu_to_mu(gt_ct_path).get_fdata()

    body_mask = nib.load(body_seg_path).get_fdata() > 0
    liver_mask = nib.load(organ_seg_path).get_fdata() == LIVER_LABEL

    # Slice thickness
    slice_thickness_mm = nib.load(pred_ct_path).header.get_zooms()[2]
    exclusion_slices = int(round((exclusion_cm * 10.0) / slice_thickness_mm))

    # Superior liver slice
    superior_slice = np.max(np.where(liver_mask)[2])

    z_min = max(0, superior_slice - exclusion_slices)
    z_max = min(pred.shape[2], superior_slice + exclusion_slices)

    exclusion_mask = np.zeros_like(body_mask, dtype=bool)
    exclusion_mask[:, :, z_min:z_max] = True

    eval_mask = body_mask & (~exclusion_mask)

    if save_mask_path is not None:
        nib.save(nib.Nifti1Image(eval_mask.astype(np.uint8), nib.load(body_seg_path).affine), save_mask_path)

    return np.mean(np.abs(pred - gt)[eval_mask])
