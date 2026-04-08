import stir
import os, subprocess
import logging

log = logging.getLogger('recon')
import nibabel as nib
import numpy as np
from scipy.ndimage import gaussian_filter


def validate_ct(pred_ct_path, ct_face_and_bed_path, hu_min_expected=-1024, hu_min_tolerance=20):
    """
    Validate the predicted CT against the ground-truth CT.

    Raises ValueError if shape or affine do not match.
    Prints a warning if the minimum HU value is unrealistically high
    (suggesting the prediction may not be in proper HU units).

    Parameters
    ----------
    pred_ct_path : str
        Path to the predicted CT NIfTI file.
    ct_face_and_bed_path : str
        Path to the ground-truth CT NIfTI file.
    hu_min_expected : float
        Expected minimum HU value (air/background), typically -1024.
    hu_min_tolerance : float
        How far the actual minimum may deviate from hu_min_expected before warning.
    """
    pred_img = nib.load(pred_ct_path)
    gt_img = nib.load(ct_face_and_bed_path)

    if pred_img.shape != gt_img.shape:
        raise ValueError(
            f"Shape mismatch between predicted CT {pred_img.shape} "
            f"and expected {gt_img.shape}."
        )

    if not np.allclose(pred_img.affine, gt_img.affine, atol=1e-3):
        raise ValueError(
            f"Affine mismatch between predicted CT {pred_img.affine} "
            f"and expected {gt_img.affine}."
        )

    pred_data = pred_img.get_fdata(dtype=np.float32)
    hu_min = pred_data.min()
    if hu_min > hu_min_expected + hu_min_tolerance:
        log.warning(
            f"Predicted CT minimum HU is {hu_min:.1f}. "
            f"Expected around {hu_min_expected} (air). "
            "The image may not be in correct HU units."
        )


def swap_face_from_gt(pred_ct_path, ct_face_and_bed_path, face_mask_path, output_path=None):
    """
    Replace the face region of a predicted CT with the ground-truth CT face.

    The face_mask defines which voxels belong to the face (non-zero = face).
    All three images must be in the same space/resolution.

    Parameters
    ----------
    pred_ct_path : str
        Path to the predicted CT NIfTI file.
    ct_face_and_bed_path : str
        Path to the ground-truth CT NIfTI file (ct.nii.gz).
    face_mask_path : str
        Path to the binary face mask NIfTI file.
    output_path : str, optional
        If given, the result is saved here. If None, the image is only returned.

    Returns
    -------
    nibabel.Nifti1Image
        CT with the face region swapped in from the ground truth.
    """
    pred_img = nib.load(pred_ct_path)
    gt_img = nib.load(ct_face_and_bed_path)
    mask_img = nib.load(face_mask_path)

    pred_data = pred_img.get_fdata(dtype=np.float32)
    gt_data = gt_img.get_fdata(dtype=np.float32)
    face_mask = mask_img.get_fdata() != 0

    if pred_data.shape != gt_data.shape or pred_data.shape != face_mask.shape:
        raise ValueError(
            f"Shape mismatch: pred {pred_data.shape}, gt {gt_data.shape}, "
            f"mask {face_mask.shape}. All inputs must be in the same space."
        )

    result = pred_data.copy()
    result[face_mask] = gt_data[face_mask]

    result_img = nib.Nifti1Image(result, pred_img.affine, pred_img.header)

    if output_path is not None:
        result_img.to_filename(output_path)
        log.debug(f"Face- and bed-swapped CT saved to {output_path}")

    return result_img

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

def smooth_image(img, fwhm_mm=4.0):
    """Gaussian smoothing. FWHM in mm, converted to sigma per axis."""
    voxel_sizes = img.header.get_zooms()[:3]
    sigma = [fwhm_mm / (2.355 * v) for v in voxel_sizes]
    smoothed = gaussian_filter(img.get_fdata(dtype=np.float32), sigma=sigma)
    return nib.Nifti1Image(smoothed, img.affine, img.header)

def save_stir_to_nifti(stir_img, output_path):
    stir.ITKOutputFileFormat().write_to_file(output_path, stir_img)

def calculate_acf(mumap_hv, reference_sinogram, output_hs, forwardprojector_par):
    with subprocess.Popen(
        ['stdbuf', '-oL', 'calculate_attenuation_coefficients', '--ACF', output_hs, mumap_hv, reference_sinogram, forwardprojector_par],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    ) as proc:
        for line in proc.stdout:
            line = line.rstrip()
            log.debug(line)
        proc.wait()
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, proc.args)
    

def mumap_to_stir(input_path, output_path, ring_spacing_mm=3.29114):
    """Zero origins, resample z to ring_spacing/2, snap z-origin to grid."""
    plane_sep = ring_spacing_mm / 2

    img = stir.FloatVoxelsOnCartesianGrid.read_from_file(input_path)
    img.set_origin(stir.FloatCartesianCoordinate3D(0.0, 0.0, 0.0))

    v = img.get_voxel_size()
    zoom_z = v.z() / plane_sep
    max_idx, min_idx = img.get_max_indices(), img.get_min_indices()
    z_size = int(max_idx[1] - min_idx[1]) + 1
    xy_size = int(max_idx[2] - min_idx[2]) + 1

    new_sizes = stir.Int3BasicCoordinate()
    new_sizes[1] = round(z_size * zoom_z)
    new_sizes[2] = xy_size
    new_sizes[3] = xy_size

    img_z = stir.zoom_image(img,
        stir.FloatCartesianCoordinate3D(zoom_z, 1.0, 1.0),
        stir.FloatCartesianCoordinate3D(0.0, 0.0, 0.0),
        new_sizes, stir.ZoomOptions(stir.ZoomOptions.preserve_values))

    o2 = img_z.get_origin()
    snapped_z = round(o2.z() / plane_sep) * plane_sep
    img_z.set_origin(stir.FloatCartesianCoordinate3D(snapped_z, 0.0, 0.0))

    stir.InterfileOutputFileFormat().write_to_file(output_path, img_z)
    log.debug(f"z-origin snapped: {o2.z():.4f} -> {snapped_z:.4f} mm, plane_sep={plane_sep:.5f} mm")


# def convert_ct_to_acf(ct_path, reference_sinogram, output_hs, ring_spacing_mm=3.29114,fwhm_mm=4.0,kvp=120):
#     root = os.path.dirname(output_hs)
#     mumap_hv_path = os.path.join(root, 'mumap_stir.hv')
#     mumap_nifti_path = os.path.join(root, 'mumap.nii.gz')
#     mumuap_smoothed_nifti_path = os.path.join(root, 'mumap_smoothed.nii.gz')
#     os.makedirs(root, exist_ok=True)

#     if not os.path.exists(mumap_nifti_path):
#         print("Converting CT HU to mu-map...")
#         mu = hu_to_mu(ct_path, kvp=kvp)
#         mu.to_filename(mumap_nifti_path)

#     if not os.path.exists(mumuap_smoothed_nifti_path):
#         print("Smoothing mu-map...")
#         mu_smoothed = smooth_image(mu, fwhm_mm=fwhm_mm)
#         mu_smoothed.to_filename(mumuap_smoothed_nifti_path)

#     if not os.path.exists(mumap_hv_path):
#         print("Converting mu-map to STIR format...")
#         mumap_to_stir(mumuap_smoothed_nifti_path, mumap_hv_path, ring_spacing_mm)
    
#     if not os.path.exists(output_hs):
#         print("Calculating ACF sinogram...")
#         calculate_acf(mumap_hv_path, reference_sinogram, output_hs)


# if __name__ == "__main__":
#     convert_ct_to_acf('ct.nii.gz', 'additive_term_SSRB.hs', 'outs/acf.hs', ring_spacing_mm=3.29114, debug=True)