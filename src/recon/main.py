import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore", UserWarning)
    import stir
    from pet_recon import *
    from ct_to_acf import *
import argparse
import json
import logging
import nibabel as nib
import time
from utils import validate_recon_dir, check_input_hashes

log = logging.getLogger('recon')


def reconstruction_pipeline(
        output_dir,
        ct_path,
        ct_face_and_bed_path,
        face_and_bed_mask_path,
        add_sino_path,
        mult_sino_path,
        prompts_sino_path,
        offset_json_path,
        recon_template,
        acf_forwardprojector,
        overwrite=False,
        verbose=False,
        ):
    
    out_pet_nifti_path         = os.path.join(output_dir, 'pet.nii.gz')
    intemediates_dir           = os.path.join(output_dir, 'intermediates')
    face_swapped_ct_path       = os.path.join(intemediates_dir, 'ct_face_swapped.nii.gz')
    mumap_nifti_path           = os.path.join(intemediates_dir, 'mumap.nii.gz')
    mumuap_smoothed_nifti_path = os.path.join(intemediates_dir, 'mumap_smoothed.nii.gz')
    mumap_hv_path              = os.path.join(intemediates_dir, 'mumap_stir.hv')
    acf_sino_path              = os.path.join(intemediates_dir, "acf.hs")
    add_acf_sino_path          = os.path.join(intemediates_dir, "add.hs")
    mult_acf_sino_path         = os.path.join(intemediates_dir, "mult.hs")
    pet_hv_path                = os.path.join(intemediates_dir, "pet_20.hv")

    log.warning("Running reconstruction requires 20GB of RAM and takes 20-120 minutes depending on CPU. The intermediates folder takes up ~50GB of additional space. Consider deleting it after succesfull reconstruction")

    os.makedirs(intemediates_dir, exist_ok=True)

    check_input_hashes(intemediates_dir, ct_path, ct_face_and_bed_path, overwrite)

    t_total = time.perf_counter()

    # Step 1: Validate CT
    if not os.path.exists(face_swapped_ct_path) or overwrite:
        t = time.perf_counter()
        log.info("[1/10] Validating CT...")
        validate_ct(ct_path, ct_face_and_bed_path)
        log.info(f"      done ({time.perf_counter()-t:.1f}s)")
    else:
        log.info("[1/10] Validating CT (skipped)")

    # Step 2: Swap face and bed region
    if not os.path.exists(face_swapped_ct_path) or overwrite:
        t = time.perf_counter()
        log.info("[2/10] Swapping face and bed region...")
        swap_face_from_gt(ct_path, ct_face_and_bed_path, face_and_bed_mask_path, output_path=face_swapped_ct_path)
        log.info(f"      done ({time.perf_counter()-t:.1f}s)")
    else:
        log.info("[2/10] Swapping face and bed region (skipped)")
    ct_path = face_swapped_ct_path

    # Offset data needed to position the reconstructed PET image correctly in space
    with open(offset_json_path) as f:
        offsets = json.load(f)

    vertical_bed_start   = offsets['vertical_bed_start']
    horizontal_bed_start = offsets['horizontal_bed_start']
    gantry_offset        = offsets['gantry_offset']

    # Step 3: HU to mu-map
    if not os.path.exists(mumap_nifti_path) or overwrite:
        t = time.perf_counter()
        log.info("[3/10] Converting CT to mu-map...")
        mu = hu_to_mu(ct_path)
        mu.to_filename(mumap_nifti_path)
        log.info(f"      done ({time.perf_counter()-t:.1f}s)")
    else:
        log.info("[3/10] Converting CT to mu-map (skipped)")
        mu = nib.load(mumap_nifti_path)

    # Step 4: Smooth mu-map
    if not os.path.exists(mumuap_smoothed_nifti_path) or overwrite:
        t = time.perf_counter()
        log.info("[4/10] Smoothing mu-map...")
        mu_smoothed = smooth_image(mu, fwhm_mm=4)
        mu_smoothed.to_filename(mumuap_smoothed_nifti_path)
        log.info(f"      done ({time.perf_counter()-t:.1f}s)")
    else:
        log.info("[4/10] Smoothing mu-map (skipped)")

    # Step 5: Convert mu-map to STIR format
    if not os.path.exists(mumap_hv_path) or overwrite:
        t = time.perf_counter()
        log.info("[5/10] Converting mu-map image from NIfTI to STIR format...")
        mumap_to_stir(mumuap_smoothed_nifti_path, mumap_hv_path, ring_spacing_mm=3.29114)
        log.info(f"      done ({time.perf_counter()-t:.1f}s)")
    else:
        log.info("[5/10] Converting mu-map image from NIfTI to STIR format (skipped)")

    # Step 6: Calculate ACF sinogram
    if not os.path.exists(acf_sino_path) or overwrite:
        t = time.perf_counter()
        log.info("[6/10] Converting mu-map image to ACF sinogram...")
        calculate_acf(mumap_hv_path, add_sino_path, acf_sino_path, acf_forwardprojector)
        log.info(f"      done ({time.perf_counter()-t:.1f}s)")
    else:
        log.info("[6/10] Converting mu-map image to ACF sinogram (skipped)")

    # Step 7: Apply ACF to additive sinogram
    if not os.path.exists(add_acf_sino_path) or overwrite:
        t = time.perf_counter()
        log.info("[7/10] Multiplying ACF on additive NAC sinogram...")
        apply_acf_to_sinogram(add_sino_path, acf_sino_path, add_acf_sino_path)
        log.info(f"      done ({time.perf_counter()-t:.1f}s)")
    else:
        log.info("[7/10] Multiplying ACF on additive sinogram (skipped)")

    # Step 8: Apply ACF to multiplicative sinogram
    if not os.path.exists(mult_acf_sino_path) or overwrite:
        t = time.perf_counter()
        log.info("[8/10] Multiplying ACF on multiplicative sinogram...")
        apply_acf_to_sinogram(mult_sino_path, acf_sino_path, mult_acf_sino_path)
        log.info(f"      done ({time.perf_counter()-t:.1f}s)")
    else:
        log.info("[8/10] Multiplying ACF on multiplicative sinogram (skipped)")

    # Step 9: Reconstruct PET
    if not os.path.exists(pet_hv_path) or overwrite:
        t = time.perf_counter()
        log.info("[9/10] Reconstructing PET (this will take some time - especially first subiteration)...")
        run_reconstruction(recon_template, add_acf_sino_path, mult_acf_sino_path, prompts_sino_path, pet_hv_path, verbose=verbose)
        log.info(f"      done ({time.perf_counter()-t:.1f}s)")
    else:
        log.info("[9/10] Reconstructing PET (skipped)")

    # Step 10: Convert to NIfTI
    if not os.path.exists(out_pet_nifti_path) or overwrite:
        t = time.perf_counter()
        log.info("[10/10] Converting PET to NIfTI...")
        stir_pet_to_nifti(vertical_bed_start, horizontal_bed_start, gantry_offset, pet_hv_path, out_pet_nifti_path)
        log.info(f"       done ({time.perf_counter()-t:.1f}s)")
    else:
        log.info("[10/10] Converting PET to NIfTI (skipped)")

    log.info(f"Output PET saved to {out_pet_nifti_path}")
    log.info(f"Total time: {time.perf_counter()-t_total:.1f}s")


if __name__ == "__main__":
    _recon_dir = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(
        description="PET reconstruction pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Expected recon_dir contents:
  add_nac_rd85.hs, add_nac_rd85.s
  mult_nac_rd85.hs, mult_nac_rd85.s
  prompts_rd85.hs, prompts_rd85.s
  offset.json
  ct_face_and_bed.nii.gz
  face_and_bed_mask.nii.gz
""",
    )
    parser.add_argument("--recon_dir",   required=True, help="Reconstruction directory (e.g. /data/sub-000/recon)")
    parser.add_argument("--ct",          required=True, help="Input CT NIfTI file to be used for attenuation correction reconstruction")
    parser.add_argument("--output_dir",  required=True, help="Output directory; pet.nii.gz and intermediates/ will be written here")
    parser.add_argument("-w", "--overwrite", action="store_true", default=False, help="Overwrite existing intermediate and output files")
    parser.add_argument("-v", "--verbose", action="store_true", default=False, help="Show output from STIR subprocess calls")
    args = parser.parse_args()

    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(os.path.join(output_dir, "intermediates"), exist_ok=True)

    log.setLevel(logging.DEBUG)
    log.propagate = False

    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    console.setFormatter(logging.Formatter("%(message)s"))
    log.addHandler(console)

    file_handler = logging.FileHandler(os.path.join(output_dir, "intermediates", "recon.log"))
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    log.addHandler(file_handler)

    recon_dir = args.recon_dir
    validate_recon_dir(recon_dir)

    ct_face_and_bed_path      = os.path.join(recon_dir,  "ct_face_and_bed.nii.gz")
    face_and_bed_mask_path    = os.path.join(recon_dir,  "face_and_bed_mask.nii.gz")
    add_sino_path     = os.path.join(recon_dir,  "add_nac_rd85.hs")
    mult_sino_path    = os.path.join(recon_dir,  "mult_nac_rd85.hs")
    prompts_sino_path = os.path.join(recon_dir,  "prompts_rd85.hs")
    offset_json_path  = os.path.join(recon_dir,  "offset.json")
    recon_template       = os.path.join(_recon_dir, "recon_OSEM_template.par")
    acf_forwardprojector = os.path.join(_recon_dir, "acf_forwardprojector.par")

    log_path = os.path.join(output_dir, "intermediates", "recon.log")
    try:
        reconstruction_pipeline(
            output_dir=output_dir,
            ct_path=args.ct,
            ct_face_and_bed_path=ct_face_and_bed_path,
            face_and_bed_mask_path=face_and_bed_mask_path,
            add_sino_path=add_sino_path,
            mult_sino_path=mult_sino_path,
            prompts_sino_path=prompts_sino_path,
            offset_json_path=offset_json_path,
            recon_template=recon_template,
            acf_forwardprojector=acf_forwardprojector,
            overwrite=args.overwrite,
            verbose=args.verbose,
        )
    except Exception:
        if not args.verbose:
            print(f"\nReconstruction failed. Check {log_path} for details, or rerun with -v/--verbose for more output.", flush=True)
        raise
