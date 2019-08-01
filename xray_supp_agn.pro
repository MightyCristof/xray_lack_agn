PRO xray_supp_agn, INFIELD = infield, $
                   DETECT = detect, $
                   SOFTX = softx, $
                   FLIM = flim, $
                   LUM = lum, $
                   CLEAN = clean, $
                   QUALITY = quality, $
                   RATIO = ratio, $
                   NHDIST = nhdist


nkeys = n_elements(infield) + $
        n_elements(detect) + $
        n_elements(softx) + $
        n_elements(flim) + $
        n_elements(lum) + $
        n_elements(clean) + $
        n_elements(quality) + $
        n_elements(ratio)
if (nkeys eq 0) then GOTO, NO_KEYS

;; load SED output
load_vars,'fits.sav','_fits'

;; re-sample sources only in X-ray fields
;; match sample to WISE AGN Catalog
if keyword_set(infield) then begin
    infield_chandra
    infield_xmm_newton
    infield_nustar
    nkeys--
endif
load_vars,'infield_cha.sav','_inf_cha'
load_vars,'infield_xmm.sav','_inf_xmm'
load_vars,'infield_nst.sav','_inf_nst'
if (nkeys eq 0) then GOTO, NO_KEYS

;; flag X-ray detections
if keyword_set(detect) then begin
    detect_chandra
    detect_xmm
    detect_nustar
    detect_wac
    nkeys--
endif
load_vars,'detections_cha.sav','_det_cha'
load_vars,'detections_xmm.sav','_det_xmm'
load_vars,'detections_nst.sav','_det_nst'
load_vars,'detections_wac.sav','_det_wac'
if (nkeys eq 0) then GOTO, NO_KEYS

;; X-ray conversion to 2-10keV
if keyword_set(softx) then begin
    gmma = 1.8
    convert_softx,gmma;,/plt
    nkeys--
endif
load_vars,'softx_conversions.sav','_softx'
if (nkeys eq 0) then GOTO, NO_KEYS

;; estimate X-ray flux limits from catalogs
if keyword_set(flim) then begin
    fit_catalog_xlim;,/plt
    nkeys--
endif
load_vars,'catalog_flux_limits.sav','_fluxlim'
if (nkeys eq 0) then GOTO, NO_KEYS

;; calculate IR luminosities and X-ray conversions
;; load SED template components
load_comp,'../comp_*.sav'
if keyword_set(lum) then begin
    agn_xray_luminosity,/dered,/combine
    nkeys--
endif
load_vars,'src_luminosity.sav','_agnlum'
if (nkeys eq 0) then GOTO, NO_KEYS

;; clean X-ray sources
if keyword_set(clean) then begin
    clean_source_chandra
    clean_source_xmm
    clean_source_nustar
    nkeys--
endif
load_vars,'cleaned_cha.sav','_clean_cha'
load_vars,'cleaned_xmm.sav','_clean_xmm'
load_vars,'cleaned_nst.sav','_clean_nst'
if (nkeys eq 0) then GOTO, NO_KEYS

;; pass quality control
if keyword_set(quality) then begin
    source_quality_cuts
    nkeys--
endif
load_vars,'quality_src.sav','_quality'
if (nkeys eq 0) then GOTO, NO_KEYS

;; luminosity ratio LX/LXIR
if keyword_set(ratio) then begin
    compute_luminosity_ratio
    nkeys--
endif
load_vars,'lum_ratio.sav','_lum_ratio'
if (nkeys eq 0) then GOTO, NO_KEYS

;; raw NH distribution
if keyword_set(nhdist) then begin
    compute_nh_distribution
    nkeys--
endif
load_Vars,'nh_dist.sav','_nhdist'
if (nkeys eq 0) then GOTO, NO_KEYS


NO_KEYS:
END






