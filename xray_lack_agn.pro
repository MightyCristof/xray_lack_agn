PRO xray_lack_agn, subdir, $
                   INFIELD = infield, $
                   DETECT = detect, $
                   CONVERT = convert, $
                   FXLIM = fxlim, $
                   AGNLUM = agnlum, $
;                   CLEAN = clean, $
                   QUALITY = quality, $
                   COMBINE = combine, $
                   NHDIST = nhdist, $
                   SURV = surv
                   

;; check for keyword commands
nkeys = n_elements(infield) + $
        n_elements(detect) + $
        n_elements(convert) + $
        n_elements(fxlim) + $
        n_elements(agnlum) + $
        n_elements(clean) + $
        n_elements(quality) + $
        n_elements(combine) + $
        n_elements(nhdist) + $
        n_elements(surv)
if (nkeys eq 0) then GOTO, NO_KEYS

;; assume current directory unless specified
if (n_elements(subdir) eq 0) then path = './' else $
                                  path = subdir+'/' & file_mkdir,path

;; print LX-LIR information to screen
if keyword_set(agnlum) then begin
    rel = 'F09'
    if (strmatch(strupcase(agnlum),'C17')) then rel = 'C17'
    print, ''
    print, '***********************************************'
    print, '********     LX-LIR RELATION: '+rel+'      ********'
    print, '***********************************************'
    print, ''
endif

;; load SED output and template components
load_vars,'fits.sav','_fits'
load_vars,'resamp.sav','_resamp'
load_comp,'../data_prep/comp*.sav'

;; directory for output
pushd,path

;; flag X-ray footprints
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

;; flag X-ray detections and WISE AGN Catalog
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

;; convert X-ray fluxes to 2-10keV
if keyword_set(convert) then begin
    gmma = 1.8
    convert_xband,gmma
    nkeys--
endif
load_vars,'xband_conversions.sav','_xconv'
if (nkeys eq 0) then GOTO, NO_KEYS

;; estimate X-ray flux limits from catalogs
if keyword_set(fxlim) then begin
    catalog_fxlim,/multi_sn
    nkeys--
endif
load_vars,'catalog_flux_limits.sav','_fxlim'
if (nkeys eq 0) then GOTO, NO_KEYS

;; calculate IR and X-ray luminosities and luminosity ratios
if keyword_set(agnlum) then begin
    agn_luminosities,/dered,rel=rel
    nkeys--
endif
load_vars,'src_luminosities.sav','_agnlum'
if (nkeys eq 0) then GOTO, NO_KEYS

;;; clean X-ray sources
;if keyword_set(clean) then begin
;    clean_source_chandra
;    clean_source_xmm
;    clean_source_nustar
;    nkeys--
;endif
;load_vars,'cleaned_cha.sav','_clean_cha'
;load_vars,'cleaned_xmm.sav','_clean_xmm'
;load_vars,'cleaned_nst.sav','_clean_nst'
;if (nkeys eq 0) then GOTO, NO_KEYS

;; pass quality cuts for analysis set
if keyword_set(quality) then begin
    source_quality_cuts
    nkeys--
endif
load_vars,'quality_src.sav','_quality'
if (nkeys eq 0) then GOTO, NO_KEYS

;; combine values by source
if keyword_set(combine) then begin
    combine_luminosities
    nkeys--
endif
load_vars,'combined_lum.sav','_combined'
if (nkeys eq 0) then GOTO, NO_KEYS

;; estimate NH distribution
if keyword_set(nhdist) then begin
    compute_nh_distribution
    nkeys--
endif
load_vars,'nh_dist.sav','_nhdist'
if (nkeys eq 0) then GOTO, NO_KEYS

;; run survival analysis on analysis set
if keyword_set(surv) then begin
    surv_analysis;,/asurv
    nkeys--        
endif


NO_KEYS:
popd



END






