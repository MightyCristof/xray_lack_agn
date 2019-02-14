PRO xray_lacking_agn, FIELD = field, $
                      DETECT = detect, $
                      SOFTX = softx, $
                      FLIM = flim, $
                      LUM = lum


;; re-sample sources only in X-ray fields
if keyword_set(field) then begin
    load_vars,'fits.sav','_fits'
    field_nustar
    field_xmm
    field_chandra
    load_vars,'xfield_nst.sav','_xray_nst'
    load_vars,'xfield_xmm.sav','_xray_xmm'
    load_vars,'xfield_cha.sav','_xray_cha'
    agn_xray_infield
endif

;; load X-ray infield results
load_vars,'infield_fits.sav','_inf_fits'
load_vars,'infield_nst.sav','_inf_nst'
load_vars,'infield_xmm.sav','_inf_xmm'
load_vars,'infield_cha.sav','_inf_cha'

if keyword_set(detect) then begin
    ;; flag infield sample sources with X-ray detections
    detect_nustar
    detect_xmm
    detect_chandra
    detect_wac
endif

;; load X-ray detection results
load_vars,'xdetect_nst.sav','_det_nst'
load_vars,'xdetect_xmm.sav','_det_xmm'
load_vars,'xdetect_cha.sav','_det_cha'
;; load WISE AGN Catalog matches
load_vars,'xdetect_wac.sav','_det_wac'

;; X-ray analysis at 2-10keV
if keyword_set(softx) then convert_xray2_10kev,1.4,/plt
load_vars,'xray_soft210.sav','_soft210'
;; X-ray analysis at 0.5-2keV
if keyword_set(softx) then convert_xray05_2kev,/plt
load_vars,'xray_soft052.sav','_soft052'

;; estimate X-ray flux limits from catalogs
if keyword_set(flim) then xcat_fluxlim,/plt
load_vars,'xray_flx_lim.sav','_flx_lim'



if keyword_set(lum) then begin
    ;; load SED template components
    load_comp,'comp_sed4.sav',/push
    ;; calculate IR luminosities and X-ray conversions
    agn_xray_lumin,/dered
endif
load_vars,'xray_lum.sav','_agn_lum'







END





;    common _inf_fits
;    common _inf_nst
;    common _inf_xmm
;    common _inf_cha
;    common _det_nst
;    common _det_xmm
;    common _det_cha
;    common _det_wac
;    common _soft210
;    common _soft052
;    common _flx_lim
;    common _agn_lum



