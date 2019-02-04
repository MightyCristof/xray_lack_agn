PRO xray_lacking_agn, prep = prep


    ;; load SED fit output
if keyword_set(prep) then load_vars,'fits.sav','_fits'

    ;; sample sources in X-ray fields
if keyword_set(prep) then xray_lacking_agn_fld
    ;; load X-ray field results
if keyword_set(prep) then load_vars,'xfield_nst.sav','_xray_nst'
if keyword_set(prep) then load_vars,'xfield_xmm.sav','_xray_xmm'
if keyword_set(prep) then load_vars,'xfield_cha.sav','_xray_cha'

    ;; flag and subset sample sources in X-ray fields
if keyword_set(prep) then xray_lacking_agn_inf
;; load X-ray infield results
load_vars,'infield_fits.sav','_inf_fits'
load_vars,'infield_nst.sav','_inf_nst'
load_vars,'infield_xmm.sav','_inf_xmm'
load_vars,'infield_cha.sav','_inf_cha'

;; flag infield sample sources with X-ray detections
if keyword_set(prep) then xray_lacking_agn_det
;; load X-ray detection results
load_vars,'xdetect_nst.sav','_det_nst'
load_vars,'xdetect_xmm.sav','_det_xmm'
load_vars,'xdetect_cha.sav','_det_cha'
;; load WISE AGN Catalog matches
load_vars,'xdetect_wac.sav','_det_wac'

;; X-ray analysis at 2-10keV
if keyword_set(prep) then xray_lacking_agn_soft210,1.4
;; load 2-10keV results
load_vars,'xray_soft210.sav','_soft210'
;; X-ray analysis at 0.5-2keV
if keyword_set(prep) then xray_lacking_agn_soft052
;; load 0.5-2keV results
load_vars,'xray_soft052.sav','_soft052'

;; load SED template components
load_comp,'comp_sed4.sav',/push
xray_lacking_agn_lum

;; Candidate AGN by X-ray source detection/limit
load_vars,'xray_lum.sav','_agn_lum'







END





;common _inf_fits
;common _inf_nst
;common _inf_xmm
;common _inf_cha
;common _det_nst
;common _det_xmm
;common _det_cha
;common _det_wac
;common _soft210
;common _soft052
;common _comp
;common _agn_lum
