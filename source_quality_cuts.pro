PRO source_quality_cuts


common _fits        
common _inf_cha     
common _inf_xmm     
common _inf_nst     
common _det_cha     
common _det_xmm     
common _det_nst     
common _det_wac     
common _softx       
common _fluxlim     
common _comp        
common _agn_lum     
common _clean_cha   
common _clean_xmm   
common _clean_nst   



;;----------------------------------------------------------------------------------------
;; Main Sample Quality Control 
;;----------------------------------------------------------------------------------------
;; redshift range 0 � z ��1
iiz = z gt 0. and z lt 0.6
;; SED chi-square goodness of fit
chi = reform(param[-2,*])
dof = reform(param[-1,*])
rchi = chi/dof
iichi = rchi le 20.

;; indices of WISE photometry
iwise = where(strmatch(band,'WISE?',/fold))
;; ensure WISE photometry exists
sn_wise = flux[iwise,*]/e_flux[iwise,*]
totsn = total(sn_wise ge 1.,1)          ;; all WISE photometry must exist and S/N � 1
iisn = totsn ge 4.                      ;; note: all non-finite sn_wise == -NaN

;; separate IR bright and weak sources
iiirb = lir ge 42.
iiirw = lir gt 0. and lir lt 42.
;; require WISE bands for high E(B-V)
iiirc = (ebv gt 50. and total(bin[iwise[2:3],*],1) eq 2.) or (ebv lt 50. and total(bin[iwise[2:3],*],1) ge 1)
;; constrain E(B-V) extrema to high AGN fraction
iiebv = (ebv lt 0.1 and agnf15.obs gt 0.7) or (ebv gt 0.1 and ebv lt 50.) or (ebv gt 50. and agnf15.obs gt 0.7)

;; passes all quality cuts, IR bright, and constrained E(B-V)
iiagn = iiz and iichi and iisn and iiirb and iiirc and iiebv

sav_vars = ['CHI','DOF','RCHI','SN_WISE','TOTSN']
sav_inds = ['IIZ','IICHI','IWISE','IISN','IIIRB','IIIRW','IIIRC','IIEBV','IIAGN']

;;----------------------------------------------------------------------------------------
;; X-ray Catalog Quality Control
;;----------------------------------------------------------------------------------------
;; combined raw X-ray detections, sample-catalog matches
det_str = strjoin('IIDET'+xfield,' or ')
re = execute('iidet = '+det_str)
;; combined cleaned X-ray detections, sample-cleaned catalog matches
clean_str = strjoin('IICLEAN'+xfield,' or ')
re = execute('iiclean = '+clean_str)

;; S/N cut on X-ray fluxes
for f = 0,nfield-1 do begin
    re = execute('SN'+xfield[f]+' = '+xray_flx[f]+'/'+xray_err[f])
    re = execute('IISN'+xfield[f]+' = SN'+xfield[f]+' ge 2.')
endfor

sav_vars = [sav_vars,'SN_CHA','SN_XMM','SN_NST']
sav_inds = [sav_inds,'IIDET','IICLEAN','IISN_CHA','IISN_XMM','IISN_NST']


;;----------------------------------------------------------------------------------------
;; Sources above the flux limit
;;----------------------------------------------------------------------------------------
iiflim_pass = 'IIFLIM_PASS'+xfield
iiflim_fail = 'IIFLIM_FAIL'+xfield
for f = 0,nfield-1 do begin
    re = execute('IIFLIM_PASS'+xfield[f]+' = iilir and FLIM'+xfield[f]+' gt 0. and fxir ge FLIM'+xfield[f])
    re = execute('IIFLIM_FAIL'+xfield[f]+' = iilir and FLIM'+xfield[f]+' gt 0. and fxir lt FLIM'+xfield[f])
endfor


sav_vars = [sav_vars]
sav_inds = [sav_inds,iiflim_pass,iiflim_fail]


;;----------------------------------------------------------------------------------------
;; Final output indices
;;----------------------------------------------------------------------------------------
;; AGN with X-ray detections/non-detections
iiagn_det = 'IIAGN_DET'+xfield
iiagn_drm = 'IIAGN_DRM'+xfield
iiagn_non = 'IIAGN_NON'+xfield
iiagn_nrm = 'IIAGN_NRM'+xfield
for f = 0,nfield-1 do begin
    re = execute(iiagn_det[f]+' = IIINF'+xfield[f]+' and IIAGN and IIDET'+xfield[f]+' and IICLEAN'+xfield[f]+' and IISN'+xfield[f])
    re = execute(iiagn_drm[f]+' = IIINF'+xfield[f]+' and IIAGN and IIDET'+xfield[f]+' and IICLEAN'+xfield[f]+' and ~IISN'+xfield[f])
    re = execute(iiagn_non[f]+' = IIINF'+xfield[f]+' and IIAGN and ~IIDET'+xfield[f]+' and IIFLIM_PASS'+xfield[f]+' and '+sdst[f]+' le 0.5*FOV'+xfield[f])
    re = execute(iiagn_nrm[f]+' = IIINF'+xfield[f]+' and IIAGN and ~IIDET'+xfield[f]+' and IIFLIM_FAIL'+xfield[f]+' or '+sdst[f]+' gt 0.5*FOV'+xfield[f])
endfor

sav_vars = [sav_vars]
sav_inds = [sav_inds,iiagn_det,iiagn_drm,iiagn_non,iiagn_nrm]

;; combined sources
iiagn_det = iiagn and iiclean
iiagn_non = iiagn and ~iidet
iiagn_rem = iiagn and iidet and ~iiclean

sav_vars = [sav_vars]
sav_inds = [sav_inds,'IIAGN_DET','IIAGN_NON','IIAGN_REM']


sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="quality_src.sav"')


END







