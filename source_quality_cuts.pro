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
common _agnlum     
common _clean_cha   
common _clean_xmm   
common _clean_nst   



;;----------------------------------------------------------------------------------------
;; Main Sample Quality Control 
;;----------------------------------------------------------------------------------------
;; redshift range 0 < z ²Ê0.6 for photo-zs
;; redshift range 0 < z ²Ê1.o for spec-zs
;; moved to PRE-FITTING
;ztype = zorig(zarr)
;iizp = z gt 0. and z le 0.6 and strmatch(ztype,'ZP')
;iizs = z gt 0. and z le 1.0 and strmatch(ztype,'ZS*')
;iiz = iizp or iizs
;; SED chi-square goodness of fit
chi = reform(param[-2,*])
dof = reform(param[-1,*])
rchi = chi/dof
iichi = rchi le 20.

;; ensure WISE photometry exists
;; MOVED TO PRE-FITTING
;; indices of WISE photometry
;iwise = where(strmatch(band,'WISE?',/fold))
;sn_wise = flux[iwise,*]/e_flux[iwise,*]
;totsn = total(sn_wise ge 1.,1)          ;; all WISE photometry must exist and S/N ³ 1
;iisnw = totsn ge 4.                      ;; note: all non-finite sn_wise == -NaN

;; separate IR bright and weak sources
iiirb = lir ge 42.
iiirw = lir gt 0. and lir lt 42.
;; require WISE bands for high E(B-V)
;iiirc = (ebv gt 50. and total(bin[iwise[2:3],*],1) eq 2.) or (ebv lt 50. and total(bin[iwise[2:3],*],1) ge 1)
;; constrain E(B-V) extremum to high AGN fraction
;iiebv = (ebv lt 0.1 and agnf15.obs gt 0.7) or (ebv gt 0.1 and ebv lt 50.) or (ebv gt 50. and agnf15.obs gt 0.7)
;; AGN fraction minimum of 70%
iifagn = agnf15.obs ge 0.7

;; passes all quality cuts, IR bright, and constrained E(B-V)
;iiagn = iiz and iichi and iisnw and iiirb and iifagn;iiirc and iiebv and iifagn
iiagn = iichi and iiirb and iifagn

;sav_vars = ['ZTYPE','CHI','DOF','RCHI','SN_WISE','TOTSN']
;sav_inds = ['IIZP','IIZS','IIZ','IICHI','IWISE','IISNW','IIIRB','IIIRW','IIFAGN','IIAGN'];'IIIRC','IIEBV','IIAGN','IIFAGN']
sav_vars = ['CHI','DOF','RCHI']
sav_inds = ['IICHI','IIIRB','IIIRW','IIFAGN','IIAGN'];'IIIRC','IIEBV','IIAGN','IIFAGN']

;;----------------------------------------------------------------------------------------
;; X-ray Catalog Quality Control
;;----------------------------------------------------------------------------------------
;; shorthand DET and CLEAN
iidet = 'IIDET'+xfield
iiclean = 'IICLEAN'+xfield

sn = 'SN'+xfield
iisn = 'IISN'+xfield
;; S/N cut on X-ray fluxes
for f = 0,nfield-1 do begin
    re = execute(sn[f]+' = FLX'+xfield[f]+'/ERR'+xfield[f])
    re = execute(iisn[f]+' = SN'+xfield[f]+' ge 3.')
endfor

sav_vars = [sav_vars,sn]
sav_inds = [sav_inds,iisn]


;; sources above the flux limit
;; requires a separate fail flag to account for FLIM exists (iifail =/= ~iipass)
iiflim_pass = 'IIFLIM_PASS'+xfield
iiflim_fail = 'IIFLIM_FAIL'+xfield
for f = 0,nfield-1 do begin
    re = execute(iiflim_pass[f]+' = iilir and FLIM'+xfield[f]+' gt 0. and fxir ge FLIM'+xfield[f])
    re = execute(iiflim_fail[f]+' = iilir and FLIM'+xfield[f]+' gt 0. and fxir lt FLIM'+xfield[f])
endfor

;; sources within 1/2 FOV
;; requires a separate fail flag to account for INF (iioffa =/= ~iicntr)
;; FOV 70% effective area
iicntr = 'IICNTR'+xfield
iioffa = 'IIOFFA'+xfield
for f = 0,nfield-1 do begin 
    re = execute(iicntr[f]+' = IIINF'+xfield[f]+' and SDST'+xfield[f]+' le EFF'+xfield[f])
    re = execute(iioffa[f]+' = IIINF'+xfield[f]+' and SDST'+xfield[f]+' gt EFF'+xfield[f])
endfor

sav_vars = [sav_vars]
sav_inds = [sav_inds,iiflim_pass,iiflim_fail,iicntr,iioffa]


;;----------------------------------------------------------------------------------------
;; Final output indices
;;----------------------------------------------------------------------------------------
;; AGN with X-ray detections/non-detections
iiagn_det = 'IIAGN_DET'+xfield
iiagn_drm = 'IIAGN_DRM'+xfield
iiagn_non = 'IIAGN_NON'+xfield
iiagn_nrm = 'IIAGN_NRM'+xfield
for f = 0,nfield-1 do begin
    re = execute(iiagn_det[f]+' = IIINF'+xfield[f]+' and IIAGN and IIDET'+xfield[f]+' and IICLEAN'+xfield[f]+' and IISN'+xfield[f]+' and IICNTR'+xfield[f])
    re = execute(iiagn_drm[f]+' = IIINF'+xfield[f]+' and IIAGN and IIDET'+xfield[f]+' and IICLEAN'+xfield[f]+' and (~IISN'+xfield[f]+' or IIOFFA'+xfield[f]+')')
    re = execute(iiagn_non[f]+' = IIINF'+xfield[f]+' and IIAGN and ~IIDET'+xfield[f]+' and IIFLIM_PASS'+xfield[f]+' and IICNTR'+xfield[f])
    re = execute(iiagn_nrm[f]+' = IIINF'+xfield[f]+' and IIAGN and ~IIDET'+xfield[f]+' and (IIFLIM_FAIL'+xfield[f]+' or IIOFFA'+xfield[f]+')')
endfor

sav_vars = [sav_vars]
sav_inds = [sav_inds,iiagn_det,iiagn_drm,iiagn_non,iiagn_nrm]


;; combined raw X-ray detections, sample-catalog matches
re = execute('iidet = '+strjoin('IIDET'+xfield,' or '))
;; combined cleaned X-ray detections, sample-cleaned catalog matches
re = execute('iiclean = '+strjoin('IICLEAN'+xfield,' or '))
;; combined final source detections
re = execute('iiagn_det = '+strjoin('IIAGN_DET'+xfield,' or '))
;; intermediate step to non-detections (MEANINGLESS AS A FULL SAMPLE COMBINED VECTOR)
re = execute('iiagn_drm = '+strjoin('IIAGN_DRM'+xfield,' or '))
re = execute('iiagn_nrm = '+strjoin('IIAGN_NRM'+xfield,' or '))
;; combined final source non-detections (in AGN_NON, no detections at all, not a non-detection that was removed)
re = execute('iiagn_non = ('+strjoin('IIAGN_NON'+xfield,' or ')+') and ~iidet and ~iiagn_nrm')

sav_vars = [sav_vars]
sav_inds = [sav_inds,'IIDET','IICLEAN','IIAGN_DET','IIAGN_NON']


sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="quality_src.sav"')


END







