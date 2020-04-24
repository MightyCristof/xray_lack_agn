PRO source_quality_cuts


common _fits      
common _resamp  
common _inf_cha     
common _inf_xmm     
common _inf_nst     
common _det_cha     
common _det_xmm     
common _det_nst     
common _det_wac     
common _softx       
common _fxlim     
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
;; separate IR bright and weak sources
iiirb = lir ge 42.
iiirw = lir gt 0. and lir lt 42.
;; AGN fraction minimum of 70%
iifagn = agnf15.obs ge 0.7

;; passes all quality cuts: chi-square, IR bright, and AGN fraction
iiqual = iichi and iiirb and iifagn

sav_vars = ['CHI','DOF','RCHI']
sav_inds = ['IICHI','IIIRB','IIIRW','IIFAGN','IIQUAL'];'IIIRC','IIEBV','IIQUAL','IIFAGN']

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
iifxlim_pass = 'IIFXLIM_PASS'+xfield
iifxlim_fail = 'IIFXLIM_FAIL'+xfield
for f = 0,nfield-1 do begin
    re = execute(iifxlim_pass[f]+' = iilir and FXLIM'+xfield[f]+' gt 0. and fxir ge FXLIM'+xfield[f])
    re = execute(iifxlim_fail[f]+' = iilir and FXLIM'+xfield[f]+' gt 0. and fxir lt FXLIM'+xfield[f])
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
sav_inds = [sav_inds,iifxlim_pass,iifxlim_fail,iicntr,iioffa]


;;----------------------------------------------------------------------------------------
;; Final output indices
;;----------------------------------------------------------------------------------------
;; AGN with X-ray detections/non-detections
iiagn_det = 'IIAGN_DET'+xfield
iiagn_drm = 'IIAGN_DRM'+xfield
iiagn_non = 'IIAGN_NON'+xfield
iiagn_nrm = 'IIAGN_NRM'+xfield
for f = 0,nfield-1 do begin
    re = execute(iiagn_det[f]+' = IIINF'+xfield[f]+' and IIQUAL and IIDET'+xfield[f]+' and IICLEAN'+xfield[f]+' and IISN'+xfield[f]+' and IICNTR'+xfield[f])
    re = execute(iiagn_drm[f]+' = IIINF'+xfield[f]+' and IIQUAL and IIDET'+xfield[f]+' and IICLEAN'+xfield[f]+' and (~IISN'+xfield[f]+' or IIOFFA'+xfield[f]+')')
    re = execute(iiagn_non[f]+' = IIINF'+xfield[f]+' and IIQUAL and ~IIDET'+xfield[f]+' and IIFXLIM_PASS'+xfield[f]+' and IICNTR'+xfield[f])
    re = execute(iiagn_nrm[f]+' = IIINF'+xfield[f]+' and IIQUAL and ~IIDET'+xfield[f]+' and (IIFXLIM_FAIL'+xfield[f]+' or IIOFFA'+xfield[f]+')')
endfor

sav_vars = [sav_vars]
sav_inds = [sav_inds,iiagn_det,iiagn_drm,iiagn_non,iiagn_nrm]

;; combined sample-catalog matches
re = execute('iix = '+strjoin('IIX'+xfield,' or '))
;; combined raw X-ray detections
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
;; combined all final source detections and non-detections
iiagn = iiagn_det or iiagn_non

sav_vars = [sav_vars]
sav_inds = [sav_inds,'IIX','IIDET','IICLEAN','IIAGN_DET','IIAGN_NON','IIAGN']


sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="quality_src.sav"')


END







