PRO source_quality_cuts


common _fits      
common _resamp 
common _comp         
common _inf_cha     
common _inf_xmm     
common _inf_nst     
common _det_cha     
common _det_xmm     
common _det_nst     
common _det_wac     
common _xconv      
common _fxlim     
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

;; passes all SED quality cuts: chi-square, IR bright, and AGN fraction
iiqual = iichi and iiirb and iifagn

sav_vars = ['CHI','DOF','RCHI']
sav_inds = ['IICHI','IIIRB','IIIRW','IIFAGN','IIQUAL']


;;----------------------------------------------------------------------------------------
;; X-ray Catalog Quality Control
;;----------------------------------------------------------------------------------------
sn = 'SN'+xfield
iisn = 'IISN'+xfield
;; S/N cut on X-ray fluxes
for i = 0,nfield-1 do begin
    re = execute(sn[i]+' = FX'+xfield[i]+'/E_FX'+xfield[i])
    re = execute(iisn[i]+' = SN'+xfield[i]+' ge 3.')
endfor

;; sources within 1/2 FOV
;; requires a separate fail flag to account for INF (iioffa =/= ~iicntr)
;; FOV_EFF where effective area ³ 70%
;; NOTE: NO DISTINCTION BETWEEN DETECTION AND NON-DETECTION
iicntr = 'IICNTR'+xfield
iioffa = 'IIOFFA'+xfield
for i = 0,nfield-1 do begin 
    re = execute(iicntr[i]+' = IIINF'+xfield[i]+' and SDST'+xfield[i]+' gt 0. and SDST'+xfield[i]+' le FOV_EFF'+xfield[i])
    re = execute(iioffa[i]+' = IIINF'+xfield[i]+' and SDST'+xfield[i]+' gt FOV_EFF'+xfield[i])
endfor

;; sources above the flux limit
;; requires a separate fail flag to account for FLIM exists (iifail =/= ~iipass)
iifxlim_pass = 'IIFXLIM_PASS'+xfield
iifxlim_fail = 'IIFXLIM_FAIL'+xfield
for i = 0,nfield-1 do begin
    re = execute(iifxlim_pass[i]+' = iilir and FXLIM'+xfield[i]+' gt 0. and fxir ge FXLIM'+xfield[i])
    re = execute(iifxlim_fail[i]+' = iilir and FXLIM'+xfield[i]+' gt 0. and fxir lt FXLIM'+xfield[i])
endfor


sav_vars = [sav_vars,sn]
sav_inds = [sav_inds,iisn,iicntr,iioffa,iifxlim_pass,iifxlim_fail]


;;----------------------------------------------------------------------------------------
;; Final output indices
;;----------------------------------------------------------------------------------------
;; AGN with X-ray detections/non-detections
iifinal_det = 'IIFINAL_DET'+xfield
iifinal_drm = 'IIFINAL_DRM'+xfield
iifinal_non = 'IIFINAL_NON'+xfield
iifinal_nrm = 'IIFINAL_NRM'+xfield
for i = 0,nfield-1 do begin
    re = execute(iifinal_det[i]+' = IIDET'+xfield[i]+' and IIQUAL and IICLEAN'+xfield[i]+' and IISN'+xfield[i]+' and IICNTR'+xfield[i])
    re = execute(iifinal_drm[i]+' = IIDET'+xfield[i]+' and IIQUAL and (IIDIRTY'+xfield[i]+' or ~IISN'+xfield[i]+' or IIOFFA'+xfield[i]+')')
    re = execute(iifinal_non[i]+' = IINON'+xfield[i]+' and IIQUAL and IIFXLIM_PASS'+xfield[i]+' and IICNTR'+xfield[i])
    re = execute(iifinal_nrm[i]+' = IINON'+xfield[i]+' and IIQUAL and (IIFXLIM_FAIL'+xfield[i]+' or IIOFFA'+xfield[i]+')')
endfor

sav_vars = [sav_vars]
sav_inds = [sav_inds,iifinal_det,iifinal_drm,iifinal_non,iifinal_nrm]

;; combined final source detections
re = execute('iifinal_det = '+strjoin('IIFINAL_DET'+xfield,' or '))
;; intermediate step to non-detections (MEANINGLESS AS A FULL SAMPLE COMBINED VECTOR)
re = execute('iifinal_drm = '+strjoin('IIFINAL_DRM'+xfield,' or '))
re = execute('iifinal_nrm = '+strjoin('IIFINAL_NRM'+xfield,' or '))
;; combined final source non-detections (in FIN_NON, no detections at all)
re = execute('iifinal_non = ('+strjoin('IIFINAL_NON'+xfield,' or ')+') and ~iix')
;; combined all final source detections and non-detections
iifinal = iifinal_det or iifinal_non

sav_vars = [sav_vars]
sav_inds = [sav_inds,'IIFINAL_DET','IIFINAL_NON','IIFINAL']


sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="quality_src.sav"')


END







