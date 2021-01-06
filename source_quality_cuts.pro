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
common _wac     
common _xconv      
common _fxlim     
common _agnlum



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
iiirb = loglir ge 42.
iiirw = loglir gt 0. and loglir lt 42.
;; AGN fraction minimum of 70%
;iifagn06 = agnf06.obs ge 0.7
;iifagn15 = agnf15.obs ge 0.7
iifagn = agnf15.obs ge 0.7

;; passes all SED quality cuts: chi-square, IR bright, and AGN fraction
iised = iichi and iiirb and iifagn

sav_vars = ['CHI','DOF','RCHI']
sav_inds = ['IICHI','IIIRB','IIIRW','IIFAGN','IISED']


;;----------------------------------------------------------------------------------------
;; X-ray Catalog Quality Control
;;----------------------------------------------------------------------------------------
;; sources within 1/2 FOV
;; requires a separate fail flag to account for INF (iioffa =/= ~iicntr)
;; FOV_EFF where effective area ³ 70%
;; NOTE: NO DISTINCTION BETWEEN DETECTION AND NON-DETECTION
;; IMPORTANT: you MUST keep SDST > 0 to filter out detections outside of the master footprint (e.g., ACIS-S)
iicntr = 'IICNTR'+xfield
iioffa = 'IIOFFA'+xfield
for i = 0,nfield-1 do begin 
    re = execute(iicntr[i]+' = IIINF'+xfield[i]+' and SDST'+xfield[i]+' le FOV_EFF'+xfield[i]+' and SDST'+xfield[i]+' ge 0.')
    re = execute(iioffa[i]+' = IIINF'+xfield[i]+' and SDST'+xfield[i]+' gt FOV_EFF'+xfield[i])
endfor

;; X-ray S/N
sn = 'SN'+xfield
iisn = 'IISN'+xfield
;; S/N cut on X-ray fluxes
for i = 0,nfield-1 do begin
    re = execute(sn[i]+' = FX'+xfield[i]+'/E_FX'+xfield[i])
    re = execute(iisn[i]+' = SN'+xfield[i]+' ge 3.')
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
;; Quality output indices for analysis
;;----------------------------------------------------------------------------------------
;; AGN with X-ray detections/non-detections
iiqual_det = 'IIQUAL_DET'+xfield
iiqual_drm = 'IIQUAL_DRM'+xfield
iiqual_non = 'IIQUAL_NON'+xfield
iiqual_nrm = 'IIQUAL_NRM'+xfield
for i = 0,nfield-1 do begin
    re = execute(iiqual_det[i]+' = IIDET'+xfield[i]+' and IISED and IICNTR'+xfield[i]+' and IISN'+xfield[i])
    re = execute(iiqual_drm[i]+' = IIDET'+xfield[i]+' and IISED and (IIOFFA'+xfield[i]+' or ~IISN'+xfield[i]+')')
;    re = execute(iiqual_det[i]+' = IIDET'+xfield[i]+' and IISED and IICLEAN'+xfield[i]+' and IISN'+xfield[i]+' and IICNTR'+xfield[i])
;    re = execute(iiqual_drm[i]+' = IIDET'+xfield[i]+' and IISED and (IIDIRTY'+xfield[i]+' or ~IISN'+xfield[i]+' or IIOFFA'+xfield[i]+')')
    re = execute(iiqual_non[i]+' = IINON'+xfield[i]+' and IISED and IICNTR'+xfield[i]+' and IIFXLIM_PASS'+xfield[i])
    re = execute(iiqual_nrm[i]+' = IINON'+xfield[i]+' and IISED and (IIOFFA'+xfield[i]+' or IIFXLIM_FAIL'+xfield[i]+')')
endfor

sav_vars = [sav_vars]
sav_inds = [sav_inds,iiqual_det,iiqual_drm,iiqual_non,iiqual_nrm]

;; combined quality source detections
re = execute('iiqual_det = '+strjoin('IIQUAL_DET'+xfield,' or '))
;; intermediate step to non-detections (MEANINGLESS AS A FULL SAMPLE COMBINED VECTOR)
re = execute('iiqual_drm = '+strjoin('IIQUAL_DRM'+xfield,' or '))
re = execute('iiqual_nrm = '+strjoin('IIQUAL_NRM'+xfield,' or '))
;; combined quality source non-detections (in FIN_NON, no detections at all)
re = execute('iiqual_non = ('+strjoin('IIQUAL_NON'+xfield,' or ')+') and ~iix')
;; combined all quality source detections and non-detections
iiqual = iiqual_det or iiqual_non

sav_vars = [sav_vars]
sav_inds = [sav_inds,'IIQUAL_DET','IIQUAL_NON','IIQUAL']


sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="quality_src.sav"')


END







