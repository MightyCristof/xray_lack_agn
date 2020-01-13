PRO compute_luminosity_ratio 


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
common _quality


;;----------------------------------------------------------------------------------------
;; COMBINE X-RAY LUMINOSITIES
;;----------------------------------------------------------------------------------------
lx = dblarr(nsrc)
lx[where(lx eq 0. and iiagn_det_cha)] = lx_cha[where(lx eq 0. and iiagn_det_cha)]
lx[where(lx eq 0. and iiagn_det_xmm)] = lx_xmm[where(lx eq 0. and iiagn_det_xmm)]
lx[where(lx eq 0. and iiagn_det_nst)] = lx_nst[where(lx eq 0. and iiagn_det_nst)]

;;----------------------------------------------------------------------------------------
;; LUMINOSITY RATIOS -- PROXY FOR OBSCURATION
;;----------------------------------------------------------------------------------------
;; detections/non-detections luminosity ratios
lldet = 'LLDET'+xfield
llnon = 'LLNON'+xfield
for f = 0,nfield-1 do begin
    re = execute(lldet[f]+' = dblarr(nsrc)-9999.')
    re = execute(llnon[f]+' = dblarr(nsrc)-9999.')
    re = execute('if (n_elements(where(IIAGN_DET'+xfield[f]+',/null)) gt 0) then '+lldet[f]+'[where(IIAGN_DET'+xfield[f]+')] = LX'+xfield[f]+'[where(IIAGN_DET'+xfield[f]+')]-lxir[where(IIAGN_DET'+xfield[f]+')]') 
    re = execute('if (n_elements(where(IIAGN_NON'+xfield[f]+',/null)) gt 0) then '+llnon[f]+'[where(IIAGN_NON'+xfield[f]+')] = LXLIM'+xfield[f]+'[where(IIAGN_NON'+xfield[f]+')]-lxir[where(IIAGN_NON'+xfield[f]+')]') 
endfor

;; combined all detections to single sample vector 
llxdet = dblarr(nsrc)-9999.
llxnon = dblarr(nsrc)-9999.
;; valid detections/non-detections (supplement with removed detections)
;; non-detections must use combined IIAGN_NON to account for valid non-detections,
;; which have been removed from other instruments
for f = 0,nfield-1 do begin
    re = execute('iidet_fill = llxdet eq -9999. and iiagn_det and IIAGN_DET'+xfield[f])
    re = execute('iinon_fill = llxnon eq -9999. and iiagn_non and IIAGN_NON'+xfield[f])
    re = execute('llxdet[where(iidet_fill,/null)] = '+lldet[f]+'[where(iidet_fill,/null)]')
    re = execute('llxnon[where(iinon_fill,/null)] = '+llnon[f]+'[where(iinon_fill,/null)]')
endfor
iixdet = llxdet ne -9999.
iixnon = llxnon ne -9999.

sav_vars = ['LX',lldet,llnon,'LLXDET','LLXNON']
sav_inds = ['IIXDET','IIXNON']


;; log E(B-V) for ploting
lebv = alog10(ebv)>(-2.5)
lebv = lebv + randomu(seed,n_elements(ebv))*0.05-0.05/2.

sav_vars = [sav_vars,'LEBV']
sav_inds = [sav_inds]


sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="lum_ratio.sav"')


END

















