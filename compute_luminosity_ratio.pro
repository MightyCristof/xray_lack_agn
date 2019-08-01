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
;; LUMINOSITY RATIOS -- PROXY FOR OBSCURATION
;;----------------------------------------------------------------------------------------
;; detections/non-detections luminosity ratios
lldet = 'LLDET'+xfield
lldrm = 'LLDRM'+xfield
llnon = 'LLNON'+xfield
llnrm = 'LLNRM'+xfield
for f = 0,nfield-1 do begin
    re = execute(lldet[f]+' = dblarr(nsrc)-9999.')
    re = execute(lldrm[f]+' = dblarr(nsrc)-9999.')
    re = execute(llnon[f]+' = dblarr(nsrc)-9999.')
    re = execute(llnrm[f]+' = dblarr(nsrc)-9999.')
    re = execute('if (n_elements(where(IIAGN_DET'+xfield[f]+',/null)) gt 0) then '+lldet[f]+'[where(IIAGN_DET'+xfield[f]+')] = LX'+xfield[f]+'[where(IIAGN_DET'+xfield[f]+')]-lxir[where(IIAGN_DET'+xfield[f]+')]') 
    re = execute('if (n_elements(where(IIAGN_DRM'+xfield[f]+',/null)) gt 0) then '+lldrm[f]+'[where(IIAGN_DRM'+xfield[f]+')] = LX'+xfield[f]+'[where(IIAGN_DRM'+xfield[f]+')]-lxir[where(IIAGN_DRM'+xfield[f]+')]') 
    re = execute('if (n_elements(where(IIAGN_NON'+xfield[f]+',/null)) gt 0) then '+llnon[f]+'[where(IIAGN_NON'+xfield[f]+')] = LXLIM'+xfield[f]+'[where(IIAGN_NON'+xfield[f]+')]-lxir[where(IIAGN_NON'+xfield[f]+')]') 
    re = execute('if (n_elements(where(IIAGN_NRM'+xfield[f]+',/null)) gt 0) then '+llnrm[f]+'[where(IIAGN_NRM'+xfield[f]+')] = LXLIM'+xfield[f]+'[where(IIAGN_NRM'+xfield[f]+')]-lxir[where(IIAGN_NRM'+xfield[f]+')]') 
endfor

;; combined all detections to single sample vector 
llxdet = dblarr(nsrc)-9999.
llxnon = dblarr(nsrc)-9999.
;; valid detections/non-detections (supplement with removed detections)
for f = 0,nfield-1 do begin
    iidfill = llxdet eq -9999.
    iinfill = llxnon eq -9999.
    re = execute('llxdet[where(IIAGN_DET'+xfield[f]+' or IIAGN_DRM'+xfield[f]+' and iidfill,/null)] = '+lldet[f]+'[where(IIAGN_DET'+xfield[f]+' or IIAGN_DRM'+xfield[f]+' and iidfill,/null)]')
    re = execute('llxnon[where(IIAGN_NON'+xfield[f]+' and iinfill,/null)] = '+llnon[f]+'[where(IIAGN_NON'+xfield[f]+' and iinfill,/null)]')
endfor
iixdet = llxdet ne -9999.
iixnon = llxnon ne -9999.


sav_vars = [lldet,lldrm,llnon,llnrm,'LLXDET','LLXNON']
sav_inds = ['IIXDET','IIXNON']

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="lum_ratio.sav"')


END

















