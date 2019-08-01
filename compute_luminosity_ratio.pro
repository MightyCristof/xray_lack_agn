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
common _agn_lum    
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
for i = 0,nfield-1 do begin
    re = execute(lldet[i]+' = dblarr(nsrc)-9999.')
    re = execute(lldrm[i]+' = dblarr(nsrc)-9999.')
    re = execute(llnon[i]+' = dblarr(nsrc)-9999.')
    re = execute(llnrm[i]+' = dblarr(nsrc)-9999.')
    re = execute('if (n_elements(where(IIAGN_DET'+xfield[i]+',/null)) gt 0) then '+lldet[i]+'[where(IIAGN_DET'+xfield[i]+')] = '+lx[i]+'[where(IIAGN_DET'+xfield[i]+')]-lxir[where(IIAGN_DET'+xfield[i]+')]') 
    re = execute('if (n_elements(where(IIAGN_DRM'+xfield[i]+',/null)) gt 0) then '+lldrm[i]+'[where(IIAGN_DRM'+xfield[i]+')] = '+lx[i]+'[where(IIAGN_DRM'+xfield[i]+')]-lxir[where(IIAGN_DRM'+xfield[i]+')]') 
    re = execute('if (n_elements(where(IIAGN_NON'+xfield[i]+',/null)) gt 0) then '+llnon[i]+'[where(IIAGN_NON'+xfield[i]+')] = '+lxlim[i]+'[where(IIAGN_NON'+xfield[i]+')]-lxir[where(IIAGN_NON'+xfield[i]+')]') 
    re = execute('if (n_elements(where(IIAGN_NRM'+xfield[i]+',/null)) gt 0) then '+llnrm[i]+'[where(IIAGN_NRM'+xfield[i]+')] = '+lxlim[i]+'[where(IIAGN_NRM'+xfield[i]+')]-lxir[where(IIAGN_NRM'+xfield[i]+')]') 
endfor

;; combined all detections to single sample vector 
llxdet = dblarr(nsrc)-9999.
llxnon = dblarr(nsrc)-9999.
;; valid detections first
for i = 0,nfield-1 do begin
    re = execute('llxdet[where(IIAGN_DET'+xfield[i]+' and llxdet eq -9999.,/null)] = '+lldet[i]+'[where(IIAGN_DET'+xfield[i]+' and llxdet eq -9999.,/null)]')
    re = execute('llxnon[where(IIAGN_NON'+xfield[i]+' and llxnon eq -9999.,/null)] = '+llnon[i]+'[where(IIAGN_NON'+xfield[i]+' and llxnon eq -9999.,/null)]')
endfor
iixdet = llxdet ne -9999.
iixnon = llxnon ne -9999.
;; supplement with removed detections
for i = 0,nfield-1 do begin
    re = execute('llxdet[where(IIAGN_DRM'+xfield[i]+' and llxdet eq -9999.,/null)] = '+lldrm[i]+'[where(IIAGN_DRM'+xfield[i]+' and llxdet eq -9999.,/null)]')
    ;re = execute('llxnrm[where(IIAGN_NRM'+xfield[i]+' and llxnon eq -9999.,/null)] = '+llnrm[i]+'[where(IIAGN_NRM'+xfield[i]+' and llxnon eq -9999.,/null)]')
endfor
iixd = llxdet ne -9999.
iixn = llxnon ne -9999.

nhxdet = ll2nh(llxdet,'2-10')
nhxnon = ll2nh(llxnon,'2-10')

sav_vars = ['LLDET','LLDRM','LLNON','LLNRM', $
            lldet,lldrm,llnon,llnrm, $
            'LLXDET','LLXNON', $
            'NHXDET','NHXNON']
sav_inds = ['IIXDET','IIXNON','IIXD','IIXN']

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="lum_ratio.sav"')


END

















