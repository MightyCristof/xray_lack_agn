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
e_lx = dblarr(nsrc)
for f = 0,nfield-1 do begin
    re = execute('ilx = where(lx eq 0. and IIAGN_DET'+xfield[f]+',lxct)')
    if (lxct gt 0.) then begin
        re = execute('lx[ilx] = lx'+xfield[f]+'[ilx]')
        re = execute('e_lx[ilx] = ERR'+xfield[f]+'[ilx]/(alog(10.)*FLX'+xfield[f]+'[ilx])')
    endif
endfor

;;----------------------------------------------------------------------------------------
;; LUMINOSITY RATIOS -- PROXY FOR OBSCURATION
;;----------------------------------------------------------------------------------------
;; detections/non-detections luminosity ratios
lldet = 'LLDET'+xfield
llnon = 'LLNON'+xfield
e_lldet = 'E_'+lldet
e_llnon = 'E_'+llnon
for f = 0,nfield-1 do begin
    re = execute(lldet[f]+' = dblarr(nsrc)-9999.')
    re = execute(llnon[f]+' = dblarr(nsrc)-9999.')
    re = execute(e_lldet[f]+' = dblarr(nsrc)-9999.')
    re = execute(e_llnon[f]+' = dblarr(nsrc)-9999.')
    ;; detections
    re = execute('idet_fill = where(IIAGN_DET'+xfield[f]+',detct)')
    if (detct gt 0.) then begin
        re = execute(lldet[f]+'[where(IIAGN_DET'+xfield[f]+')] = LX'+xfield[f]+'[where(IIAGN_DET'+xfield[f]+')]-lxir[where(IIAGN_DET'+xfield[f]+')]')
        re = execute(e_lldet[f]+'[where(IIAGN_DET'+xfield[f]+')] = sqrt((ERR'+xfield[f]+'[where(IIAGN_DET'+xfield[f]+')]/(alog(10.)*FLX'+xfield[f]+'[where(IIAGN_DET'+xfield[f]+')]))^2. + (e_fir[where(IIAGN_DET'+xfield[f]+')])^2.)')
    endif
    ;; non-detections
    re = execute('inon_fill = where(IIAGN_NON'+xfield[f]+',nonct)')
    if (nonct gt 0.) then begin
        re = execute(llnon[f]+'[where(IIAGN_NON'+xfield[f]+')] = LXLIM'+xfield[f]+'[where(IIAGN_NON'+xfield[f]+')]-lxir[where(IIAGN_NON'+xfield[f]+')]')
        re = execute(e_llnon[f]+'[where(IIAGN_NON'+xfield[f]+')] = sqrt((median(CAT_ERR'+xfield[f]+'/(alog(10.)*CAT_FLX'+xfield[f]+')))^2. + (e_fir[where(IIAGN_NON'+xfield[f]+')])^2.)')
    endif
endfor

;; combined all detections to single sample vector 
llxdet = dblarr(nsrc)-9999.
llxnon = dblarr(nsrc)-9999.
e_llxdet = dblarr(nsrc)-9999.
e_llxnon = dblarr(nsrc)-9999.
;; valid detections/non-detections (supplement with removed detections)
;; non-detections must use combined IIAGN_NON to account for valid non-detections,
;; which have been removed from other instruments
for f = 0,nfield-1 do begin
    ;; detections
    re = execute('idet_fill = where(llxdet eq -9999. and iiagn_det and IIAGN_DET'+xfield[f]+',detct)')
    if (detct gt 0.) then begin
        re = execute('llxdet[idet_fill] = '+lldet[f]+'[idet_fill]')
        re = execute('e_llxdet[idet_fill] = '+e_lldet[f]+'[idet_fill]')
    endif
    ;; non-detections
    re = execute('inon_fill = where(llxnon eq -9999. and iiagn_non and IIAGN_NON'+xfield[f]+',nonct)')
    if (nonct gt 0.) then begin
        re = execute('llxnon[inon_fill] = '+llnon[f]+'[inon_fill]')
        re = execute('e_llxnon[inon_fill] = '+e_llnon[f]+'[inon_fill]')        
    endif
endfor
iixdet = llxdet ne -9999.
iixnon = llxnon ne -9999.

sav_vars = ['LX','E_LX',lldet,llnon,e_lldet,e_llnon, $
            'LLXDET','LLXNON','E_LLXDET','E_LLXNON']
sav_inds = ['IIXDET','IIXNON']


;; log E(B-V) for ploting
lebv = alog10(ebv)>(-2.5)
lebv = lebv + randomu(seed,nsrc)*0.05-0.05/2.

sav_vars = [sav_vars,'LEBV']
sav_inds = [sav_inds]


sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="lum_ratio.sav"')


END

















