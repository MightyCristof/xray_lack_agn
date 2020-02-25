PRO compute_luminosity_ratio 


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
common _quality


;;----------------------------------------------------------------------------------------
;; COMBINE X-RAY LUMINOSITIES
;;----------------------------------------------------------------------------------------
lx = dblarr(nsrc)
e_lx = dblarr(nsrc)
loglx = dblarr(nsrc)
e_loglx = dblarr(nsrc)
lxlim = dblarr(nsrc)
e_lxlim = dblarr(nsrc)
loglxlim = dblarr(nsrc)
e_loglxlim = dblarr(nsrc)
for f = 0,nfield-1 do begin
    re = execute('idet_fill = where(lx eq 0. and IIAGN_DET'+xfield[f]+',detct)')
    if (detct gt 0.) then begin
        re = execute('lx[idet_fill] = lx'+xfield[f]+'[idet_fill]')
        re = execute('e_lx[idet_fill] = e_lx'+xfield[f]+'[idet_fill]')
        re = execute('loglx[idet_fill] = loglx'+xfield[f]+'[idet_fill]')
        re = execute('e_loglx[idet_fill] = e_loglx'+xfield[f]+'[idet_fill]')
    endif
    re = execute('inon_fill = where(lxlim eq 0. and IIAGN_NON'+xfield[f]+',nonct)')
    if (nonct gt 0.) then begin
        re = execute('lxlim[inon_fill] = lxlim'+xfield[f]+'[inon_fill]')
        re = execute('e_lxlim[inon_fill] = e_lxlim'+xfield[f]+'[inon_fill]')
        re = execute('loglxlim[inon_fill] = loglxlim'+xfield[f]+'[inon_fill]')
        re = execute('e_loglxlim[inon_fill] = e_loglxlim'+xfield[f]+'[inon_fill]')
    endif
endfor

sav_vars = ['LX','E_LX','LOGLX','E_LOGLX','LXLIM','E_LXLIM','LOGLXLIM','E_LOGLXLIM']
sav_inds = []


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
        ;; start in linear space
        re = execute(lldet[f]+'[where(IIAGN_DET'+xfield[f]+')] = lx'+xfield[f]+'[where(IIAGN_DET'+xfield[f]+')]/(lxir[where(IIAGN_DET'+xfield[f]+')]*lxir_scat[where(IIAGN_DET'+xfield[f]+')])')
        re = execute(e_lldet[f]+'[where(IIAGN_DET'+xfield[f]+')] = '+lldet[f]+'[where(IIAGN_DET'+xfield[f]+')] * sqrt((ERR'+xfield[f]+'[where(IIAGN_DET'+xfield[f]+')]/(FLX'+xfield[f]+'[where(IIAGN_DET'+xfield[f]+')]))^2. + (flx_sigm[1,where(IIAGN_DET'+xfield[f]+')]/flx_sigm[0,where(IIAGN_DET'+xfield[f]+')])^2.)')
        ;; convert to log space
        re = execute(e_lldet[f]+'[where(IIAGN_DET'+xfield[f]+')] /= (alog(10.)*'+lldet[f]+'[where(IIAGN_DET'+xfield[f]+')])')
        re = execute(lldet[f]+'[where(IIAGN_DET'+xfield[f]+')] = alog10('+lldet[f]+'[where(IIAGN_DET'+xfield[f]+')])')
    endif
    ;; non-detections
    re = execute('inon_fill = where(IIAGN_NON'+xfield[f]+',nonct)')
    if (nonct gt 0.) then re = execute(llnon[f]+'[where(IIAGN_NON'+xfield[f]+')] = LOGLXLIM'+xfield[f]+'[where(IIAGN_NON'+xfield[f]+')]-loglxir[where(IIAGN_NON'+xfield[f]+')]')
endfor

sav_vars = [sav_vars,lldet,llnon,e_lldet,e_llnon]
sav_inds = [sav_inds]


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
    if (nonct gt 0.) then re = execute('llxnon[inon_fill] = '+llnon[f]+'[inon_fill]')
endfor
iixdet = llxdet ne -9999.
iixnon = llxnon ne -9999.

sav_vars = [sav_vars,'LLXDET','LLXNON','E_LLXDET','E_LLXNON']
sav_inds = [sav_inds,'IIXDET','IIXNON']


;; log E(B-V) for ploting
lebv = alog10(ebv)
lebv = lebv + randomu(seed,nsrc)*0.05-0.05/2.
lebv[where(~finite(lebv),/null)] = -2.5 

sav_vars = [sav_vars,'LEBV']
sav_inds = [sav_inds]


sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="lum_ratio.sav"')


END

















