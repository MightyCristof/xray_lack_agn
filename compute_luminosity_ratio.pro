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
common _xconv      
common _fxlim    
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
llnon = 'LLNON'+xfield
e_lldet = 'E_'+lldet
e_llnon = 'E_'+llnon
for i = 0,nfield-1 do begin
    re = execute(lldet[i]+' = dblarr(nsrc)-9999.')
    re = execute(llnon[i]+' = dblarr(nsrc)-9999.')
    re = execute(e_lldet[i]+' = dblarr(nsrc)-9999.')
    re = execute(e_llnon[i]+' = dblarr(nsrc)-9999.')
    ;; detections
    re = execute('iivalid = IIFINAL_DET'+xfield[i])
    ivalid = where(iivalid,detct)
    if (detct gt 0.) then begin
        ;; start in linear space
        re = execute(lldet[i]+'[ivalid] = lx'+xfield[i]+'[ivalid]/(lxir[ivalid]*lxir_scat[ivalid])')
        ;; errors attributed to X-ray flux and IR flux
        re = execute(e_lldet[i]+'[ivalid] = '+lldet[i]+'[ivalid] * sqrt((ERR'+xfield[i]+'[ivalid]/FLX'+xfield[i]+'[ivalid])^2. + (flx_sigm[1,ivalid]/flx_sigm[0,ivalid])^2.)')
        ;; check resamp distribution for 0 (MEDABSDEV==0), -1 (only one source), or -9999 (should not ever be the case where AGN component; sanity check)
        iest = where(iivalid and flx_sigm[1,*] le 0.,estct)
        if (estct gt 0) then stop
        ;; convert to log space
        re = execute(e_lldet[i]+'[ivalid] /= (alog(10.)*'+lldet[i]+'[ivalid])>(-9999.)')
        re = execute(lldet[i]+'[ivalid] = alog10('+lldet[i]+'[ivalid])>(-9999.)')
    endif
    ;; non-detections
    re = execute('iivalid = IIFINAL_NON'+xfield[i])
    ivalid = where(iivalid,nonct)
    if (nonct gt 0.) then re = execute(llnon[i]+'[ivalid] = LOGLXLIM'+xfield[i]+'[ivalid]-loglxir[ivalid]')
endfor

sav_vars = [lldet,e_lldet,llnon,e_llnon]
sav_inds = []


;;----------------------------------------------------------------------------------------
;; COMBINE X-RAY LUMINOSITIES AND RATIOS
;;----------------------------------------------------------------------------------------
;; luminosities
lx = dblarr(nsrc)
e_lx = dblarr(nsrc)
loglx = dblarr(nsrc)-9999.
e_loglx = dblarr(nsrc)-9999.
;; non-detections limits
lxlim = dblarr(nsrc)
e_lxlim = dblarr(nsrc)
loglxlim = dblarr(nsrc)-9999.
e_loglxlim = dblarr(nsrc)-9999.
;; ratios
lldet = dblarr(nsrc)-9999.
e_lldet = dblarr(nsrc)-9999.
llnon = dblarr(nsrc)-9999.
e_llnon = dblarr(nsrc)-9999.
for i = 0,nfield-1 do begin
    re = execute('idet_fill = where(lx eq 0. and IIFINAL_DET'+xfield[i]+',detct)')
    if (detct gt 0.) then begin
        re = execute('lx[idet_fill] = lx'+xfield[i]+'[idet_fill]')
        re = execute('e_lx[idet_fill] = e_lx'+xfield[i]+'[idet_fill]')
        re = execute('loglx[idet_fill] = loglx'+xfield[i]+'[idet_fill]')
        re = execute('e_loglx[idet_fill] = e_loglx'+xfield[i]+'[idet_fill]')
        re = execute('lldet[idet_fill] = LLDET'+xfield[i]+'[idet_fill]')
        re = execute('e_lldet[idet_fill] = E_LLDET'+xfield[i]+'[idet_fill]')
    endif
    re = execute('inon_fill = where(lxlim eq 0. and iifinal_non and IIFINAL_NON'+xfield[i]+',nonct)')
    if (nonct gt 0.) then begin
        re = execute('lxlim[inon_fill] = lxlim'+xfield[i]+'[inon_fill]')
        re = execute('e_lxlim[inon_fill] = e_lxlim'+xfield[i]+'[inon_fill]')
        re = execute('loglxlim[inon_fill] = loglxlim'+xfield[i]+'[inon_fill]')
        re = execute('e_loglxlim[inon_fill] = e_loglxlim'+xfield[i]+'[inon_fill]')
        re = execute('llnon[inon_fill] = LLNON'+xfield[i]+'[inon_fill]')
        re = execute('e_llnon[inon_fill] = E_LLNON'+xfield[i]+'[inon_fill]')
    endif
endfor

sav_vars = [sav_vars,'LX','E_LX','LOGLX','E_LOGLX', $
                     'LXLIM','E_LXLIM','LOGLXLIM','E_LOGLXLIM', $
                     'LLDET','E_LLDET','LLNON','E_LLNON']
sav_inds = [sav_inds]


;;----------------------------------------------------------------------------------------
;; LOG E(B-V)
;;----------------------------------------------------------------------------------------
lebv = alog10(ebv)
lebv = lebv + randomu(seed,nsrc)*0.05-0.05/2.
lebv[where(~finite(lebv),/null)] = -2.5 

sav_vars = [sav_vars,'LEBV']
sav_inds = [sav_inds]


sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="lum_ratio.sav"')


END

















