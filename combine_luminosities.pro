PRO combine_luminosities


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
;; COMBINE X-RAY LUMINOSITIES, LIMITS, AND RATIOS
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

sav_vars = ['LX','E_LX','LOGLX','E_LOGLX', $
            'LXLIM','E_LXLIM','LOGLXLIM','E_LOGLXLIM', $
            'LLDET','E_LLDET','LLNON','E_LLNON']
sav_inds = []


;;----------------------------------------------------------------------------------------
;; LOG E(B-V)
;;----------------------------------------------------------------------------------------
lebv = alog10(ebv)
lebv = lebv + randomu(seed,nsrc)*0.05-0.05/2.
lebv[where(~finite(lebv),/null)] = -2.5 

sav_vars = [sav_vars,'LEBV']
sav_inds = [sav_inds]


sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="combined_lum.sav"')


END

















