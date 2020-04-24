PRO agn_xray_luminosity, DERED = dered, $
                         PLT = plt


common _fits
common _resamp
common _inf_cha
common _inf_xmm
common _inf_nst
common _det_nst
common _det_xmm
common _det_cha
common _softx
common _fxlim
common _comp


;; fit output 
ebv = reform(param[0,*])
e_ebv = reform(ebv_sigm[1,*])
c_agn = reform(param[2,*])
iilir = c_agn gt 0.                 ;; AGN SED contribution exists
ilir = where(iilir,nagn,ncomplement=ngal,/null)

;;----------------------------------------------------------------------------------------
;; AGN FRACTION -- 
;;----------------------------------------------------------------------------------------
;; calculate AGN fraction at 6- and 15-micron
agnf6 = f_agn(6.,param,model=agnm6)
agnf15 = f_agn(15.,param,model=agnm15)

sav_vars = ['EBV','E_EBV','NAGN','NGAL','AGNF6','AGNM6','AGNF15','AGNM15']
sav_inds = ['IILIR']


;;----------------------------------------------------------------------------------------
;; LUMINOSITIES -- L(IR) AND LX(LIR)
;;----------------------------------------------------------------------------------------
;; IR 6-micron AGN luminosity calculated from SED model parameters
lir = dblarr(nsrc)
e_lir = dblarr(nsrc)
loglir = dblarr(nsrc)
e_loglir = dblarr(nsrc)
if keyword_set(dered) then begin
    lir[ilir] = l_agn(6.,dblarr(nagn),z[ilir],c_agn[ilir])    ;; 6-micron intrinsic
endif else $
    lir[ilir] = l_agn(6.,ebv[ilir],z[ilir],c_agn[ilir])    ;; 6-micron observed

;; correct AGN luminosity where template over- or underestimates beyond uncertainties
;; luminosity corrected
lcorr = correct_agn_lum(wave,flux,e_flux,param)
iicorr = lcorr ne 0.
lir[ilir] += lcorr[ilir]
e_lir[ilir] = reform(lir_sigm[1,ilir])
loglir[ilir] = alog10(lir[ilir])
e_loglir[ilir] = e_lir[ilir]/(alog(10.)*lir[ilir])

sav_vars = [sav_vars,'LCORR','LIR','E_LIR','LOGLIR','E_LOGLIR']
sav_inds = [sav_inds,'IICORR']


;; LX(LIR) from LX-LIR relationship
lxir = dblarr(nsrc)
lxir_scat = dblarr(nsrc)-9999.
loglxir = dblarr(nsrc)
loglxir_scat = dblarr(nsrc)-9999.
;lxir[where(iilir)] = lxir_chen(lir[where(iilir)],/scatter)
lxir_rel = lxir_fiore(loglir[ilir],/scatter)
loglxir[ilir] = lxir_rel[0,*]
loglxir_scat[ilir] = lxir_rel[1,*]
lxir[ilir] = 10.^loglxir[ilir]
lxir_scat[ilir] = 10.^loglxir_scat[ilir]

sav_vars = [sav_vars,'LXIR','LXIR_SCAT','LOGLXIR','LOGLXIR_SCAT']
sav_inds = [sav_inds]


;; convert LX(LIR) to FX( LX(LIR) ) for flux-limit plots
;; erg/s to erg/s/cm^2
;; luminosity distance
dl2 = dlum(z,/sq)
fxir = dblarr(nsrc)
fxir_scat = dblarr(nsrc)
logfxir = dblarr(nsrc)
logfxir_scat = dblarr(nsrc)
fxir[ilir] = lxir[ilir]/(4.*!const.pi*dl2[ilir])
fxir_scat[ilir] = lxir_scat[ilir]
logfxir[ilir] = alog10(fxir[ilir])
logfxir_scat[ilir] = loglxir_scat[ilir]

sav_vars = [sav_vars,'DL2','FXIR','FXIR_SCAT','LOGFXIR','LOGFXIR_SCAT']
sav_inds = [sav_inds]



;;----------------------------------------------------------------------------------------
;; LUMINOSITIES -- LX
;;----------------------------------------------------------------------------------------
;; 2-10keV luminosity arrays
lx = 'LX'+xfield
e_lx = 'E_'+lx
loglx = 'LOG'+lx
e_loglx = 'E_'+loglx
;; X-ray luminosity converted from flux, K-corrected to rest-frame Fkcor = Fobs*(1+z)^(Γ-2)
;; Chandra, XMM, NuSTAR
cat_gamma = [1.8,1.8,1.8]
for i = 0,nfield-1 do begin
    re = execute(lx[i]+' = dblarr(nsrc)')
    re = execute(e_lx[i]+' = dblarr(nsrc)')
    re = execute(loglx[i]+' = dblarr(nsrc)-9999.')
    re = execute(e_loglx[i]+' = dblarr(nsrc)')
    re = execute('idet = where(IIDET'+xfield[i]+' and FLX'+xfield[i]+' gt 0. and ERR'+xfield[i]+' gt 0.,detct)')
    if (detct gt 0.) then begin
        re = execute('flx_kcorr = FLX'+xfield[i]+'[idet]*(1+z[idet])^(cat_gamma[i]-2.)')
        re = execute(lx[i]+'[idet] = 4.*!const.pi*dl2[idet]*flx_kcorr')
        re = execute(e_lx[i]+'[idet] = '+lx[i]+'[idet] * sqrt((ERR'+xfield[i]+'[idet]/flx_kcorr)^2. + (red_sigm[1,idet]/red_sigm[0,idet])^2.)')
        re = execute(loglx[i]+'[idet] = alog10('+lx[i]+'[idet])')
        re = execute(e_loglx[i]+'[idet] = '+e_lx[i]+'[idet]/(alog(10.)*'+lx[i]+'[idet])')
    endif
endfor

sav_vars = [sav_vars,'CAT_GAMMA',lx,e_lx,loglx,e_loglx]
sav_inds = [sav_inds]


;;----------------------------------------------------------------------------------------
;; INTERPOLATE FLUX LIMIT FLUX & LUMINOSITY
;;----------------------------------------------------------------------------------------
;; 2-10keV
fxlim = 'FXLIM'+xfield        ;; flux limit at source
e_fxlim = 'E_'+fxlim
logfxlim = 'LOG'+fxlim
e_logfxlim = 'E_'+logfxlim
lxlim = 'LXLIM'+xfield      ;; luminosity at flux limit
e_lxlim = 'E_'+lxlim
loglxlim = 'LOG'+lxlim
e_loglxlim = 'E_'+loglxlim
fxlim_cs = 'FXLIM_CS'+xfield  ;; flux-limit function coefficients 
degr = [6,6,1]              ;; degree of polynomial to fit flux-limit
for i = 0,nfield-1 do begin
    re = execute(fxlim[i]+' = dblarr(nsrc)')
    re = execute(e_fxlim[i]+' = dblarr(nsrc)')
    re = execute(logfxlim[i]+' = dblarr(nsrc)-9999.')
    re = execute(e_logfxlim[i]+' = dblarr(nsrc)')
    re = execute(lxlim[i]+' = dblarr(nsrc)')
    re = execute(e_lxlim[i]+' = dblarr(nsrc)')
    re = execute(loglxlim[i]+' = dblarr(nsrc)-9999.')
    re = execute(e_loglxlim[i]+' = dblarr(nsrc)')
    re = execute(fxlim_cs[i]+' = dblarr(degr[i])')
    re = execute('isrc = where(iiinf'+xfield[i]+')')
    re = execute(fxlim[i]+'[isrc] = extrapolate_flim(CAT_LIM'+xfield[i]+',CAT_EXP'+xfield[i]+',TEXP'+xfield[i]+'[isrc],degr[i],FLIM_CS='+fxlim_cs[i]+')')
    re = execute(e_fxlim[i]+'[isrc] = '+fxlim[i]+'[isrc] / 3d')
    re = execute(logfxlim[i]+'[isrc] = alog10('+fxlim[i]+'[isrc])')
    re = execute(e_logfxlim[i]+'[isrc] = '+e_fxlim[i]+'[isrc]/(alog(10)*'+fxlim[i]+'[isrc])')
    re = execute(lxlim[i]+'[isrc] = 4.*!const.pi*dl2[isrc]*'+fxlim[i]+'[isrc]')
    re = execute(e_lxlim[i]+'[isrc] = '+lxlim[i]+'[isrc] * sqrt(('+e_fxlim[i]+'[isrc]/'+fxlim[i]+'[isrc])^2. + (red_sigm[1,isrc]/red_sigm[0,isrc])^2.)')
    re = execute(loglxlim[i]+'[isrc] = alog10('+lxlim[i]+'[isrc])')
    re = execute(e_loglxlim[i]+'[isrc] = '+e_lxlim[i]+'[isrc]/(alog(10)*'+lxlim[i]+'[isrc])')
endfor

sav_vars = [sav_vars,'DEGR',fxlim_cs,fxlim,e_fxlim,logfxlim,e_logfxlim,lxlim,e_lxlim,loglxlim,e_loglxlim]
sav_inds = [sav_inds]


sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="src_luminosity.sav"')





END


