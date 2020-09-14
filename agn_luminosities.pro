PRO agn_luminosities, DERED = dered, $
                      REL = rel           


common _fits
common _resamp
common _comp
common _inf_cha
common _inf_xmm
common _inf_nst
common _det_nst
common _det_xmm
common _det_cha
common _xconv
common _fxlim


;; check LX-LIR RELATION
rel = strupcase(rel)
if (total([strmatch(rel,'C17'),strmatch(rel,'F09')]) eq 0) then begin
    print, "INVALID LX-LIR RELATION"
    return
endif

;; SED model output
;; where AGN component
c_agn = reform(param[2,*])
iilir = c_agn gt 0.                 ;; AGN SED contribution exists
ilir = where(iilir,nagn,ncomplement=ngal,/null)
;; extinction parameter
ebv = reform(param[0,*])
e_ebv = reform(ebv_sigm[1,*])
;; check resamp distribution for 0 (MEDABSDEV==0), -1 (only one source), -9999 (no AGN)
icest = where(iilir and ebv_sigm[1,*] le 0,cestct)
if (cestct gt 0.) then e_ebv[icest] = ebv[icest]*0.1
;; SED best fit redshift parameter
zsed = reform(param[1,*])
e_zsed = reform(red_sigm[1,*])
;; check resamp distribution for -1 (only one source), use original redshift error
izest = where(e_zsed le 0,zestct)
if (zestct gt 0.) then e_zsed[izest] = zerr[izest]
;; luminosity distance
dl2 = dlum(z,/sq)

;;----------------------------------------------------------------------------------------
;; AGN FRACTION -- 
;;----------------------------------------------------------------------------------------
;; calculate AGN fraction at 6- and 15-micron
agnf6 = f_agn(6.,param,model=agnm6)
agnf15 = f_agn(15.,param,model=agnm15)

sav_vars = ['EBV','E_EBV','ZSED','E_ZSED','NAGN','NGAL','AGNF6','AGNM6','AGNF15','AGNM15','DL2']
sav_inds = ['IILIR']


;;----------------------------------------------------------------------------------------
;; L(IR) AND LX(LIR)
;;----------------------------------------------------------------------------------------

;; IR 6-micron AGN luminosity and flux calculated from SED model parameters
lir = dblarr(nsrc)
e_lir = dblarr(nsrc)
loglir = dblarr(nsrc)-9999.
e_loglir = dblarr(nsrc)-9999.
fir = dblarr(nsrc)
e_fir = dblarr(nsrc)
logfir = dblarr(nsrc)
e_logfir = dblarr(nsrc)
if keyword_set(dered) then begin
    lir[ilir] = l_agn(6.,dblarr(nagn),z[ilir],c_agn[ilir])    ;; 6-micron intrinsic
endif else $
    lir[ilir] = l_agn(6.,ebv[ilir],z[ilir],c_agn[ilir])    ;; 6-micron observed

;; correct AGN luminosity where template over- or underestimates beyond uncertainties
;; luminosity corrected
lcorr = correct_agn_lum(wave,flux,e_flux,param)
iicorr = fix(lcorr ne 0.)
iicorr[where(iicorr and lcorr lt 1.,/null)] = -1
icorr = where(iicorr ne 0.,corrct)
if (corrct gt 0) then lir[icorr] *= lcorr[icorr]
e_lir[ilir] = reform(lir_sigm[1,ilir])
fir[ilir] = lir[ilir]/(4.*!const.pi*dl2[ilir])
e_fir[ilir] = reform(fir_sigm[1,ilir])
;; check resamp distribution for 0 (MEDABSDEV==0), -1 (only one source), -9999 (no AGN)
iest = where(iilir and lir_sigm[1,*] le 0,estct)
if (estct gt 0) then begin 
    e_lir[iest] = lir[iest]*0.1
    e_fir[iest] = fir[iest]*0.1
endif
loglir[ilir] = alog10(lir[ilir])>(-9999.)
e_loglir[ilir] = e_lir[ilir]/(alog(10.)*lir[ilir])>(-9999.)
logfir[ilir] = alog10(fir[ilir])>(-9999.)
e_logfir[ilir] = e_fir[ilir]/(alog(10.)*fir[ilir])>(-9999.)


sav_vars = [sav_vars,'LCORR','LIR','E_LIR','LOGLIR','E_LOGLIR', $
                             'FIR','E_FIR','LOGFIR','E_LOGFIR']
sav_inds = [sav_inds,'IICORR']


;; LX(LIR) from LX-LIR relationship
;; LXIR_SCAT is not uncertainty, it is added scatter about LX(LIR)
;; REMEMBER!: add LOGLXIR_SCAT, multiply LXIR_SCAT
lxir = dblarr(nsrc)
lxir_scat = dblarr(nsrc)
loglxir = dblarr(nsrc)-9999.
loglxir_scat = dblarr(nsrc)-9999.
case rel of 
    'F09': loglxir[ilir] = lxir_f09(loglir[ilir],scatter=scat)
    'C17': loglxir[ilir] = lxir_c17(loglir[ilir],scatter=scat)
endcase
                         
loglxir_scat[ilir] = scat
lxir[ilir] = 10.^loglxir[ilir]
lxir_scat[ilir] = 10.^loglxir_scat[ilir]

sav_vars = [sav_vars,'LXIR','LXIR_SCAT','LOGLXIR','LOGLXIR_SCAT']
sav_inds = [sav_inds]


;; convert LX(LIR) to FX( LX(LIR) ) for flux-limit plots
;; erg/s to erg/s/cm^2
fxir = dblarr(nsrc)
fxir_scat = dblarr(nsrc)
logfxir = dblarr(nsrc)-9999.
logfxir_scat = dblarr(nsrc)-9999.
fxir[ilir] = lxir[ilir]/(4.*!const.pi*dl2[ilir])
fxir_scat[ilir] = lxir_scat[ilir]
logfxir[ilir] = alog10(fxir[ilir])>(-9999.)
logfxir_scat[ilir] = loglxir_scat[ilir]

sav_vars = [sav_vars,'FXIR','FXIR_SCAT','LOGFXIR','LOGFXIR_SCAT']
sav_inds = [sav_inds]



;;----------------------------------------------------------------------------------------
;; LUMINOSITIES
;;----------------------------------------------------------------------------------------
;; 2-10keV luminosity arrays
lx = 'LX'+xfield
e_lx = 'E_'+lx
loglx = 'LOG'+lx
e_loglx = 'E_'+loglx
;; Chandra, XMM, NuSTAR
cat_gamma = [1.8,1.8,1.8]
for i = 0,nfield-1 do begin
    re = execute(lx[i]+' = dblarr(nsrc)')
    re = execute(e_lx[i]+' = dblarr(nsrc)')
    re = execute(loglx[i]+' = dblarr(nsrc)-9999.')
    re = execute(e_loglx[i]+' = dblarr(nsrc)-9999.')
    re = execute('iivalid = IIDET'+xfield[i])
    ivalid = where(iivalid,validct)
    if (validct gt 0.) then begin
        ;; K correct flux to rest-frame, f_kcorr = f_obs*(1+z)^(Γ-2)
        re = execute('fx_kcorr = FX'+xfield[i]+'[ivalid]*(1+z[ivalid])^(cat_gamma[i]-2.)')
        re = execute(lx[i]+'[ivalid] = 4.*!const.pi*dl2[ivalid]*fx_kcorr')
        re = execute(e_lx[i]+'[ivalid] = '+lx[i]+'[ivalid] * sqrt((E_FX'+xfield[i]+'[ivalid]/fx_kcorr)^2. + (red_sigm[1,ivalid]/z[ivalid])^2.)')
        ;; check resamp distribution for 0 (MEDABSDEV==0), -1 (only one source), or -9999 (should not ever be the case; sanity check)
        iest = where(iivalid and red_sigm[1,*] le 0.,estct)
        if (estct gt 0) then re = execute(e_lx[i]+'[iest] = '+lx[i]+'[iest] * sqrt((E_FX'+xfield[i]+'[iest]/fx_kcorr)^2. + (zerr[iest]/z[iest])^2.)')
        ;if (estct gt 0) then re = execute(e_lx[i]+'[iest] = '+lx[i]+'[iest]*0.1')
        re = execute(loglx[i]+'[ivalid] = alog10('+lx[i]+'[ivalid])>(-9999.)')
        re = execute(e_loglx[i]+'[ivalid] = '+e_lx[i]+'[ivalid]/(alog(10.)*'+lx[i]+'[ivalid])>(-9999.)')
    endif
endfor

sav_vars = [sav_vars,'CAT_GAMMA',lx,e_lx,loglx,e_loglx]
sav_inds = [sav_inds]


;;----------------------------------------------------------------------------------------
;; LIMITS
;;----------------------------------------------------------------------------------------
;; flux limit at source
fxlim = 'FXLIM'+xfield
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
    re = execute(e_logfxlim[i]+' = dblarr(nsrc)-9999.')
    re = execute(lxlim[i]+' = dblarr(nsrc)')
    re = execute(e_lxlim[i]+' = dblarr(nsrc)')
    re = execute(loglxlim[i]+' = dblarr(nsrc)-9999.')
    re = execute(e_loglxlim[i]+' = dblarr(nsrc)-9999.')
    re = execute(fxlim_cs[i]+' = dblarr(degr[i])')
    re = execute('iivalid = IINON'+xfield[i])
    ivalid = where(iivalid,validct)
    if (validct gt 0) then begin
        re = execute(fxlim[i]+'[ivalid] = extrapolate_flim(CAT_LIM'+xfield[i]+',CAT_EXP'+xfield[i]+',TEXP'+xfield[i]+'[ivalid],degr[i],FLIM_CS='+fxlim_cs[i]+')')
        re = execute(e_fxlim[i]+'[ivalid] = '+fxlim[i]+'[ivalid] / 3d')
        re = execute(logfxlim[i]+'[ivalid] = alog10('+fxlim[i]+'[ivalid])>(-9999.)')
        re = execute(e_logfxlim[i]+'[ivalid] = '+e_fxlim[i]+'[ivalid]/(alog(10)*'+fxlim[i]+'[ivalid])>(-9999.)')
        re = execute(lxlim[i]+'[ivalid] = 4.*!const.pi*dl2[ivalid]*'+fxlim[i]+'[ivalid]')
        re = execute(e_lxlim[i]+'[ivalid] = '+lxlim[i]+'[ivalid] * sqrt(('+e_fxlim[i]+'[ivalid]/'+fxlim[i]+'[ivalid])^2. + (red_sigm[1,ivalid]/z[ivalid])^2.)')
        ;; check resamp distribution for 0 (MEDABSDEV==0), -1 (only one source), or -9999 (should not ever be the case; sanity check)
        iest = where(iivalid and red_sigm[1,*] le 0.,estct)
        if (estct gt 0.) then re = execute(e_lxlim[i]+'[iest] = '+lxlim[i]+'[iest] * sqrt(('+e_fxlim[i]+'[iest]/'+fxlim[i]+'[iest])^2. + (zerr[iest]/z[iest])^2.)')
        ;re = execute(e_lxlim[i]+'[iest] = '+lxlim[i]+'[iest]*0.1')
        re = execute(loglxlim[i]+'[ivalid] = alog10('+lxlim[i]+'[ivalid])>(-9999.)')
        re = execute(e_loglxlim[i]+'[ivalid] = '+e_lxlim[i]+'[ivalid]/(alog(10)*'+lxlim[i]+'[ivalid])>(-9999.)')
    endif
endfor

sav_vars = [sav_vars,'DEGR',fxlim_cs,fxlim,e_fxlim,logfxlim,e_logfxlim,lxlim,e_lxlim,loglxlim,e_loglxlim]
sav_inds = [sav_inds]


;;----------------------------------------------------------------------------------------
;; RATIOS
;;----------------------------------------------------------------------------------------
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
    re = execute('iivalid = IIDET'+xfield[i]+' and lxir gt 0.')
    ivalid = where(iivalid,detct)
    if (detct gt 0.) then begin
        ;; start in linear space
        re = execute(lldet[i]+'[ivalid] = lx'+xfield[i]+'[ivalid]/(lxir[ivalid]*lxir_scat[ivalid])')
        ;; errors attributed to X-ray flux and IR flux
        re = execute(e_lldet[i]+'[ivalid] = '+lldet[i]+'[ivalid] * sqrt((E_FX'+xfield[i]+'[ivalid]/FX'+xfield[i]+'[ivalid])^2. + (fir_sigm[1,ivalid]/fir_sigm[0,ivalid])^2.)')
        ;; check resamp distribution for 0 (MEDABSDEV==0), -1 (only one source), or -9999 (should not ever be the case where AGN component; sanity check)
        iest = where(iivalid and fir_sigm[1,*] le 0.,estct)
        if (estct gt 0) then stop
        ;; convert to log space
        re = execute(e_lldet[i]+'[ivalid] /= (alog(10.)*'+lldet[i]+'[ivalid])>(-9999.)')
        re = execute(lldet[i]+'[ivalid] = alog10('+lldet[i]+'[ivalid])>(-9999.)')
    endif
    ;; non-detections
    re = execute('iivalid = IINON'+xfield[i]+' and lxir gt 0.')
    ivalid = where(iivalid,nonct)
    if (nonct gt 0.) then re = execute(llnon[i]+'[ivalid] = LOGLXLIM'+xfield[i]+'[ivalid]-loglxir[ivalid]')
endfor

sav_vars = [sav_vars,lldet,e_lldet,llnon,e_llnon]
sav_inds = [sav_inds]


sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="src_luminosities.sav"')


END


