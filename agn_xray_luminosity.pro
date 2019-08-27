PRO agn_xray_luminosity, DERED = dered, $
                         PLT = plt, $
                         COMBINE = combine


common _fits
common _inf_cha
common _inf_xmm
common _inf_nst
common _det_nst
common _det_xmm
common _det_cha
common _softx
common _fluxlim
common _comp


;; fit output 
ebv = reform(param[0,*])
c_agn = reform(param[2,*])
iilir = c_agn gt 0.                 ;; AGN SED contribution exists
ilir = where(iilir,nagn,ncomplement=ngal,/null)

;;----------------------------------------------------------------------------------------
;; AGN FRACTION -- 
;;----------------------------------------------------------------------------------------
;; calculate AGN fraction at 6- and 15-micron
agnf6 = f_agn(6.,param,model=agnm6)
agnf15 = f_agn(15.,param,model=agnm15)


sav_vars = ['EBV','NAGN','NGAL','AGNF6','AGNM6','AGNF15','AGNM15']
sav_inds = ['IILIR']



;;----------------------------------------------------------------------------------------
;; LUMINOSITIES -- L(IR) AND LX(LIR)
;;----------------------------------------------------------------------------------------
;; IR 6-micron AGN luminosity calculated from SED model parameters
lir = dblarr(nsrc)
if keyword_set(dered) then begin
    l06int = l_agn(6.,dblarr(nagn),z[ilir],c_agn[ilir],/log)    ;; 6-micron intrinsic
    lir[ilir] = l06int                                          ;; 6-micron intrinsic
endif else $
    lir[ilir] = l_agn(6.,ebv[ilir],z[ilir],c_agn[ilir],/log)    ;; 6-micron observed

;; correct AGN luminosity where template over- or underestimates beyond uncertainties
;; luminosity corrected
;iilumc = iilir and iiirc and iichi; and agnf6.obs gt 0.7;((ebv lt 0.2 and agnf6.obs gt 0.9) or (ebv gt 50. and agnf6.obs gt 0.9))
;ilumc = where(iilumc)
;lir[ilumc] = correct_agn_lum(lir[ilumc],wave,flux[*,ilumc],e_flux[*,ilumc],param[*,ilumc],z[ilumc],agnf6[ilumc].obs,/over,/under,NCORR=ncorr)
lir = correct_agn_lum(lir,wave,flux,e_flux,param,z,agnf6.obs,/over,/under,IICORR=iicorr)
;; LX as a function of LIR, LX-LIR relationship
lxir = dblarr(nsrc)                     ;; unobscured LX given L6um
lcut = 44.79                                ;; luminosity turnover based on LX-L6um relationship
ilo = where(iilir and lir lt lcut)          ;; LIR exists and is below turnover
ihi = where(iilir and lir ge lcut)          ;; LIR exists and is gt turnover, redundant but sanity check
lxir[ilo] = 0.84*(lir[ilo]-45.)+44.60   ;; for LIR < 44.79 erg s-1
lxir[ihi] = 0.40*(lir[ihi]-45.)+44.51   ;; for LIR > 44.79 erg s-1

;; convert LX(LIR) to FX( LX(LIR) ) for flux-limit plots
;; erg/s to erg/s/cm^2
;; luminosity distance
dl2 = dlum(z,/squared)
fxir = dblarr(nsrc)
fxir[ilir] = 10.^(lxir[ilir])/(4.*!const.pi*dl2[ilir])


sav_vars = [sav_vars,'LIR','LXIR','DL2','FXIR']
sav_inds = [sav_inds,'IICORR']



;;----------------------------------------------------------------------------------------
;; LUMINOSITIES -- LX
;;----------------------------------------------------------------------------------------
;; X-ray luminosity converted from flux, K-corrected to rest-frame Fkcor = Fobs*(1+z)^(Γ-2)
;; Chandra, XMM, NuSTAR
cat_gamma = [2.0,1.8,1.8]
;; 2-10keV
lx = 'LX'+xfield
for i = 0,nfield-1 do begin
    re = execute(lx[i]+' = dblarr(nsrc)')
    re = execute('idet = where(IIDET'+xfield[i]+',ctdet)')
    if (ctdet eq 0.) then continue
    re = execute(lx[i]+'[idet] = alog10((4.*!const.pi*dl2[idet])*(FLX'+xfield[i]+'[idet]*(1+z[idet])^(cat_gamma[i]-2.)))>0.')
endfor


sav_vars = [sav_vars,'CAT_GAMMA',lx]
sav_inds = [sav_inds]


;;----------------------------------------------------------------------------------------
;; INTERPOLATE FLUX LIMIT FLUX & LUMINOSITY
;;----------------------------------------------------------------------------------------
;; 2-10keV
flim = 'FLIM'+xfield        ;; flux limit at source
lxlim = 'LXLIM'+xfield      ;; luminosity at flux limit
flim_cs = 'FLIM_CS'+xfield  ;; flux-limit function coefficients 
degr = [6,6,1]              ;; degree of polynomial to fit flux-limit
for i = 0,nfield-1 do begin
    re = execute(flim[i]+' = dblarr(nsrc)')
    re = execute(lxlim[i]+' = dblarr(nsrc)')
    re = execute(flim_cs[i]+' = dblarr(degr[i])')
    re = execute('isrc = where(iiinf'+xfield[i]+')')
    re = execute(flim[i]+'[isrc] = extrapolate_flim(CAT_LIM'+xfield[i]+',CAT_EXP'+xfield[i]+',TEXP'+xfield[i]+'[isrc],degr[i],FLIM_CS='+flim_cs[i]+')')
    re = execute(lxlim[i]+'[isrc] = alog10(4.*!const.pi*dl2[isrc]*FLIM'+xfield[i]+'[isrc])>0.')
endfor


sav_vars = [sav_vars,flim,lxlim,flim_cs,'DEGR']
sav_inds = [sav_inds]


sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="src_luminosity.sav"')





END


