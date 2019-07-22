PRO agn_xray_lumin, DERED = dered, $
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
sav_inds = ['IILIR','ILIR']



;;----------------------------------------------------------------------------------------
;; LUMINOSITIES -- L(IR) AND LX(LIR)
;;----------------------------------------------------------------------------------------
;; IR 6-micron AGN luminosity calculated from SED model parameters
l15obs = l_agn(15.,ebv,z,c_agn,/log) > 0.   ;; 15-micron observed
if keyword_set(dered) then begin
    l06int = l_agn(6.,dblarr(nsrc),z,c_agn,/log)    ;; 6-micron intrinsic
    lir = l06int                                    ;; 6-micron intrinsic
endif else $
    lir = l_agn(6.,ebv,z,c_agn,/log)                ;; 6-micron observed

;; correct AGN luminosity where template over- or underestimates beyond uncertainties
;; luminosity corrected
;iilumc = iilir and iiirc and iichi; and agnf6.obs gt 0.7;((ebv lt 0.2 and agnf6.obs gt 0.9) or (ebv gt 50. and agnf6.obs gt 0.9))
;ilumc = where(iilumc)
;lir[ilumc] = correct_agn_lum(lir[ilumc],wave,flux[*,ilumc],e_flux[*,ilumc],param[*,ilumc],z[ilumc],agnf6[ilumc].obs,/over,/under,NCORR=ncorr)
lir = correct_agn_lum(lir,wave,flux,e_flux,param,z,agnf6.obs,/over,/under,NCORR=ncorr)

;; LX as a function of LIR, LX-LIR relationship
lxir_210 = dblarr(nsrc)                     ;; unobscured LX given L6um
lcut = 44.79                                ;; luminosity turnover based on LX-L6um relationship
ilo = where(iilir and lir lt lcut)          ;; LIR exists and is below turnover
ihi = where(iilir and lir ge lcut)          ;; LIR exists and is gt turnover, redundant but sanity check
lxir_210[ilo] = 0.84*(lir[ilo]-45.)+44.60   ;; for LIR < 44.79 erg s-1
lxir_210[ihi] = 0.40*(lir[ihi]-45.)+44.51   ;; for LIR > 44.79 erg s-1

;; convert LX(LIR) to FX( LX(LIR) ) for flux-limit plots
;; luminosity distance
dl2 = dlum(z,/squared)
fxir_210 = 10.^(lxir_210)/(4.*!const.pi*dl2)


sav_vars = [sav_vars,'NCORR','LIR','LXIR_210','DL2','FXIR_210']
sav_inds = [sav_inds]



;;----------------------------------------------------------------------------------------
;; LUMINOSITIES -- LX
;;----------------------------------------------------------------------------------------
;; X-ray luminosity converted from flux, K-corrected to rest-frame Fkcor = Fobs*(1+z)^(Γ-2)
;; Chandra, XMM, NuSTAR
cat_gamma = [2.0,1.8,1.8]
;; 2-10keV
lx_210 = 'LX'+xfield+'_210'
for i = 0,nfield-1 do begin
    re = execute(lx_210[i]+' = dblarr(nsrc)')
    re = execute('idet = where(iidet'+xfield[i]+',ctdet)')
    if (ctdet eq 0.) then continue
    re = execute(lx_210[i]+'[idet] = alog10((4.*!const.pi*dl2[idet])*('+xray_flx_210[i]+'[idet]*(1+z[idet])^(cat_gamma[i]-2.)))>0.')
endfor


sav_vars = [sav_vars,'CAT_GAMMA','LX_210', $
                                 lx_210]
sav_inds = [sav_inds]


;;----------------------------------------------------------------------------------------
;; INTERPOLATE FLUX LIMIT FLUX & LUMINOSITY
;;----------------------------------------------------------------------------------------
;; 2-10keV
flim_210 = 'FLIM'+xfield+'_210'     ;; flux limit at source
lxlim_210 = 'LXLIM'+xfield+'_210'   ;; luminosity at flux limit
for i = 0,nfield-1 do begin
    re = execute(flim_210[i]+' = dblarr(nsrc)')
    re = execute(lxlim_210[i]+' = dblarr(nsrc)')
    re = execute('isrc = where(iiinf'+xfield[i]+')')
    re = execute(flim_210[i]+'[isrc] = extrap_flim('+cat_lim_210[i]+','+cat_exp_210[i]+','+texp[i]+'[isrc])')
    re = execute(lxlim_210[i]+'[isrc] = alog10(4.*!const.pi*dl2[isrc]*'+flim_210[i]+'[isrc])>0.')
endfor


sav_vars = [sav_vars,'FLIM_210','LXLIM_210', $
                     flim_210,lxlim_210]
sav_inds = [sav_inds]


sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="xray_lum.sav"')





END


