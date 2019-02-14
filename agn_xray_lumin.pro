PRO agn_xray_lumin, DERED = dered, $
                    PLT = plt, $
                    COMBINE = combine


common _inf_fits
common _inf_nst
common _inf_xmm
common _inf_cha
common _det_nst
common _det_xmm
common _det_cha
common _det_wac
common _soft210
common _soft052
common _flx_lim
common _comp


;; fit output 
ebv = reform(param_inf[0,*])
c_agn = reform(param_inf[2,*])
iilir = c_agn gt 0.                 ;; AGN SED contribution exists
ilir = where(iilir)

sav_vars = ['EBV']
sav_inds = ['IILIR','ILIR']


;;----------------------------------------------------------------------------------------
;; LUMINOSITIES -- L(IR) AND LX(LIR)
;;----------------------------------------------------------------------------------------
;; IR 6-micron AGN luminosity calculated from SED model parameters
l15obs = l_agn(15.,ebv,z_inf,c_agn,/log) > 0.   ;; 15-micron observed
if keyword_set(dered) then begin
    l15int = l_agn(15.,dblarr(nsrc),z_inf,c_agn,/log) > 0.  ;; 15-micron intrinsic
    l06int = l_agn(6.,dblarr(nsrc),z_inf,c_agn,/log) > 0.   ;; 6-micron intrinsic
    lir = l06int                                            ;; 6-micron intrinsic
    ;conv = median(l15int[ilir] - l06int[ilir])              ;; slope 15- to 6-micron
    ;lir[ilir] -= conv
endif else $
    lir = l_agn(6.,ebv,z_inf,c_agn,/log) > 0    ;; 6-micron observed

lxir_210 = dblarr(nsrc)                     ;; unobscured LX given L6um
lhi_cut = 44.79                             ;; luminosity turnover based on LX-L6um relationship
ilo = where(iilir and lir lt lhi_cut)       ;; LIR exists and is below turnover
ihi = where(iilir and lir gt lhi_cut)       ;; LIR exists and is gt turnover, redundant but sanity check
lxir_210[ilo] = 0.84*(lir[ilo]-45.)+44.60   ;; for LIR lt 44.79 erg s-1
lxir_210[ihi] = 0.40*(lir[ihi]-45.)+44.51   ;; for LIR ge 44.79 erg s-1

;; LX(LIR) luminosity conversion 2-10 to 0.5-2 keV
lxir_052 = dblarr(nsrc)
out052_in210 = 3.300E-01
lxir_052[ilir] = lxir_210[ilir] + alog10(out052_in210)

sav_vars = [sav_vars,'LIR','LXIR_210','LXIR_052','OUT052_IN210']
sav_inds = [sav_inds]



;;----------------------------------------------------------------------------------------
;; LUMINOSITIES -- LX
;;----------------------------------------------------------------------------------------
;; luminosity distance
dl = dlum(z_inf)

;; X-ray luminosity converted from flux
lx_210 = 'LX'+xfield+'_210'
for i = 0,nfield-1 do begin
    re = execute(lx_210[i]+' = dblarr(nsrc)')
    re = execute(lx_210[i]+'['+idet_210[i]+'] = alog10((4.*!const.pi*dl['+idet_210[i]+']^2)*'+xray_flx_210[i]+'['+idet_210[i]+'])>0.')
endfor

lx_052 = 'LX'+xfield+'_052'
for i = 0,nfield-1 do begin
    re = execute(lx_052[i]+' = dblarr(nsrc)')
    re = execute(lx_052[i]+'['+idet_052[i]+'] = alog10((4.*!const.pi*dl['+idet_052[i]+']^2)*'+xray_flx_052[i]+'['+idet_052[i]+'])>0.')
endfor

sav_vars = [sav_vars,'LX_210',lx_210, $
                     'LX_052',lx_052]
sav_inds = [sav_inds]


;;----------------------------------------------------------------------------------------
;; X-ray flux limit luminosities at each source -- (includes upper limit luminosities)
;;----------------------------------------------------------------------------------------
;; 2-10keV
flim_210 = 'FLIM'+xfield+'_210'     ;; flux_inf limit at source
lxlim_210 = 'LXLIM'+xfield+'_210'   ;; luminosity at flux_inf limit
for i = 0,nfield-1 do begin
    re = execute(flim_210[i]+' = dblarr(nsrc)')
    re = execute(lxlim_210[i]+' = dblarr(nsrc)')
    re = execute('isrc = where(iiinf'+xfield[i]+')')
    re = execute(flim_210[i]+'[isrc] = extrap_flim('+cat_lim_210[i]+','+cat_exp_210[i]+','+texp[i]+'[isrc])')
    re = execute(lxlim_210[i]+'[isrc] = alog10(4.*!const.pi*dl[isrc]^2*'+flim_210[i]+'[isrc])>0.')
endfor

;; 0.5-2keV
flim_052 = 'FLIM'+xfield+'_052'     ;; flux_inf limit at source
lxlim_052 = 'LXLIM'+xfield+'_052'   ;; luminosity at flux_inf limit
for i = 0,nfield-1 do begin
    re = execute(flim_052[i]+' = dblarr(nsrc)')
    re = execute(lxlim_052[i]+' = dblarr(nsrc)')
    re = execute('isrc = where(iiinf'+xfield[i]+')')
    re = execute(flim_052[i]+'[isrc] = extrap_flim('+cat_lim_052[i]+','+cat_exp_052[i]+','+texp[i]+'[isrc])')
    re = execute(lxlim_052[i]+'[isrc] = alog10(4.*!const.pi*dl[isrc]^2*'+flim_052[i]+'[isrc])>0.')
endfor

sav_vars = [sav_vars,'FLIM_210','LXLIM_210',flim_210,lxlim_210, $
                     'FLIM_052','LXLIM_052',flim_052,lxlim_052]
sav_inds = [sav_inds]


;;----------------------------------------------------------------------------------------
;; SOURCE QUALITY CUTS
;;----------------------------------------------------------------------------------------
iiz = z_inf gt 0. and z_inf lt 0.6
;; ensure WISE photometry
sn_wise = flux_inf[where(strmatch(band_inf,'WISE*')),*]/e_flux_inf[where(strmatch(band_inf,'WISE*')),*]
totsn = total(sn_wise ge 1.,1)          ;; all WISE photometry must exist and S/N ³ 1
iisn = totsn eq 4.                      ;; note: all non-finite sn_wise == -NaN
;; SED chi-square 
chi = reform(param_inf[-2,*])
dof = reform(param_inf[-1,*])
rchi = chi/dof
iichi = rchi le 20.
;; passes all quality cuts
iiq = iiz and iisn; and iichi

;; INDEXING FOR ANALYSIS
iiirb = l15obs ge 44. and lir ge 44.    ;; restrict to IR-bright sources, constrain to IR-bright 15-micron observed
iiagn = iiirb and iiq                   ;; IR-bright AGN that pass quality cuts

iiagn_det = 'IIAGN_DET'+xfield             ;; Candidate AGN w.  X-ray detections
iiagn_lim = 'IIAGN_LIM'+xfield             ;; Candidate AGN w.o X-ray detections
iagn_det = 'IAGN_DET'+xfield
iagn_lim = 'IAGN_LIM'+xfield

iiagn_det_210 = iiagn_det+'_210'
iiagn_lim_210 = iiagn_lim+'_210'
iagn_det_210 = iagn_det+'_210'
iagn_lim_210 = iagn_lim+'_210'

for i = 0,nfield-1 do begin
    re = execute(iiagn_det_210[i]+' = iiagn and '+iidet_210[i])
    re = execute(iiagn_lim_210[i]+' = iiagn and '+iiinf[i])
    re = execute(iagn_det_210[i]+' = where('+iiagn_det_210[i]+')')
    re = execute(iagn_lim_210[i]+' = where('+iiagn_lim_210[i]+')')
endfor

iiagn_det_052 = iiagn_det+'_052'
iiagn_lim_052 = iiagn_lim+'_052'
iagn_det_052 = iagn_det+'_052'
iagn_lim_052 = iagn_lim+'_052'

for i = 0,nfield-1 do begin
    re = execute(iiagn_det_052[i]+' = iiagn and '+iidet_052[i])
    re = execute(iiagn_lim_052[i]+' = iiagn and '+iiinf[i])
    re = execute(iagn_det_052[i]+' = where('+iiagn_det_052[i]+')')
    re = execute(iagn_lim_052[i]+' = where('+iiagn_lim_052[i]+')')
endfor

sav_vars = [sav_vars,'SN_WISE','TOTSN','CHI','DOF','RCHI']
sav_inds = [sav_inds,'IIZ','IISN','IICHI','IIQ','IIIRB','IIAGN', $
                     'IIAGN_DET','IIAGN_LIM','IAGN_DET','IAGN_LIM', $
                     'IIAGN_DET_210','IIAGN_LIM_210','IAGN_DET_210','IAGN_LIM_210', $
                     iiagn_det_210,iiagn_lim_210,iagn_det_210,iagn_lim_210, $
                     'IIAGN_DET_052','IIAGN_LIM_052','IAGN_DET_052','IAGN_LIM_052', $
                     iiagn_det_052,iiagn_lim_052,iagn_det_052,iagn_lim_052]

;;----------------------------------------------------------------------------------------
;; LUMINOSITY RATIOS -- PROXY FOR OBSCURATION
;;----------------------------------------------------------------------------------------
;; detections/limits luminosity ratios
lldet_210 = 'LLDET'+xfield+'_210'
lllim_210 = 'LLLIM'+xfield+'_210'
for i = 0,nfield-1 do begin
    re = execute(lldet_210[i]+' = dblarr(nsrc)-9999.')
    re = execute(lllim_210[i]+' = dblarr(nsrc)-9999.')
    re = execute(lldet_210[i]+'['+iagn_det_210[i]+'] = '+lx_210[i]+'['+iagn_det_210[i]+']-lxir_210['+iagn_det_210[i]+']')
    re = execute(lllim_210[i]+'['+iagn_lim_210[i]+'] = '+lxlim_210[i]+'['+iagn_lim_210[i]+']-lxir_210['+iagn_lim_210[i]+']')
endfor

lldet_052 = 'LLDET'+xfield+'_052'
lllim_052 = 'LLLIM'+xfield+'_052'
for i = 0,nfield-1 do begin
    re = execute(lldet_052[i]+' = dblarr(nsrc)-9999.')
    re = execute(lllim_052[i]+' = dblarr(nsrc)-9999.')
    re = execute(lldet_052[i]+'['+iagn_det_052[i]+'] = '+lx_052[i]+'['+iagn_det_052[i]+']-lxir_052['+iagn_det_052[i]+']')
    re = execute(lllim_052[i]+'['+iagn_lim_052[i]+'] = '+lxlim_052[i]+'['+iagn_lim_052[i]+']-lxir_052['+iagn_lim_052[i]+']')
endfor


sav_vars = [sav_vars,'LLDET_210','LLLIM_210',lldet_210,lllim_210, $
                     'LLDET_052','LLLIM_052',lldet_052,lllim_052]
sav_inds = [sav_inds]

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="xray_lum.sav"')



;;----------------------------------------------------------------------------------------
;; Luminosity data subset by field
;;----------------------------------------------------------------------------------------
if keyword_set(combine) then begin

    lldet = lldet_xmm_210
    lldet[where(lldet eq -9999.)] = lldet_cha_210[where(lldet eq -9999.)]
    lldet[where(lldet eq -9999.)] = lldet_nst_210[where(lldet eq -9999.)]
    lllim = lllim_xmm_210
    lllim[where(lllim eq -9999.)] = lllim_cha_210[where(lllim eq -9999.)]
    lllim[where(lllim eq -9999.)] = lllim_nst_210[where(lllim eq -9999.)]
    ;; combined detections/limits indices
    ixdet = where(lldet gt -9999.)
    ixlim = where(lllim gt -9999.)
    ;; AGN Catalog matches
    iwdet = where(lldet gt -9999. and iiwac)
    iwlim = where(lllim gt -9999. and iiwac)
    ;; log E(B-V)    
    lebv = alog10(ebv)>(-2.5)
    lebv = lebv + randomu(seed,n_elements(ebv))*0.05-0.05/2.

    sav_vars = [sav_vars,'LLDET','LLLIM','LEBV']
    sav_inds = [sav_inds,'IXDET','IXLIM','IWDET','IWLIM']

    sav_str = strjoin([sav_vars,sav_inds],',')
    re = execute('save,'+sav_str+',/compress,file="xray_lum.sav"')
 
 endif

;;----------------------------------------------------------------------------------------
;; PLOTS!!
;;----------------------------------------------------------------------------------------
if keyword_set(plt) then begin
    
    e = {symbol:'o',sym_size:0.5,sym_filled:1,sym_transparency:50,linestyle:'', $
         ;xra:[5e-3,1e2],xlog:1,yra:[-4,2],ylog:0, $
         xra:[-3.5,2.5],xlog:0,yra:[-6,4],ylog:0, $
         xtitle:'$log E(B-V)_{AGN}$',ytitle:'$log( L_{2-10} keV} / L_{2-10 keV}(L_{IR}) )$', $
         buffer:0}
    plim = plot(lebv[ixlim],lllim[ixlim],col='orange',_extra=e,name='AGN X-ray non-det.')
    pdet = plot(lebv[ixdet],lldet[ixdet],col='dodger blue',/ov,_extra=e,name='AGN + X-ray')
    pwac = plot(lebv[iwlim],lllim[iwlim],'or',sym_size=0.25,sym_filled=1,/ov)
    
    plim = plot(lebv[iwlim],lllim[iwlim],col='orange',_extra=e,name='AGN X-ray non-det.')
    pdet = plot(lebv[iwdet],lldet[iwdet],col='dodger blue',/ov,_extra=e,name='AGN + X-ray')
    
    
    lebv = alog10(ebv)>(-2.5)
    lebv = lebv + randomu(seed,n_elements(ebv))*0.05-0.05/2.
    ebvd_210 = {NST:lebv[iagn_det_nst_210],XMM:lebv[iagn_det_xmm_210],CHA:lebv[iagn_det_cha_210]}
    ebvl_210 = {NST:lebv[iagn_lim_nst_210],XMM:lebv[iagn_lim_xmm_210],CHA:lebv[iagn_lim_cha_210]}
    ebvw_210 = {NST:lebv[where(iiagn_lim_nst_210 and iiwac)],XMM:lebv[where(iiagn_lim_xmm_210 and iiwac)],CHA:lebv[where(iiagn_lim_cha_210 and iiwac)]}
    ld_210 = {NST:lldet_nst_210[iagn_det_nst_210],XMM:lldet_xmm_210[iagn_det_xmm_210],CHA:lldet_cha_210[iagn_det_cha_210]}
    ll_210 = {NST:lllim_nst_210[iagn_lim_nst_210],XMM:lllim_xmm_210[iagn_lim_xmm_210],CHA:lllim_cha_210[iagn_lim_cha_210]}
    lw_210 = {NST:lllim_nst_210[where(iiagn_lim_nst_210 and iiwac)],XMM:lllim_xmm_210[where(iiagn_lim_xmm_210 and iiwac)],CHA:lllim_cha_210[where(iiagn_lim_cha_210 and iiwac)]}
    ebvsub = ['EBVD_210','EBVL_210','EBVW_210']
    llsub = ['LD_210','LL_210','LW_210']

    e = {symbol:'o',sym_size:0.5,sym_filled:1,sym_transparency:50,linestyle:'', $
         ;xra:[5e-3,1e2],xlog:1,yra:[-6,4],ylog:0, $
         xra:[-3.5,2.5],xlog:0,yra:[-6,4],ylog:0, $
         xtitle:'$log E(B-V)_{AGN}$',ytitle:'$log( L_{0.5-2 keV} / L_{0.5-2 keV}(L_{IR}) )$', $
         buffer:0}
    plim = plot(ebvl_210.nst,ll_210.nst,col='orange',_extra=e,name='AGN X-ray non-det.')
    !NULL = plot(ebvl_210.xmm,ll_210.xmm,col='orange',/ov,_extra=e)
    !NULL = plot(ebvl_210.cha,ll_210.cha,col='orange',/ov,_extra=e)
    pdx = plot(ebvd_210.xmm,ld_210.xmm,col='medium blue',/ov,_extra=e,name='AGN + XMM')
    pdc = plot(ebvd_210.cha,ld_210.cha,col='dodger blue',/ov,_extra=e,name='AGN + Chandra')
    pdn = plot(ebvd_210.nst,ld_210.nst,col='sky blue',/ov,_extra=e,name='AGN + NuSTAR')
    !NULL = plot(e.xra,[0,0],'--',/ov,_extra=e)
    pwac = plot(ebvw_210.nst,lw_210.nst,'or',sym_size=0.75,/ov,name='WISE AGN 17')
    !NULL = plot(ebvw_210.xmm,lw_210.xmm,'or',sym_size=0.75,/ov)
    !NULL = plot(ebvw_210.cha,lw_210.cha,'or',sym_size=0.75,/ov)
    leg = legend(target=[plim,pdn,pdx,pdc,pwac],position=[1.05,0.87],/auto_text_color,sample_width=0,font_size=8)

    xnhlines = rebin([[-3.3],[-2.8]],4,2)
    ynhlines = rebin([-0.264805,-1.03190,-1.52835,-2.06607],4,2)
    nhlabel = ['$N_H=10^{24}cm^{-2}$','$    =5x10^{24}cm^{-2}$','$    =10^{25}cm^{-2}$','$    =5x10^{25}cm^{-2}$']
    for i = 0,3 do p = plot(xnhlines[i,*],ynhlines[i,*],':',thick=4,/ov)
    e = {font_size:8,font_style:'bold',data:1}
    for i = 0,3 do t = text(xnhlines[i,0]-0.15,ynhlines[i,0]-0.4,nhlabel[i],_extra=e)

    
endif


END







;;----------------------------------------------------------------------------------------
;; NH DISTRIBUTIONS -- RAW UNMODELED OBSCURATION
;;----------------------------------------------------------------------------------------
if keyword_set(nh) then begin
nhdet_210 = 'NHDET'+xfield+'_210'
nhlim_210 = 'NHLIM'+xfield+'_210'
;nhdet_052 = 'NHDET'+xfield+'_052'
;nhlim_052 = 'NHLIM'+xfield+'_052'
for i = 0,nfield-1 do begin
    re = execute(nhdet_210[i]+' = dblarr(nsrc)')
    re = execute(nhlim_210[i]+' = dblarr(nsrc)')
    eband = '2-10'
    re = execute(nhdet_210[i]+'['+iagn_det[i]+'] = ll2nh('+lldet_210[i]+'['+iagn_det[i]+'],eband)')
    re = execute(nhlim_210[i]+'['+iagn_lim[i]+'] = ll2nh('+lllim_210[i]+'['+iagn_lim[i]+'],eband)')
;    re = execute(nhdet_052[i]+' = dblarr(nsrc)')
;    re = execute(nhlim_052[i]+' = dblarr(nsrc)')
;   eband = '0.5-2'
;    re = execute(nhdet_052[i]+'['+iagn_det[i]+'] = ll2nh('+lldet_052[i]+'['+iagn_det[i]+'],eband)')
;    re = execute(nhlim_052[i]+'['+iagn_lim[i]+'] = ll2nh('+lllim_052[i]+'['+iagn_lim[i]+'],eband)')
endfor

sav_vars = [sav_vars,'NHDET_210','NHLIM_210',nhdet_210,nhlim_210];,'NHDET_052','NHLIM_052',nhdet_052,nhlim_052]
sav_inds = [sav_inds]
endif







