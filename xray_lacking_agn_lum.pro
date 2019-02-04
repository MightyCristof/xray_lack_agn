PRO xray_lacking_agn_lum, DERED = dered, $
                          PLT = plt


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
common _comp


;;----------------------------------------------------------------------------------------
;; X-RAY AND IR LUMINOSITIES
;;----------------------------------------------------------------------------------------
;; IR 6-micron AGN luminosity calculated from SED model parameters
;; fit output 
ebv = reform(param_inf[0,*])
c_agn = reform(param_inf[2,*])
iilir = c_agn gt 0.         ;; AGN SED contribution exists
ilir = where(iilir)

l15obs = l_agn(15.,ebv,z_inf,c_agn,/log) > 0.   ;; 15-micron observed
l06obs = l_agn(6.,ebv,z_inf,c_agn,/log) > 0.    ;; 6-micron observed
if keyword_set(dered) then begin
    l15int = l_agn(15.,dblarr(nsrc),z_inf,c_agn,/log) > 0.  ;; 15-micron intrinsic
    l06int = l_agn(6.,dblarr(nsrc),z_inf,c_agn,/log) > 0.   ;; 6-micron intrinsic
    conv = median(l15int[ilir] - l06int[ilir])          ;; slope 15- to 6-micron
    lir = l06int                                        ;; 6-micron intrinsic
    lir[ilir] -= conv
endif else $
    lir = l_agn(6.,ebv,z_inf,c_agn,/log) > 0                ;; 6-micron redenned

;; X-ray luminosity converted from flux_inf
;; calculate LX(L6um) from Chen+17
dl = dlum(z_inf)
;lx = 'LX_'+['210','052']
lx_210 = 'LX_210'+xfield
lx_052 = 'LX_052'+xfield
for i = 0,nfield-1 do begin
    re = execute(lx_210[i]+' = dblarr(nsrc)')
    re = execute(lx_210[i]+'['+xdet[i]+'] = alog10((4.*!const.pi*dl['+xdet[i]+']^2)*'+xray_flx_210[i]+'['+xdet[i]+'])>0.')
    re = execute(lx_052[i]+' = dblarr(nsrc)')
    re = execute(lx_052[i]+'['+xdet[i]+'] = alog10((4.*!const.pi*dl['+xdet[i]+']^2)*'+xray_flx_052[i]+'['+xdet[i]+'])>0.')
endfor

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

sav_vars = ['EBV','LIR','DL','LX_210','LX_052',lx_210,lx_052,'LXIR_210','LXIR_052']
sav_inds = ['IILIR','ILIR']


;;----------------------------------------------------------------------------------------
;; X-RAY flux_inf LIMITS
;;----------------------------------------------------------------------------------------
cat_fld = (strsplit(xfield,'_',/extract)).ToArray()
raw_exp = cat_fld+'.'+use_exp
raw_flx = cat_fld+'.'+(strsplit(use_flx,"_210",/regex,/extract)).ToArray()
raw_err = cat_fld+'.'+(strsplit(use_err,"_210",/regex,/extract)).ToArray()
cat_exp = 'CAT_'+use_exp
cat_flx = 'CAT_'+(strsplit(use_flx,"_210",/regex,/extract)).ToArray()
cat_err = 'CAT_'+(strsplit(use_err,"_210",/regex,/extract)).ToArray()
cat_lim = 'CAT_LIM'+xfield
for i = 0,nfield-1 do begin
    re = execute('cat_ind = where('+raw_exp[i]+' gt 0. and '+raw_flx[i]+' gt 0. and '+raw_err[i]+' gt 0.)')
    re = execute(cat_exp[i]+' = '+strjoin(strsplit(raw_exp[i],".",/extract),"[cat_ind]."))
    re = execute(cat_flx[i]+' = '+strjoin(strsplit(raw_flx[i],".",/extract),"[cat_ind]."))
    re = execute(cat_err[i]+' = '+strjoin(strsplit(raw_err[i],".",/extract),"[cat_ind]."))
endfor
;; construct Chandra error as flux_inf-lolim
icha = where(cat_fld eq 'CHA')
re = execute(cat_err[icha[0]]+' = '+cat_flx[icha[0]]+'-'+cat_err[icha[0]])

;; convert catalog to 2-10keV
cat_exp_210 = cat_exp
cat_flx_210 = cat_flx+'_210'
cat_err_210 = cat_err+'_210'
cat_lim_210 = cat_lim+'_210'
for i = 0,nfield-1 do begin
    re = execute(cat_flx_210[i]+' = '+cat_flx[i]+' * '+use_cnv[i])
    re = execute(cat_err_210[i]+' = '+cat_err[i]+' * '+use_cnv[i])
endfor

;; calculate catalog flux_inf limit at 2-10keV
dex = [4.,4.,4.]
root = [0.5,2./3.,1.]
cn = [1.,1e-2,1.]
for i = 0,nfield-1 do re = execute(cat_lim_210[i]+' = xray_flim('+cat_exp_210[i]+','+cat_flx_210[i]+','+cat_err_210[i]+',dex=dex[i],root=root[i],cn=cn[i])')
;; plot flux_inf limits
if keyword_set(plt) then begin
    xtitle = '$'+['t_{exp} NuSTAR','t_{exp} XMM','t_{exp} Chandra']+'$'
    ytitle = '$'+['F_{X,2-10keV}','F_{X,2-10keV}','F_{X,2-10keV}']+'$'
    for i = 0,nfield-1 do begin
        ;re = execute('p = errorplot('+cat_exp_210[i]+','+cat_flx_210[i]+','+cat_err_210[i]+',".",linestyle="",/xlog,/ylog,xtitle=xtitle[i],ytitle=ytitle[i])')
        re = execute('p = plot('+cat_exp_210[i]+','+cat_flx_210[i]+',".",/xlog,/ylog,xtitle=xtitle[i],ytitle=ytitle[i])')
        re = execute('p = plot('+cat_exp_210[i]+','+cat_lim_210[i]+',"--r",/ov)')
    endfor
endif

;; 2-10keV flux_inf limit at each source
flim_210 = 'FLIM'+xfield+'_210'     ;; flux_inf limit at source
lxlim_210 = 'LXLIM'+xfield+'_210'   ;; luminosity at flux_inf limit
for i = 0,nfield-1 do begin
    re = execute(flim_210[i]+' = dblarr(nsrc)')
    re = execute(lxlim_210[i]+' = dblarr(nsrc)')
    re = execute('isrc = where(iiinf'+xfield[i]+')')
    re = execute(flim_210[i]+'[isrc] = extrap_flim('+cat_lim_210[i]+','+cat_exp_210[i]+','+texp[i]+'[isrc])')
    re = execute(lxlim_210[i]+'[isrc] = alog10(4.*!const.pi*dl[isrc]^2*'+flim_210[i]+'[isrc])>0.')
endfor

;; convert 2-10keV to 0.5-2keV
flim_052 = (strsplit(flim_210,'_210',/regex,/extract)).ToArray()+'_052'
lxlim_052 = (strsplit(lxlim_210,'_210',/regex,/extract)).ToArray()+'_052'
for i = 0,nfield-1 do begin
    re = execute(flim_052[i]+' = '+flim_210[i]+' * '+econv)
    re = execute(lxlim_052[i]+' = '+lxlim_210[i]+' + alog10('+econv+')')
endfor

;sav_vars = [sav_vars,'cat_exp_210','cat_flx_210','cat_lim_210','FLIM','LXLIM', $
;                     cat_exp_210,cat_flx_210,cat_err_210,cat_lim_210,flim,lxlim]
;sav_inds = [sav_inds]
sav_vars = [sav_vars,'FLIM_210','LXLIM_210','FLIM_052','LXLIM_052', $
                     flim_210,lxlim_210,flim_052,lxlim_052]
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

for i = 0,nfield-1 do begin
    re = execute(iiagn_det[i]+' = iiagn and '+xxdet[i])
    re = execute(iiagn_lim[i]+' = iiagn and '+xxinf[i]+' and finite('+lxlim_210[i]+')')
    re = execute(iagn_det[i]+' = where('+iiagn_det[i]+')')
    re = execute(iagn_lim[i]+' = where('+iiagn_lim[i]+')')
endfor

sav_vars = [sav_vars,'SN_WISE','TOTSN','CHI','DOF','RCHI']
sav_inds = [sav_inds,'IIZ','IISN','IICHI','IIQ','IIIRB','IIAGN', $
                     'IIAGN_DET','IIAGN_LIM','IAGN_DET','IAGN_LIM', $
                     iiagn_det,iiagn_lim,iagn_det,iagn_lim]


;;----------------------------------------------------------------------------------------
;; LUMINOSITY RATIOS -- PROXY FOR OBSCURATION
;;----------------------------------------------------------------------------------------
;; detections/limits luminosity ratios
lldet_210 = 'LLDET'+xfield+'_210'
lllim_210 = 'LLLIM'+xfield+'_210'
lldet_052 = 'LLDET'+xfield+'_052'
lllim_052 = 'LLLIM'+xfield+'_052'
for i = 0,nfield-1 do begin
    re = execute(lldet_210[i]+' = dblarr(nsrc)')
    re = execute(lllim_210[i]+' = dblarr(nsrc)')
    re = execute(lldet_210[i]+'['+iagn_det[i]+'] = '+lx_210[i]+'['+iagn_det[i]+']-lxir_210['+iagn_det[i]+']')
    re = execute(lllim_210[i]+'['+iagn_lim[i]+'] = '+lxlim_210[i]+'['+iagn_lim[i]+']-lxir_210['+iagn_lim[i]+']')
    re = execute(lldet_052[i]+' = dblarr(nsrc)')
    re = execute(lllim_052[i]+' = dblarr(nsrc)')
    re = execute(lldet_052[i]+'['+iagn_det[i]+'] = '+lx_052[i]+'['+iagn_det[i]+']-lxir_052['+iagn_det[i]+']')
    re = execute(lllim_052[i]+'['+iagn_lim[i]+'] = '+lxlim_052[i]+'['+iagn_lim[i]+']-lxir_052['+iagn_lim[i]+']')
endfor

sav_vars = [sav_vars,'LLDET_210','LLLIM_210',lldet_210,lllim_210,'LLDET_052','LLLIM_052',lldet_052,lllim_052]
sav_inds = [sav_inds]

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="xray_lum.sav"')


;;----------------------------------------------------------------------------------------
;; NH DISTRIBUTIONS -- RAW UNMODELED OBSCURATION
;;----------------------------------------------------------------------------------------
if keyword_set(run) then begin
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

;;----------------------------------------------------------------------------------------
;; Luminosity data subset by field
;;----------------------------------------------------------------------------------------
lebv = alog10(ebv)>(-2.5)
lebv = lebv + randomu(seed,n_elements(ebv))*0.05-0.05/2.
ebvdet = {NST:lebv[iagn_det_nst],XMM:lebv[iagn_det_xmm],CHA:lebv[iagn_det_cha]}
ebvlim = {NST:lebv[iagn_lim_nst],XMM:lebv[iagn_lim_xmm],CHA:lebv[iagn_lim_cha]}
ebvwac = {NST:lebv[where(iiagn_lim_nst and iiwac)],XMM:lebv[where(iiagn_lim_xmm and iiwac)],CHA:lebv[where(iiagn_lim_cha and iiwac)]}
lldet = {NST:lldet_nst_052[iagn_det_nst],XMM:lldet_xmm_052[iagn_det_xmm],CHA:lldet_cha_052[iagn_det_cha]}
lllim = {NST:lllim_nst_052[iagn_lim_nst],XMM:lllim_xmm_052[iagn_lim_xmm],CHA:lllim_cha_052[iagn_lim_cha]}
llwac = {NST:lllim_nst_052[where(iiagn_lim_nst and iiwac)],XMM:lllim_xmm_052[where(iiagn_lim_xmm and iiwac)],CHA:lllim_cha_052[where(iiagn_lim_cha and iiwac)]}
ebvsub = ['EBVDET','EBVLIM','EBVWAC']
llsub = ['LLDET','LLLIM','LLWAC']

sav_vars = [sav_vars,'EBVSUB',ebvsub,'LLSUB',llsub]
sav_inds = [sav_inds]

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="xray_lum.sav"')


;;----------------------------------------------------------------------------------------
;; PLOTS!!
;;----------------------------------------------------------------------------------------
if keyword_set(plt) then begin

    ;; X-ray flux limits
    e = {symb:'o',linestyle:' ',sym_filled:1,sym_size:0.5,sym_transparency:75, $
         xr:[1e3,1e6],xlog:1,yr:[1e-14,1e-10],ylog:1, $
         xtitle:'$t_{exp} [s]$',ytitle:'$log( F_{2-10 keV} / erg s^{-1} cm^{-2} )$', $
         font_name:'Times',buffer:1}
    p1 = plot(cat_exp_nst,cat_flx_nst,col='purple',_extra=e,NAME='NuSTAR Catalogs',/nodata)
    p2 = plot(texp_nst[iagn_lim_nst],fxir[iagn_lim_nst],col='dark orange',_extra=e,/ov,NAME='AGN X-ray non-det.')
    p3 = plot(expt_nst[iagn_det_nst],fx_nst[iagn_det_nst],col='dodger blue','S',sym_filled=1,/ov,NAME='AGN + X-ray det.')
    ivalid = where(texp_nst gt 0.)
    p5 = plot([min(texp_nst[ivalid],imin),max(cat_exp_nst,imax)],[flim_nst[ivalid[imin]],cat_lim_nst[imax]],'--',/ov)     ;; flux limit
    l = legend(target=[p1,p2,p3],sample_width=0,position=[1.02,0.95],/relative,/auto_text_color)
    p1.save,'flux_limit_nst.png'
    ;; XMM
    e = {symb:'o',linestyle:' ',sym_filled:1,sym_size:0.5,sym_transparency:75, $
         xr:[3e2,3e5],yr:[1e-22,1e-6],xlog:1,ylog:1, $
         xtitle:'$t_{exp} [s]$',ytitle:'$log( F_{2-10 keV} / erg s^{-1} cm^{-2} )$', $
         font_name:'Times',buffer:1}
    p1 = plot(cat_exp_xmm,cat_flx_xmm,col='purple',_extra=e,NAME='3XMM-DR8 Catalog')
    p2 = plot(texp_xmm[iagn_lim_xmm],fxir[iagn_lim_xmm],col='dark orange',_extra=e,/ov,NAME='AGN X-ray non-det.')
    p3 = plot(expt_xmm[iagn_det_xmm],fx_xmm[iagn_det_xmm],col='dodger blue',_extra=e,/ov,NAME='AGN + X-ray det.')
    ivalid = where(texp_xmm gt 0.)
    p5 = plot([min(texp_xmm[ivalid],imin),max(cat_exp_xmm,imax)],[flim_xmm[ivalid[imin]],cat_lim_xmm[imax]],'--',/ov)     ;; flux limit
    l = legend(target=[p1,p2,p3],sample_width=0,position=[1.02,0.95],/relative,/auto_text_color)
    p1.save,'flux_limit_xmm.png'
    ;; Chanrda
    e = {symb:'o',linestyle:' ',sym_filled:1,sym_size:0.5,sym_transparency:75, $
         xr:[5e2,1e7],yr:[1e-18,1e-6],xlog:1,ylog:1, $
         xtitle:'$t_{exp} [s]$',ytitle:'$log( F_{2-10 keV} / erg s^{-1} cm^{-2} )$', $
         font_name:'Times',buffer:1}
    p1 = plot(cat_exp_cha,cat_flx_cha,col='purple',_extra=e,NAME='CSC2 Catalog')
    p2 = plot(texp_cha[iagn_lim_cha],fxir[iagn_lim_cha],col='dark orange',_extra=e,/ov,NAME='AGN X-ray non-det.')
    p3 = plot(expt_cha[iagn_det_cha],fx_cha[iagn_det_cha],col='dodger blue',_extra=e,/ov,NAME='AGN + X-ray det.')
    ivalid = where(texp_cha gt 0.)
    p5 = plot([min(texp_cha[ivalid],imin),max(cat_exp_cha,imax)],[flim_cha[ivalid[imin]],cat_lim_cha[imax]],'--',/ov)     ;; flux limit
    l = legend(target=[p1,p2,p3],sample_width=0,position=[1.02,0.95],/relative,/auto_text_color)
    p1.save,'flux_limit_cha.png'






    e = {symbol:'o',sym_size:0.5,sym_filled:1,sym_transparency:50,linestyle:'', $
         ;xra:[5e-3,1e2],xlog:1,yra:[-6,4],ylog:0, $
         xra:[-3.5,2.5],xlog:0,yra:[-6,4],ylog:0, $
         xtitle:'$log E(B-V)_{AGN}$',ytitle:'$log( L_{0.5-2 keV} / L_{0.5-2 keV}(L_{IR}) )$', $
         buffer:0}
    !NULL = plot(findgen(2),findgen(2),/nodata,_extra=e)
    plim = plot(ebvlim.nst,lllim.nst,col='orange',/ov,_extra=e,name='AGN X-ray non-det.')
    !NULL = plot(ebvlim.xmm,lllim.xmm,col='orange',/ov,_extra=e)
    !NULL = plot(ebvlim.cha,lllim.cha,col='orange',/ov,_extra=e)
    pdx = plot(ebvdet.xmm,lldet.xmm,col='medium blue',/ov,_extra=e,name='AGN + XMM')
    pdc = plot(ebvdet.cha,lldet.cha,col='dodger blue',/ov,_extra=e,name='AGN + Chandra')
    pdn = plot(ebvdet.nst,lldet.nst,col='sky blue',/ov,_extra=e,name='AGN + NuSTAR')
    !NULL = plot(e.xra,[0,0],'--',/ov,_extra=e)
    pwac = plot(ebvwac.nst,llwac.nst,'or',sym_size=0.75,/ov,name='WISE AGN 17')
    !NULL = plot(ebvwac.xmm,llwac.xmm,'or',sym_size=0.75,/ov)
    !NULL = plot(ebvwac.cha,llwac.cha,'or',sym_size=0.75,/ov)
    leg = legend(target=[plim,pdn,pdx,pdc,pwac],position=[1.05,0.87],/auto_text_color,sample_width=0,font_size=8)

    xnhlines = rebin([[-3.3],[-2.8]],4,2)
    ynhlines = rebin([-0.264805,-1.03190,-1.52835,-2.06607],4,2)
    nhlabel = ['$N_H=10^{24}cm^{-2}$','$    =5x10^{24}cm^{-2}$','$    =10^{25}cm^{-2}$','$    =5x10^{25}cm^{-2}$']
    for i = 0,3 do p = plot(xnhlines[i,*],ynhlines[i,*],':',thick=4,/ov)
    e = {font_size:8,font_style:'bold',data:1}
    for i = 0,3 do t = text(xnhlines[i,0]-0.15,ynhlines[i,0]-0.4,nhlabel[i],_extra=e)

    ;; 2D Histograms
    ct = COLORTABLE(55)     ;; Table 55 == CB-Oranges)
    hn = HIST_2D(ebvlim.nst, lllim.nst, bin1=0.2, bin2=0.2)
    ;hn = CONGRID(hn, MAX(ebvlim.nst), MAX(lllim.nst))
    ;hn = BYTSCL(ALOG10(hn > 1))
    gn = IMAGE(hn, AXIS_STYLE=2, MARGIN=0.1, rgb_table=ct)    
    hx = HIST_2D(ebvlim.xmm, lllim.xmm, bin1=0.2, bin2=0.2)
    ;hx = CONGRID(hx, MAX(ebvlim.xmm), MAX(lllim.xmm))
    ;hx = BYTSCL(ALOG10(hx > 1))
    gx = IMAGE(hx, AXIS_STYLE=2, MARGIN=0.1, rgb_table=ct)
    hc = HIST_2D(ebvlim.cha, lllim.cha, bin1=0.2, bin2=0.2)
    ;hc = CONGRID(hc, MAX(ebvlim.cha), MAX(lllim.cha))
    ;hc = BYTSCL(ALOG10(hc > 1))
    gc = IMAGE(hc, AXIS_STYLE=2, MARGIN=0.1, rgb_table=ct)
    
endif














END










