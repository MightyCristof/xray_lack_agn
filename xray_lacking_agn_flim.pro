PRO xray_lacking_agn_flim


common _inf_nst
common _inf_xmm
common _inf_cha
common _det_nst
common _det_xmm
common _det_cha
common _soft210
common _soft052


;;----------------------------------------------------------------------------------------
;; X-RAY flux_inf LIMITS
;;----------------------------------------------------------------------------------------

;; Energy band 2-10 keV
used_exp = ['SEXP','EP_ONTIME','ACIS_TIME']
used_flx = ['SBF','SC_EP_4_FLUX','FLUX_POWLAW_APER90_H']
used_err = ['E_SBF','SC_EP_4_FLUX_ERR','FLUX_POWLAW_APER90_LOLIM_H']
used_cnv = ['OUT210_NST_S','OUT210_XMM_2','OUT210_CHA_S']

cat_exp_210 = 'CAT'+xfield+'_EXP_210'
cat_flx_210 = 'CAT'+xfield+'_FLX_210'
cat_err_210 = 'CAT'+xfield+'_ERR_210'
cat_lim_210 = 'CAT'+xfield+'_LIM_210'

for i = 0,nfield-1 do begin
    cat_fld = (strsplit(xfield[i],'_',/extract))[0]
    raw_exp = cat_fld+'.'+used_exp[i]
    raw_flx = cat_fld+'.'+used_flx[i]
    if (cat_fld eq 'CHA') then raw_err = cat_fld+'.'+used_flx[i]+' - '+cat_fld+'.'+used_err[i] else $
                               raw_err = cat_fld+'.'+used_err[i]
    re = execute(cat_exp_210[i]+' = '+raw_exp)
    re = execute(cat_flx_210[i]+' = '+raw_flx)
    re = execute(cat_err_210[i]+' = '+raw_err)
    re = execute('raw_sn = '+cat_flx_210[i]+'/'+cat_err_210[i])
    re = execute('cat_ind = where('+cat_exp_210[i]+' gt 0. and raw_sn gt 3.,ng)')
    if (ng lt 1.) then stop
    re = execute(cat_exp_210[i]+' = '+cat_exp_210[i]+'[cat_ind]')
    re = execute(cat_flx_210[i]+' = '+cat_flx_210[i]+'[cat_ind] * '+used_cnv[i])
    re = execute(cat_err_210[i]+' = '+cat_err_210[i]+'[cat_ind] * '+used_cnv[i])
endfor

;; calculate catalog flux limit
dex = [4.,4.,4.]
root = [0.5,2./3.,1.]
cn = [1.,1.,1.]
for i = 0,nfield-1 do re = execute(cat_lim_210[i]+' = xray_flim('+cat_exp_210[i]+','+cat_flx_210[i]+','+cat_err_210[i]+',dex=dex[i],root=root[i],cn=cn[i])')


sav_vars = ['CAT_EXP_210','CAT_FLX_210','CAT_ERR_210','CAT_LIM_210', $
            cat_exp_210,cat_flx_210,cat_err_210,cat_lim_210]
sav_inds = []


;; Energy band 0.5-2 keV
used_exp = ['SEXP','EP_ONTIME','ACIS_TIME']
used_flx = ['SBF','SC_EP_2_FLUX','FLUX_POWLAW_APER90_S']
used_err = ['E_SBF','SC_EP_2_FLUX_ERR','FLUX_POWLAW_APER90_LOLIM_S']
used_cnv = ['OUT052_NST_S','OUT052_XMM_2','OUT052_CHA_S']

cat_exp_052 = 'CAT'+xfield+'_EXP_052'
cat_flx_052 = 'CAT'+xfield+'_FLX_052'
cat_err_052 = 'CAT'+xfield+'_ERR_052'
cat_lim_052 = 'CAT'+xfield+'_LIM_052'

for i = 0,nfield-1 do begin
    cat_fld = (strsplit(xfield[i],'_',/extract))[0]
    raw_exp = cat_fld+'.'+used_exp[i]
    raw_flx = cat_fld+'.'+used_flx[i]
    if (cat_fld eq 'CHA') then raw_err = cat_fld+'.'+used_flx[i]+' - '+cat_fld+'.'+used_err[i] else $
                               raw_err = cat_fld+'.'+used_err[i]
    re = execute(cat_exp_052[i]+' = '+raw_exp)
    re = execute(cat_flx_052[i]+' = '+raw_flx)
    re = execute(cat_err_052[i]+' = '+raw_err)
    re = execute('raw_sn = '+cat_flx_052[i]+'/'+cat_err_210[i])
    re = execute('cat_ind = where('+cat_exp_052[i]+' gt 0. and raw_sn gt 3.,ng)')
    if (ng lt 1.) then stop
    re = execute(cat_exp_052[i]+' = '+cat_exp_052[i]+'[cat_ind]')
    re = execute(cat_flx_052[i]+' = '+cat_flx_052[i]+'[cat_ind] * '+used_cnv[i])
    re = execute(cat_err_052[i]+' = '+cat_err_052[i]+'[cat_ind] * '+used_cnv[i])
endfor

;; calculate catalog flux limit
dex = [4.,4.,4.]
root = [0.5,2./3.,1.]
cn = [1.,1.,1.]
for i = 0,nfield-1 do re = execute(cat_lim_210[i]+' = xray_flim('+cat_exp_210[i]+','+cat_flx_210[i]+','+cat_err_210[i]+',dex=dex[i],root=root[i],cn=cn[i])')


sav_vars = ['CAT_EXP_210','CAT_FLX_210','CAT_ERR_210','CAT_LIM_210', $
            cat_exp_210,cat_flx_210,cat_err_210,cat_lim_210]
sav_inds = []



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







;; plot flux limits
if keyword_set(plt) then begin
    xtitle = '$'+['t_{exp} NuSTAR','t_{exp} XMM','t_{exp} Chandra']+'$'
    ytitle = '$'+['F_{X,2-10keV}','F_{X,2-10keV}','F_{X,2-10keV}']+'$'
    for i = 0,nfield-1 do begin
        ;re = execute('p = errorplot('+cat_exp_210[i]+','+cat_flx_210[i]+','+cat_err_210[i]+',".",linestyle="",/xlog,/ylog,xtitle=xtitle[i],ytitle=ytitle[i])')
        re = execute('p = plot('+cat_exp_210[i]+','+cat_flx_210[i]+',".",/xlog,/ylog,xtitle=xtitle[i],ytitle=ytitle[i])')
        re = execute('p = plot('+cat_exp_210[i]+','+cat_lim_210[i]+',"--r",/ov)')
    endfor
endif














END




