PRO fit_catalog_xlim, PLT = plt


common _inf_cha
common _inf_xmm
common _inf_nst
common _det_cha
common _det_xmm
common _det_nst
common _softx


;;----------------------------------------------------------------------------------------
;; X-RAY flux_inf LIMITS
;;----------------------------------------------------------------------------------------
;clean_cha = clean_source_chandra(cha)
;clean_xmm = clean_source_xmm(xmm)
;clean_nst = clean_source_nustar(nst)
;;----------------------------------------------------------------------------------------
;; Energy band 2-10 keV
cat_exp_210 = 'CAT_EXP'+xfield+'_210'
cat_flx_210 = 'CAT_FLX'+xfield+'_210'
cat_err_210 = 'CAT_ERR'+xfield+'_210'
cat_lim_210 = 'CAT_LIM'+xfield+'_210'

for f = 0,nfield-1 do begin
    cat_fld = (strsplit(xfield[f],'_',/extract))[0]
    re = execute('ncat = n_elements('+cat_fld+')')
    re = execute(cat_exp_210[f]+' = dblarr(ncat)')
    re = execute(cat_flx_210[f]+' = dblarr(ncat)')
    re = execute(cat_err_210[f]+' = dblarr(ncat)')
    conv = xconv.(f)
    iorder = sort(abs((conv-conv[ixband[f]])/conv[ixband[f]]))
    for b = 0,nxband[f]-1 do begin
        re = execute('iuse = where('+cat_fld+'.'+tt.(f)[iorder[b]]+' gt 0. and '+cat_fld+'.'+ff.(f)[iorder[b]]+' gt 0. and '+cat_exp_210[f]+' eq 0.,nuse)')
        if (nuse gt 0.) then begin
            re = execute(cat_exp_210[f]+'[iuse] = '+cat_fld+'[iuse].'+tt.(f)[iorder[b]])
            re = execute(cat_flx_210[f]+'[iuse] = '+cat_fld+'[iuse].'+ff.(f)[iorder[b]]+' * conv[iorder[b]]')
            re = execute(cat_err_210[f]+'[iuse] = '+cat_fld+'[iuse].'+ee.(f)[iorder[b]]+' * conv[iorder[b]]')
        endif
    endfor
    re = execute('isort = sort('+cat_exp_210[f]+')')
    re = execute(cat_exp_210[f]+' = '+cat_exp_210[f]+'[isort]')
    re = execute(cat_flx_210[f]+' = '+cat_flx_210[f]+'[isort]')
    re = execute(cat_err_210[f]+' = '+cat_err_210[f]+'[isort]')
    re = execute('cat_sn = '+cat_flx_210[f]+'/'+cat_err_210[f])
    re = execute('icat_sn = where(cat_sn ge 3. and '+cat_exp_210[f]+' gt 0.,ng)')
    if (ng eq 0.) then stop
    re = execute(cat_exp_210[f]+' = '+cat_exp_210[f]+'[icat_sn]')
    re = execute(cat_flx_210[f]+' = '+cat_flx_210[f]+'[icat_sn]')
    re = execute(cat_err_210[f]+' = '+cat_err_210[f]+'[icat_sn]')
endfor

for f = 0,nfield-1 do re = execute(cat_lim_210[f]+' = xray_flim('+cat_exp_210[f]+','+cat_flx_210[f]+','+cat_err_210[f]+')')

sav_vars = ['CAT_EXP_210','CAT_FLX_210','CAT_ERR_210','CAT_LIM_210', $
            cat_exp_210,cat_flx_210,cat_err_210,cat_lim_210]
sav_inds = []

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="catalog_flux_limits.sav"')


;; plot flux limits
if keyword_set(plt) then begin
    e = {color:'light grey',transparency:0, $
         xra:[2,7],yra:[-16,-10],xlog:0,ylog:0, $
         aspect_ratio:1,dimension:[750,330], $
         buffer:1}
    current = 0
    for i = 0,nfield-1 do begin
        re = execute('p = plot(alog10('+cat_exp_210[i]+'),alog10('+cat_flx_210[i]+'),".",_extra=e,layout=[3,1,i+1],current=current)')
        re = execute('p = plot(alog10('+cat_exp_210[i]+'),alog10('+cat_lim_210[i]+'),"--r",thick=2,/ov)')
        p.xtitle = ('$t_{exp} '+['Chandra','XMM','NuSTAR']+' [ks]$')[i]
        p.ytitle = ('$log '+['F_{X,2-10keV}','F_{X,2-10keV}','F_{X,2-10keV}']+' [erg s^{-1} cm^{-2}]$')[i]
        if (current eq 0) then current = 1
    endfor
    p.save,'plot_flux_limit.png'
endif


END

;; calculate catalog flux limit
;root = [0.5,5./6.,1.]
;cn = [3e-13,1.5e-13,3.5e-13]
;for i = 0,nfield-1 do re = execute(cat_lim_210[i]+' = xray_flim('+cat_exp_210[i]+','+cat_flx_210[i]+','+cat_err_210[i]+',root=root[i],cn=cn[i])')
;root1 = dblarr(nfield)
;root2 = dblarr(nfield)
;cn = dblarr(nfield)
;for f = 0,nfield-1 do begin
;    re = execute(cat_lim_210[f]+' = xray_flim('+cat_exp_210[f]+','+cat_flx_210[f]+','+cat_err_210[f]+',phot_k=exp1,back_k=exp2,nrml=nrml)')
;    root1[f] = exp1
;    root2[f] = exp2
;    cn[f] = nrml
;endfor



;;----------------------------------------------------------------------------------------
;; Energy band 0.5-2 keV
;used_exp = tt[ixband_052]
;used_flx = ff[ixband_052]
;used_err = ee[ixband_052]
;used_cnv = xray_cnv_052
;
;cat_exp_052 = 'CAT_EXP'+xfield+'_052'
;cat_flx_052 = 'CAT_FLX'+xfield+'_052'
;cat_err_052 = 'CAT_ERR'+xfield+'_052'
;cat_lim_052 = 'CAT_LIM'+xfield+'_052'
;
;for i = 0,nfield-1 do begin
;    cat_fld = (strsplit(xfield[i],'_',/extract))[0]
;    raw_exp = cat_fld+'[sort('+cat_fld+'.'+used_exp[i]+')].'+used_exp[i]
;    raw_flx = cat_fld+'[sort('+cat_fld+'.'+used_exp[i]+')].'+used_flx[i]
;    raw_err = cat_fld+'[sort('+cat_fld+'.'+used_exp[i]+')].'+used_err[i]
;    re = execute(cat_exp_052[i]+' = '+raw_exp)
;    re = execute(cat_flx_052[i]+' = '+raw_flx)
;    re = execute(cat_err_052[i]+' = '+raw_err)
;    re = execute('cat_sn = '+cat_flx_052[i]+'/'+cat_err_052[i])
;    re = execute('cat_ind = where('+cat_exp_052[i]+' gt 0. and cat_sn ge 3.,ng)')
;    if (ng eq 0.) then stop
;    re = execute(cat_exp_052[i]+' = '+cat_exp_052[i]+'[cat_ind]')
;    re = execute(cat_flx_052[i]+' = '+cat_flx_052[i]+'[cat_ind] * '+used_cnv[i])
;    re = execute(cat_err_052[i]+' = '+cat_err_052[i]+'[cat_ind] * '+used_cnv[i])
;endfor
;
;;; calculate catalog flux limit
;;dex = [3.,3.,3.]
;;root = [0.5,5./6.,1.]
;;cn = [9e-14,3.6e-14,8.8e-14]
;;for i = 0,nfield-1 do re = execute(cat_lim_052[i]+' = xray_flim('+cat_exp_052[i]+','+cat_flx_052[i]+','+cat_err_052[i]+',root=root[i],cn=cn[i])')
;
;root1 = dblarr(nfield)
;root2 = dblarr(nfield)
;cn = dblarr(nfield)
;for i = 0,nfield-1 do begin
;    re = execute(cat_lim_052[i]+' = xray_flim('+cat_exp_052[i]+','+cat_flx_052[i]+','+cat_err_052[i]+',phot_k=exp1,back_k=exp2,nrml=nrml)')
;    root1[i] = exp1
;    root2[i] = exp2
;    cn[i] = nrml
;endfor
;
;sav_vars = [sav_vars,'CAT_EXP_052','CAT_FLX_052','CAT_ERR_052','CAT_LIM_052', $
;                     cat_exp_052,cat_flx_052,cat_err_052,cat_lim_052]
;sav_inds = [sav_inds]