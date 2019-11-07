PRO fit_catalog_flim, PLT = plt


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
cat_exp = 'CAT_EXP'+xfield
cat_flx = 'CAT_FLX'+xfield
cat_err = 'CAT_ERR'+xfield
cat_lim = 'CAT_LIM'+xfield

for f = 0,nfield-1 do begin
    cat_fld = (strsplit(xfield[f],'_',/extract))[0]
    re = execute('ncat = n_elements('+cat_fld+')')
    re = execute(cat_exp[f]+' = dblarr(ncat)')
    re = execute(cat_flx[f]+' = dblarr(ncat)')
    re = execute(cat_err[f]+' = dblarr(ncat)')
    conv = xconv.(f)
    iorder = sort(abs((conv-conv[ixband[f]])/conv[ixband[f]]))
    for b = 0,nxband[f]-1 do begin
        re = execute('iuse = where('+cat_fld+'.'+tt.(f)[iorder[b]]+' gt 0. and '+cat_fld+'.'+ff.(f)[iorder[b]]+' gt 0. and '+cat_exp[f]+' eq 0.,nuse)')
        if (nuse gt 0.) then begin
            re = execute(cat_exp[f]+'[iuse] = '+cat_fld+'[iuse].'+tt.(f)[iorder[b]])
            re = execute(cat_flx[f]+'[iuse] = '+cat_fld+'[iuse].'+ff.(f)[iorder[b]]+' * conv[iorder[b]]')
            re = execute(cat_err[f]+'[iuse] = '+cat_fld+'[iuse].'+ee.(f)[iorder[b]]+' * conv[iorder[b]]')
        endif
    endfor
    re = execute('isort = sort('+cat_exp[f]+')')
    re = execute(cat_exp[f]+' = '+cat_exp[f]+'[isort]')
    re = execute(cat_flx[f]+' = '+cat_flx[f]+'[isort]')
    re = execute(cat_err[f]+' = '+cat_err[f]+'[isort]')
    re = execute('cat_sn = '+cat_flx[f]+'/'+cat_err[f])
    re = execute('icat_sn = where(cat_sn ge 2. and '+cat_exp[f]+' gt 0.,ng)')
    if (ng eq 0.) then stop
    re = execute(cat_exp[f]+' = '+cat_exp[f]+'[icat_sn]')
    re = execute(cat_flx[f]+' = '+cat_flx[f]+'[icat_sn]')
    re = execute(cat_err[f]+' = '+cat_err[f]+'[icat_sn]')
endfor

for f = 0,nfield-1 do re = execute(cat_lim[f]+' = xray_flux_limit('+cat_exp[f]+','+cat_flx[f]+','+cat_err[f]+')')

sav_vars = [cat_exp,cat_flx,cat_err,cat_lim]
sav_inds = []

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="catalog_flux_limits.sav"')


;; plot flux limits
if keyword_set(plt) then begin
    e = {color:'light grey',transparency:0, $
         xra:[2.5,6.5],yra:[-16,-12], $
         xlog:0,ylog:0, $
         aspect_ratio:0,dimension:[750,330], $
         buffer:0}
    ;xra = [[2.,7.],[2.,6.],[3.,7.]]
    ;yra = [[-20.,-10.],[-20.,-10.],[-20.,-10.]]
    ;xr=xra[*,i],yr=yra[*,i],
    current = 0
    for i = 0,nfield-1 do begin
        re = execute('len = n_elements('+cat_exp[i]+')')
        re = execute('irand = round(randomu(seed,len/4.)*n_elements('+cat_exp[i]+'))')
        irand = irand[sort(irand)]
        re = execute('p = plot(alog10('+cat_exp[i]+'[irand]),alog10('+cat_flx[i]+'[irand]),".",_extra=e,layout=[3,1,i+1],current=current)')
        re = execute('p = plot(alog10('+cat_exp[i]+'[irand]),alog10('+cat_lim[i]+'[irand]),"--r",thick=2,/ov)')
        p.xtitle = ('$t_{exp} '+['Chandra','XMM','NuSTAR']+' [ks]$')[i]
        p.ytitle = ('$log '+['F_{X,2-10keV}','F_{X,2-10keV}','F_{X,2-10keV}']+' [erg s^{-1} cm^{-2}]$')[i]
        if (current eq 0) then current = 1
    endfor
    ;p.save,'plot_flux_limit.png'
endif


END




