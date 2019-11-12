PRO fit_catalog_flim, MULTI_SN = multi_sn, $
                      PLT = plt


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

if keyword_set(multi_sn) then begin
    ;; S/N ³ 3
    cat_exp3 = 'CAT_EXP3'+xfield
    cat_flx3 = 'CAT_FLX3'+xfield
    cat_err3 = 'CAT_ERR3'+xfield
    cat_lim3 = 'CAT_LIM3'+xfield

    for f = 0,nfield-1 do begin
        cat_fld = (strsplit(xfield[f],'_',/extract))[0]
        re = execute('ncat = n_elements('+cat_fld+')')
        re = execute(cat_exp3[f]+' = dblarr(ncat)')
        re = execute(cat_flx3[f]+' = dblarr(ncat)')
        re = execute(cat_err3[f]+' = dblarr(ncat)')
        conv = xconv.(f)
        iorder = sort(abs((conv-conv[ixband[f]])/conv[ixband[f]]))
        for b = 0,nxband[f]-1 do begin
            re = execute('iuse = where('+cat_fld+'.'+tt.(f)[iorder[b]]+' gt 0. and '+cat_fld+'.'+ff.(f)[iorder[b]]+' gt 0. and '+cat_exp3[f]+' eq 0.,nuse)')
            if (nuse gt 0.) then begin
                re = execute(cat_exp3[f]+'[iuse] = '+cat_fld+'[iuse].'+tt.(f)[iorder[b]])
                re = execute(cat_flx3[f]+'[iuse] = '+cat_fld+'[iuse].'+ff.(f)[iorder[b]]+' * conv[iorder[b]]')
                re = execute(cat_err3[f]+'[iuse] = '+cat_fld+'[iuse].'+ee.(f)[iorder[b]]+' * conv[iorder[b]]')
            endif
        endfor
        re = execute('isort = sort('+cat_exp3[f]+')')
        re = execute(cat_exp3[f]+' = '+cat_exp3[f]+'[isort]')
        re = execute(cat_flx3[f]+' = '+cat_flx3[f]+'[isort]')
        re = execute(cat_err3[f]+' = '+cat_err3[f]+'[isort]')
        re = execute('cat_sn3 = '+cat_flx3[f]+'/'+cat_err3[f])
        re = execute('icat_sn3 = where(cat_sn3 ge 3. and '+cat_exp3[f]+' gt 0.,ng)')
        if (ng eq 0.) then stop
        re = execute(cat_exp3[f]+' = '+cat_exp3[f]+'[icat_sn3]')
        re = execute(cat_flx3[f]+' = '+cat_flx3[f]+'[icat_sn3]')
        re = execute(cat_err3[f]+' = '+cat_err3[f]+'[icat_sn3]')
    endfor

    for f = 0,nfield-1 do re = execute(cat_lim3[f]+' = xray_flux_limit('+cat_exp3[f]+','+cat_flx3[f]+','+cat_err3[f]+')')

    sav_vars = [sav_vars,cat_exp3,cat_flx3,cat_err3,cat_lim3]
    sav_inds = [sav_inds]

    ;; S/N ³ 5
    cat_exp5 = 'CAT_EXP5'+xfield
    cat_flx5 = 'CAT_FLX5'+xfield
    cat_err5 = 'CAT_ERR5'+xfield
    cat_lim5 = 'CAT_LIM5'+xfield

    for f = 0,nfield-1 do begin
        cat_fld = (strsplit(xfield[f],'_',/extract))[0]
        re = execute('ncat = n_elements('+cat_fld+')')
        re = execute(cat_exp5[f]+' = dblarr(ncat)')
        re = execute(cat_flx5[f]+' = dblarr(ncat)')
        re = execute(cat_err5[f]+' = dblarr(ncat)')
        conv = xconv.(f)
        iorder = sort(abs((conv-conv[ixband[f]])/conv[ixband[f]]))
        for b = 0,nxband[f]-1 do begin
            re = execute('iuse = where('+cat_fld+'.'+tt.(f)[iorder[b]]+' gt 0. and '+cat_fld+'.'+ff.(f)[iorder[b]]+' gt 0. and '+cat_exp5[f]+' eq 0.,nuse)')
            if (nuse gt 0.) then begin
                re = execute(cat_exp5[f]+'[iuse] = '+cat_fld+'[iuse].'+tt.(f)[iorder[b]])
                re = execute(cat_flx5[f]+'[iuse] = '+cat_fld+'[iuse].'+ff.(f)[iorder[b]]+' * conv[iorder[b]]')
                re = execute(cat_err5[f]+'[iuse] = '+cat_fld+'[iuse].'+ee.(f)[iorder[b]]+' * conv[iorder[b]]')
            endif
        endfor
        re = execute('isort = sort('+cat_exp5[f]+')')
        re = execute(cat_exp5[f]+' = '+cat_exp5[f]+'[isort]')
        re = execute(cat_flx5[f]+' = '+cat_flx5[f]+'[isort]')
        re = execute(cat_err5[f]+' = '+cat_err5[f]+'[isort]')
        re = execute('cat_sn5 = '+cat_flx5[f]+'/'+cat_err5[f])
        re = execute('icat_sn5 = where(cat_sn5 ge 5. and '+cat_exp5[f]+' gt 0.,ng)')
        if (ng eq 0.) then stop
        re = execute(cat_exp5[f]+' = '+cat_exp5[f]+'[icat_sn5]')
        re = execute(cat_flx5[f]+' = '+cat_flx5[f]+'[icat_sn5]')
        re = execute(cat_err5[f]+' = '+cat_err5[f]+'[icat_sn5]')
    endfor

    for f = 0,nfield-1 do re = execute(cat_lim5[f]+' = xray_flux_limit('+cat_exp5[f]+','+cat_flx5[f]+','+cat_err5[f]+')')

    sav_vars = [sav_vars,cat_exp5,cat_flx5,cat_err5,cat_lim5]
    sav_inds = [sav_inds]
endif

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




