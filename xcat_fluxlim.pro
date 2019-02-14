PRO xcat_fluxlim, PLT = plt


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

;;----------------------------------------------------------------------------------------
;; Energy band 2-10 keV
used_exp = tt[ixband_210]
used_flx = ff[ixband_210]
used_err = ee[ixband_210]
used_cnv = xray_cnv_210

cat_exp_210 = 'CAT_EXP'+xfield+'_210'
cat_flx_210 = 'CAT_FLX'+xfield+'_210'
cat_err_210 = 'CAT_ERR'+xfield+'_210'
cat_lim_210 = 'CAT_LIM'+xfield+'_210'

for i = 0,nfield-1 do begin
    cat_fld = (strsplit(xfield[i],'_',/extract))[0]
    raw_exp = cat_fld+'.'+used_exp[i]
    raw_flx = cat_fld+'.'+used_flx[i]
    if (cat_fld eq 'CHA') then raw_err = cat_fld+'.'+used_flx[i]+' - '+cat_fld+'.'+used_err[i] else $
                               raw_err = cat_fld+'.'+used_err[i]
    re = execute(cat_exp_210[i]+' = '+raw_exp)
    re = execute(cat_flx_210[i]+' = '+raw_flx)
    re = execute(cat_err_210[i]+' = '+raw_err)
    re = execute('cat_sn = '+cat_flx_210[i]+'/'+cat_err_210[i])
    re = execute('cat_ind = where('+cat_exp_210[i]+' gt 0. and cat_sn gt 3.,ng)')
    if (ng lt 1.) then stop
    re = execute(cat_exp_210[i]+' = '+cat_exp_210[i]+'[cat_ind]')
    re = execute(cat_flx_210[i]+' = '+cat_flx_210[i]+'[cat_ind] * '+used_cnv[i])
    re = execute(cat_err_210[i]+' = '+cat_err_210[i]+'[cat_ind] * '+used_cnv[i])
endfor

;; calculate catalog flux limit
dex = [4.,4.,4.]
root = [0.5,5./6.,1.]
cn = [1e-13,2.5e-14,4e-14]
for i = 0,nfield-1 do re = execute(cat_lim_210[i]+' = xray_flim('+cat_exp_210[i]+','+cat_flx_210[i]+','+cat_err_210[i]+',dex=dex[i],root=root[i],cn=cn[i])')
    

sav_vars = ['CAT_EXP_210','CAT_FLX_210','CAT_ERR_210','CAT_LIM_210', $
            cat_exp_210,cat_flx_210,cat_err_210,cat_lim_210]
sav_inds = []


;;----------------------------------------------------------------------------------------
;; Energy band 0.5-2 keV
used_exp = tt[ixband_052]
used_flx = ff[ixband_052]
used_err = ee[ixband_052]
used_cnv = xray_cnv_052

cat_exp_052 = 'CAT_EXP'+xfield+'_052'
cat_flx_052 = 'CAT_FLX'+xfield+'_052'
cat_err_052 = 'CAT_ERR'+xfield+'_052'
cat_lim_052 = 'CAT_LIM'+xfield+'_052'

for i = 0,nfield-1 do begin
    cat_fld = (strsplit(xfield[i],'_',/extract))[0]
    raw_exp = cat_fld+'.'+used_exp[i]
    raw_flx = cat_fld+'.'+used_flx[i]
    if (cat_fld eq 'CHA') then raw_err = cat_fld+'.'+used_flx[i]+' - '+cat_fld+'.'+used_err[i] else $
                               raw_err = cat_fld+'.'+used_err[i]
    re = execute(cat_exp_052[i]+' = '+raw_exp)
    re = execute(cat_flx_052[i]+' = '+raw_flx)
    re = execute(cat_err_052[i]+' = '+raw_err)
    re = execute('cat_sn = '+cat_flx_052[i]+'/'+cat_err_052[i])
    re = execute('cat_ind = where('+cat_exp_052[i]+' gt 0. and cat_sn gt 3.,ng)')
    if (ng lt 1.) then stop
    re = execute(cat_exp_052[i]+' = '+cat_exp_052[i]+'[cat_ind]')
    re = execute(cat_flx_052[i]+' = '+cat_flx_052[i]+'[cat_ind] * '+used_cnv[i])
    re = execute(cat_err_052[i]+' = '+cat_err_052[i]+'[cat_ind] * '+used_cnv[i])
endfor

;; calculate catalog flux limit
dex = [4.,4.,4.]
root = [0.5,5./6.,1.]
cn = [3e-14,6e-15,1e-14]
for i = 0,nfield-1 do re = execute(cat_lim_052[i]+' = xray_flim('+cat_exp_052[i]+','+cat_flx_052[i]+','+cat_err_052[i]+',dex=dex[i],root=root[i],cn=cn[i])')

sav_vars = [sav_vars,'CAT_EXP_052','CAT_FLX_052','CAT_ERR_052','CAT_LIM_052', $
                     cat_exp_052,cat_flx_052,cat_err_052,cat_lim_052]
sav_inds = [sav_inds]

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="xray_flx_lim.sav"')


;; plot flux limits
if keyword_set(plt) then begin
    e = {color:'dodger blue', $
         ;xra:[5e2,5e6],yra:[1e-16,1e-10],xlog:1,ylog:1, $
         xra:[2.5,6.5],yra:[-16,-10],xlog:0,ylog:0, $
         aspect_ratio:1,dimension:[1250,550], $
         buffer:0}
    ;; 2-10 keV
    xtitle = '$'+['t_{exp} NuSTAR','t_{exp} XMM','t_{exp} Chandra']+' [ks]$'
    ytitle = '$log '+['F_{X,2-10keV}','F_{X,2-10keV}','F_{X,2-10keV}']+' [erg s^{-1} cm^{-2}]$'
    title = '$X-ray: '+['2-10','0.5-2']+' keV$'
    current = 0
    for i = 0,nfield-1 do begin
        ;re = execute('p = plot('+cat_exp_210[i]+','+cat_flx_210[i]+',".",_extra=e,layout=[3,1,i+1],current=current)')
        ;re = execute('p = plot('+cat_exp_210[i]+','+cat_lim_210[i]+',"--r",/ov)')
        re = execute('p = plot(alog10('+cat_exp_210[i]+'),alog10('+cat_flx_210[i]+'),".",_extra=e,layout=[3,1,i+1],current=current)')
        re = execute('p = plot(alog10('+cat_exp_210[i]+'),alog10('+cat_lim_210[i]+'),"--r",/ov)')
        p.xtitle = xtitle[i]
        p.ytitle = ytitle[i]
        if (current eq 0) then current = 1
    endfor
    t = text(0.5,0.9,'$2-10 keV$',font_style='bold',alignment=0.5,/normal)
    p.save,'flux_limit_soft210.png'
    ;; 0.5-2 keV
    xtitle = '$'+['t_{exp} NuSTAR','t_{exp} XMM','t_{exp} Chandra']+'$'
    ytitle = '$'+['F_{X,0.5-2keV}','F_{X,0.5-2keV}','F_{X,0.5-2keV}']+'$'
    current = 0
    for i = 0,nfield-1 do begin
        re = execute('p = plot('+cat_exp_052[i]+','+cat_flx_052[i]+',".",_extra=e,layout=[3,1,i+1],current=current)')
        re = execute('p = plot('+cat_exp_052[i]+','+cat_lim_052[i]+',"--r",/ov)')
        p.xtitle = xtitle[i]
        p.ytitle = ytitle[i]
        if (current eq 0) then current = 1
    endfor
    t = text(0.5,0.9,'$0.5-2 keV$',font_style='bold',alignment=0.5,/normal)
    p.save,'flux_limit_soft052.png'
endif


END




