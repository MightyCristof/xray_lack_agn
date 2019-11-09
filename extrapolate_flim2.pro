FUNCTION extrapolate_flim2, cat_lim, $
                           cat_exp, $
                           src_exp
                                 

nexp = n_elements(src_exp)
iexp = where(src_exp gt 0.)

lcat_lim = alog10(cat_lim[sort(cat_exp)])
lcat_exp = alog10(cat_exp[sort(cat_exp)])
lsrc_exp = alog10(src_exp[iexp])

src_lim = dblarr(nexp)
src_lim[iexp] = interpol(lcat_lim,lcat_exp,lsrc_exp)


return, src_lim


END



;cat_exp = cat_exp_cha
;cat_flx = cat_flx_cha

;lcat_exp = alog10(cat_exp[sort(cat_exp)])
;lcat_flx = alog10(cat_flx[sort(cat_exp)])
;fitc = poly_fit(lcat_exp,lcat_flx,degr,yfit=yfit)
;p = plot(lcat_exp,yfit,'-r',/ov)
;p = plot(lcat_exp,yfit-1.,'-r',/ov)