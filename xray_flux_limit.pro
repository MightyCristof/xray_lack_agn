FUNCTION xray_flux_limit, xexp, $
                          xflx, $
                          xerr, $
                          nbins


;; convert to log-log
isort = sort(xexp)
lexp = alog10(xexp[isort])
lflx = alog10(xflx[isort])
lerr = sqrt((xerr[isort]/(xflx[isort]*alog(10)))^2)
;; number of sources in catalog
ncat = n_elements(lexp)
;binsz = floor(ncat*0.05)
mmlexp = minmax(lexp)
;; DIVIDE DATA IN HALVES
;middle = total(mmlexp)/2.
;binsz = (middle-mmlexp[0])/nbins
;; DIVIDE DATA IN THIRDS
thirds = diff(mmlexp)/3.
binsz = thirds/nbins
;fourths = diff(mmlexp)/4.
;binsz = fourths/nbins
frac = 0.05

;; PHOTON-LIMITED CASE
bins1 = []
mins1 = []
errs1 = []
for i = 0,nbins-1 do begin
    binra = [mmlexp[0]+binsz*i,mmlexp[0]+binsz*(i+1)]
    ibin = where(lexp ge binra[0] and lexp lt binra[1],binct)
    ;print, binct
    if (binct lt 10 < ncat*0.02) then continue
    isrt = sort(lflx[ibin])
    arrlen = binct*frac < 500
    bins1 = [bins1,median(lexp[ibin[isrt[0:arrlen]]])]
    mins1 = [mins1,median(lflx[ibin[isrt[0:arrlen]]])]
    ;mins1 = [mins1,median(lflx[ibin[isrt[0:arrlen]]]-lerr[ibin[isrt[0:arrlen]]])]
endfor

;; CONNECTING REGION TO SMOOTH TRANSITION
bins3 = []
mins3 = []
errs3 = []
for i = nbins,2*nbins-1 do begin
    binra = [mmlexp[0]+binsz*i,mmlexp[0]+binsz*(i+1)]
    ibin = where(lexp ge binra[0] and lexp lt binra[1],binct)
    ;print, binct
    if (binct lt 10 < ncat*0.02) then continue
    isrt = sort(lflx[ibin])
    arrlen = binct*frac < 500
    bins3 = [bins3,median(lexp[ibin[isrt[0:arrlen]]])]
    mins3 = [mins3,median(lflx[ibin[isrt[0:arrlen]]])]
    ;mins3 = [mins3,median(lflx[ibin[isrt[0:arrlen]]]-lerr[ibin[isrt[0:arrlen]]])]
endfor

;; BACKGROUND-LIMITED CASE
bins2 = []
mins2 = []
errs2 = []
for i = 2*nbins,3*nbins-1 do begin
    binra = [mmlexp[0]+binsz*i,mmlexp[0]+binsz*(i+1)]
    ibin = where(lexp ge binra[0] and lexp lt binra[1],binct)
    ;print, binct
    if (binct lt 10 < ncat*0.02) then continue
    isrt = sort(lflx[ibin])
    arrlen = binct*frac < 500
    bins2 = [bins2,median(lexp[ibin[isrt[0:arrlen]]])]
    mins2 = [mins2,median(lflx[ibin[isrt[0:arrlen]]])]
    ;mins2 = [mins2,median(lflx[ibin[isrt[0:arrlen]]]-lerr[ibin[isrt[0:arrlen]]])]
endfor

;; construct photon- and background-limited functions
;re1 = poly_fit(bins1,mins1,1,yfit=yfit1)
;re2 = poly_fit(bins2,mins2,1,yfit=yfit2)
;photlim = re1[0]+re1[1]*lexp
;backlim = re2[0]+re2[1]*lexp
;xlim = alog10(10.^photlim+10.^backlim)

;p = plot(lexp,lflx,'.')
;p = plot(bins1,mins1,'or',sym_filled=1,/ov)
;p = plot(bins3,mins3,'og',sym_filled=1,/ov)
;p = plot(bins2,mins2,'oy',sym_filled=1,/ov)
;p = plot(lexp,photlim,'--r',thick=4,/ov)
;p = plot(lexp,backlim,'--b',thick=4,/ov)
;p = plot(lexp,xlim,'--g',thick=4,/ov)

bins = [bins1,bins3,bins2]
mins = [mins1,mins3,mins2]
re = poly_fit(bins,mins,2,yfit=yfit)
xlim = re[0]+re[1]*lexp+re[2]*lexp^2
;p = plot(lexp,xlim,'--',col='orange',thick=4,/ov)

return, 10.^xlim



END


;re = xray_flux_limit2(cat_exp_cha,cat_flx_cha,cat_err_cha,5)
;re = xray_flux_limit2(cat_exp_xmm,cat_flx_xmm,cat_err_xmm,5)
;re = xray_flux_limit2(cat_exp_nst,cat_flx_nst,cat_err_nst,5)







