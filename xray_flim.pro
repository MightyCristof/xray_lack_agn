FUNCTION xray_flim, expt, $
                    xflx, $
                    xerr;, $
                    ;PHOT_K = k1, $
                    ;BACK_K = k2, $
                    ;NRML = c2


;; convert to log-log
isort = sort(expt)
lexp = alog10(expt[isort])
lflx = alog10(xflx[isort])
lerr = sqrt((xerr[isort]/(xflx[isort]*alog(10)))^2)
;; number of sources in catalog
ncat = n_elements(lexp)
;binsz = floor(ncat*0.05)
mmlexp = minmax(lexp)
if (ncat lt 1000.) then nbins = 6. else $
                        nbins = 20.
binsz = diff(mmlexp)/nbins
;; PHOTON-LIMITED CASE
bins1 = []
mins1 = []
minerr1 = []

nobj = 0.
;; loop through bins
for i = 0,nbins-2 do begin
    binra = [mmlexp[0]+binsz*i,mmlexp[0]+binsz*(i+1)]
    ibin = where(lexp ge binra[0] and lexp le binra[1],binct)
    if (binct lt nobj*0.75) then break
    nobj = binct
    if (binct eq 0) then continue
    resistant_mean,lexp[ibin],2.,mn,sigmn,nrej,goodvec=ig
    binexp = lexp[ibin[ig]]
    binflx = lflx[ibin[ig]]
    binerr = lerr[ibin[ig]]
    ilo = [0:floor(n_elements(binflx)*0.1)<n_elements(binflx)-1]
    clipexp = binexp[ilo]
    clipflx = binflx[ilo]
    cliperr = binerr[ilo]
    loflx = min(clipflx-cliperr,imin)
    loexp = clipexp[imin]
    loerr = cliperr[imin]
    ;; ihi = [-(floor(n_elements(binflx)*0.1)<n_elements(binflx)):-1]
    clipexp = binexp[-(floor(n_elements(binflx)*0.1)<n_elements(binflx)):-1]
    clipflx = binflx[-(floor(n_elements(binflx)*0.1)<n_elements(binflx)):-1]
    cliperr = binerr[-(floor(n_elements(binflx)*0.1)<n_elements(binflx)):-1]
    hiflx = min(clipflx-cliperr,imax)
    hiexp = clipexp[imax]
    hierr = cliperr[imax]
    bins1 = [bins1,mn]
    if (hiflx-loflx eq 0.) then begin
        mins1 = [mins1,loflx]
        minerr1 = [minerr1,loerr]
    endif else begin
        mins1 = [mins1,interpol([loflx,hiflx],[loexp,hiexp],mn)]
        minerr1 = [minerr1,interpol([loerr,hierr],[loexp,hiexp],mn)]
    endelse
endfor

if (n_elements(bins1) lt 3.) then stop
dy1 = deriv(bins1,mins1)
resistant_mean,dy1,2.,mn1,sigm1,nrej1,goodvec=ig1
re1 = poly_fit(bins1[ig1],mins1[ig1],1)


if (i lt nbins-2) then begin
    bins2 = []
    mins2 = []
    minerr2 = []

    ;bintr = mmlexp[0]+binsz*i
    ;nbins = nbins-i
    i--
    for i = i,nbins-2 do begin
        binra = [mmlexp[0]+binsz*i,mmlexp[0]+binsz*(i+1)]
        ibin = where(lexp ge binra[0] and lexp le binra[1],binct)
        ;if (binct lt nobj*0.15) then break
        ;nobj = binct
        if (binct eq 0) then continue
        resistant_mean,lexp[ibin],2.,mn,sigmn,nrej,goodvec=ig
        binexp = lexp[ibin[ig]]
        binflx = lflx[ibin[ig]]
        binerr = lerr[ibin[ig]]
        ilo = [0:floor(n_elements(binflx)*0.1)<n_elements(binflx)-1]
        clipexp = binexp[ilo]
        clipflx = binflx[ilo]
        cliperr = binerr[ilo]
        loflx = min(clipflx-cliperr,imin)
        loexp = clipexp[imin]
        loerr = cliperr[imin]
        ;; ihi = [-(floor(n_elements(binflx)*0.1)<n_elements(binflx)):-1]
        clipexp = binexp[-(floor(n_elements(binflx)*0.1)<n_elements(binflx)):-1]
        clipflx = binflx[-(floor(n_elements(binflx)*0.1)<n_elements(binflx)):-1]
        cliperr = binerr[-(floor(n_elements(binflx)*0.1)<n_elements(binflx)):-1]
        hiflx = min(clipflx-cliperr,imax)
        hiexp = clipexp[imax]
        hierr = cliperr[imax]
        bins2 = [bins2,mn]
        mins2 = [mins2,interpol([loflx,hiflx],[loexp,hiexp],mn)]
        minerr2 = [minerr2,interpol([loerr,hierr],[loexp,hiexp],mn)]
    endfor
    dy2 = deriv(bins2,mins2)
    resistant_mean,dy2,1.0,mn2,sigm2,nrej2,goodvec=ig2
    re2 = poly_fit(bins2[ig2],mins2[ig2],1)
endif else begin
    bins2 = mmlexp[0]+binsz##(i+[1,2,3])
    mins2 = mins1[ig1[-1]]*[1,1,1]-sigm1
    minerr2 = minerr1[ig1[-1]]*[1,1,1]
    ig2 = [0,1,2]
    re2 = poly_fit(bins2,mins2,1)
endelse

if (ncat lt 1000.) then begin
    bins2 = mmlexp[0]+binsz##(i+[1,2,3])
    mins2 = mins1[ig1[-1]]*[1,1,1]-100.
    minerr2 = minerr1[ig1[-1]]*[1,1,1]
    ig2 = [0,1,2]
    re2 = poly_fit(bins2,mins2,1)
endif

;; flux limit
photlim = re1[0]+re1[1]*lexp
backlim = re2[0]+re2[1]*lexp
xlim = 10.^photlim+10.^backlim

return, xlim


if keyword_set(plt) then begin
    p = plot(lexp,lflx,'.',col='light grey')
    p = errorplot(bins1,mins1,minerr1,'o',linestyle='',/ov)
    p = errorplot(bins1[ig1],mins1[ig1],minerr1[ig1],'or',sym_filled=1,linestyle='',/ov)
    p = errorplot(bins2,mins2,minerr2,'o',linestyle='',/ov)
    p = errorplot(bins2[ig2],mins2[ig2],minerr2[ig2],'ob',sym_filled=1,linestyle='',/ov)
    p = plot(lexp,photlim,'-r',/ov)
    p = plot(lexp,backlim,'-b',/ov)
    p = plot(lexp,alog10(10.^photlim+10.^backlim),'-',col='purple',/ov)
endif


END
;       expt = cat_exp_cha_210
;       xflx = cat_flx_cha_210
;       xerr = cat_err_cha_210
;       
;       expt = cat_exp_xmm_210
;       xflx = cat_flx_xmm_210
;       xerr = cat_err_xmm_210
;       
;       expt = cat_exp_nst_210
;       xflx = cat_flx_nst_210
;       xerr = cat_err_nst_210
;       




;   p = plot(lexp,lflx,'.',col='light grey')
;   p = plot(bins1,mins1,'or',/ov)
;   p = plot(lexp,photlim,'--r',/ov)
;   p = plot(lexp,backlim,'--b',/ov)














;while nobj lt ncat-1 do begin
;    binra = [mmlexp[0]+binsz*counter,mmlexp[0]+binsz*(counter+1)]
;    counter++
;    ibin = where(lexp ge binra[0] and lexp le binra[1],binct)
;    if (binct eq 0) then continue
;    resistant_mean,lexp[ibin],2.,mn,sigmn,nrej,goodvec=ig
;    if (n_elements(ig) eq 0) then continue
;    ;; increase indices
;    nobj += binct
;    ;ibin = [counter:(counter+binsz-1)<(ncat-1):1]
;    ;if (n_elements(ig) lt 0.05*binsz) then stop;break
;    binexp = lexp[ibin[ig]]
;    binflx = lflx[ibin[ig]]
;    binerr = lerr[ibin[ig]]
;    ;mmexp = minmax(binexp)
;    ;; left-hand side
;    ilo = [0:floor(n_elements(binflx)*0.01)<n_elements(binflx)-1]
;    clipexp = binexp[ilo]
;    clipflx = binflx[ilo]
;    cliperr = binerr[ilo]
;    loflx = min(clipflx-cliperr,imin)
;    loexp = clipexp[imin]
;    loerr = cliperr[imin]
;    ;iup = [-(floor(n_elements(binflx)*0.01)<n_elements(binflx)):-1]
;    clipexp = binexp[-(floor(n_elements(binflx)*0.01)<n_elements(binflx)):-1]
;    clipflx = binflx[-(floor(n_elements(binflx)*0.01)<n_elements(binflx)):-1]
;    cliperr = binerr[-(floor(n_elements(binflx)*0.01)<n_elements(binflx)):-1]
;    hiflx = min(clipflx-cliperr,imax)
;    hiexp = clipexp[imax]
;    hierr = cliperr[imax]
;    ;; mid-point fit
;    bins1 = [bins1,mn1]
;    mins1 = [mins1,interpol([loflx,hiflx],[loexp,hiexp],mn1)]
;    minerr1 = [minerr1,interpol([loerr,hierr],[loexp,hiexp],mn1)]
;    ;; increase indices
;    ;counter = (counter+binsz)<(ncat-1)
;endwhile

;; BACKGROUND-LIMITED CASE
;bins2 = []
;mins2 = []
;minerr2 = []
;while counter lt ncat-1 do begin
;    ibin = [counter:(counter+binsz-1)<(ncat-1):1]
;    if (n_elements(ibin) lt 0.1*binsz) then break
;    binexp = lexp[ibin]
;    binflx = lflx[ibin]
;    binerr = lerr[ibin]
;    midpt = max(binexp)-width(minmax(binexp))/2.
;    ;sdev = stddev(abs(binexp-midpt))
;    ;itrim = where(binexp gt midpt-2.*sdev and binexp lt midpt+2.*sdev)
;    itrim = indgen(n_elements(binflx))
;    isort = sort(binflx[itrim])
;    midexp = binexp[itrim[isort]]
;    midflx = binflx[itrim[isort]]
;    miderr = binerr[itrim[isort]]
;    itrim = ceil(n_elements(midflx)*0.1)<n_elements(midflx)-1
;    finexp = midexp[0:itrim]
;    finflx = midflx[0:itrim]
;    finerr = miderr[0:itrim]
;    finval = finflx-finerr
;    bins2 = [bins2,midpt]
;    mins2 = [mins2,median(finval)]
;    minerr2 = [minerr2,sqrt(total(finerr^2))]
;    ;; increase indices
;    counter = (counter+binsz)<(ncat-1)
;endwhile
;dy2 = deriv(bins2,mins2)
;resistant_mean,dy2,1.0,mn2,sigm2,nrej2,goodvec=ig2
;re2 = poly_fit(bins2[ig2],mins2[ig2],1)



;bins1 = [bins1,mn1]
;mins1 = [mins1,interpol([loflx,hiflx],[loexp,hiexp],mn)]
;minerr1 = [minerr1,interpol([loflx,hiflx],[loexp,hiexp],mn)]]
;midpt = max(binexp)-width(minmax(binexp))/2.
;maxflx = min(binflx,imax)
;minexp = max(binexp,imin)
;minflx = binflx[imin]
;maxflx = binflx[imax]
;nbin = n_elements(ibin)
;
;;sdev = stddev(abs(binexp-midpt))
;;itrim = where(binexp gt midpt-2.*sdev and binexp lt midpt+2.*sdev)
;itrim = indgen(n_elements(binflx))
;isort = sort(binflx[itrim])
;midexp = binexp[itrim[isort]]
;midflx = binflx[itrim[isort]]
;miderr = binerr[itrim[isort]]
;itrim = ceil(n_elements(midflx)*.1)<n_elements(midflx)-1
;finexp = midexp[0:itrim]
;finflx = midflx[0:itrim]
;finerr = miderr[0:itrim]
;finval = finflx-finerr
;bins1 = [bins1,midpt]
;mins1 = [mins1,median(finval)]
;minerr1 = [minerr1,sqrt(total(finerr^2))]
;;; see if cumulative bin size is increasing
;if (counter eq 0) then begin
;    cumbinsz = bins1[0]
;    oldbinsz = cumbinsz
;endif else begin
;    cumbinsz = mean(diff(bins1))
;endelse
;;; if cumulative bin size is increasing, sources are becoming sparse; break loop        
;if (cumbinsz gt oldbinsz) then break else oldbinsz = cumbinsz 




;; initial fit
;; bin exposure times
;bins0 = [];[mm[0]:mm[1]:diff(mm)/15.]
;mins0 = [];dblarr(n_elements(bins))
;minerr0 = [];mins
;;; loop counter
;counter0 = 0
;;; loop through bins
;while counter0 lt ncat-1 do begin
;    ibin = [counter0:(counter0+binsz-1)<(ncat-1)]
;    binexp = lexp[ibin]
;    binflx = lflx[ibin]
;    binerr = lerr[ibin]
;    midpt = max(binexp)-width(minmax(binexp))/2.
;    ;sdev = stddev(abs(binexp-midpt))
;    ;itrim = where(binexp gt midpt-sdev and binexp lt midpt+sdev,trimct)
;    itrim = indgen(n_elements(binflx)) & trimct = n_elements(itrim)
;    if (trimct eq 0) then stop
;    midexp = binexp[itrim]
;    midflx = binflx[itrim]
;    miderr = binerr[itrim] 
;    ;; sort by flux
;    isort = sort(midflx)
;    midexp = midexp[isort]
;    midflx = midflx[isort]
;    miderr = miderr[isort]
;    ;; lower 10% of fluxes
;    ilower = [0:ceil(0.1*n_elements(midflx))]
;    bins0 = [bins0,midpt]
;    mins0 = [mins0,median(midflx[ilower])]
;    minerr0 = [minerr0,sqrt(total(miderr[ilower]^2))]
;;    if (n_elements(midflx) ge 50) then begin
;;        ;; clip outliers
;;        yhist = histogram(midflx,bin=scott(midflx),locations=xhist)
;;        ygauss = gaussfit(xhist,yhist,coeff,nterms=6)
;;        ikeep = where(midflx ge coeff[1]-coeff[2] and midflx lt coeff[1]+coeff[2],complement=iclip,ncomplement=nclip)
;;        midexp = midexp[ikeep]
;;        midflx = midflx[ikeep]
;;        miderr = miderr[ikeep]        
;;    endif
;;    itrim = ceil(n_elements(midflx)*0.02)
;;    finexp = midexp[0:itrim]
;;    finflx = midflx[0:itrim]
;;    finerr = miderr[0:itrim]
;;    bins0 = [bins0,midpt]
;;    mins0 = [mins0,finflx[0]]
;;    minerr0 = [minerr0,sqrt(total(finerr^2))]
;    ;; increase indices
;    counter0 = (counter0+binsz)<(ncat-1)
;endwhile
;dy0 = deriv(bins0,mins0)
;resistant_mean,dy0,2.0,mn0,sigm0,nrej0,goodvec=ig0
;re0 = poly_fit(bins0[ig0],mins0[ig0],2)
;;re0 = poly_fit(bins0,mins0,2)
;
;xlim = re0[0]+re0[1]*lexp+re0[2]*lexp^2
;if (ncat-counter0 lt binsz/2.) then return, 10.^xlim
;
;print, 'SHOULD HAVE RETURNED'








;dy1 = deriv(bins,mins)
;resistant_mean,dy1,1.5,mn1,sigm1,nrej1,goodvec=ig1
;k1 = mean(dy1[ig1[0:-2]])
;
;
;iturn = n_elements(bins)-ceil(n_elements(bins)/6.)
;dy1 = deriv(bins[0:iturn],mins[0:iturn])
;resistant_mean,dy1,1.,mn,sig_mn,nrej,goodvec=ig
;k1 = mean(dy1[1:ig[-1]/2.])
;
;
;
;for i = 0,n_elements(bins)-2 do begin
;    ;; find sources in bin
;    ibin = where(lexp ge bins[i] and lexp lt bins[i+1],binct)
;    if (binct eq 0.) then continue
;    ;print, binct
;    ;; bin source flux
;    binflx = lflx[ibin]
;    binerr = lerr[ibin]
;    ;; assume lower 2%
;    isort = sort(binflx)
;    iend = ceil(n_elements(binflx)*0.02)
;    binflx = (binflx[isort])[0:iend]
;    binerr = (binerr[isort])[0:iend]
;    ;; clip outliers
;    resistant_mean,binflx,2.,mn,sigmn,nrej,goodvec=ig
;    mins[i] = mean(binflx[ig])
;    mins_err[i] = mean(binerr[ig])
;endfor
;bins = bins[0:-2] & mins = mins[0:-2] & mins_err = mins_err[0:-2]
;;; add 1/2 step to bins
;bins -= width(bins)/2.
;;; remove any bins that contain no sources
;izero = where(mins eq 0.,complement=ival)
;bins = bins[ival]
;mins = mins[ival]
;if (n_elements(bins) lt 5) then goto, SINGLE_LIMIT
;;; take the derivative of this flux limit line
;iturn = n_elements(bins)-ceil(n_elements(bins)/6.)
;dy1 = deriv(bins[0:iturn],mins[0:iturn])
;
;;; iteratre through bins, starting from the least exposure time 
;;; and adding sources until one falls out
;;; assume this is a turning point in the slope and where 
;;; the background-limited case starts to impact
;resistant_mean,dy1,1.,mn,sig_mn,nrej,goodvec=ig
;k1 = mean(dy1[1:ig[-1]/2.])
;
;;for i = 0,n_elements(ig)-1 do begin
;;    resistant_mean,dy1[ig[0:i]],1.5,mn,sig_mn,nrej,goodvec=igphot
;;    if (n_elements(igphot) gt 2 and nrej gt 0) then break
;;endfor
;;; photon-limited flux-limit only slope
;;k1 = median(dy1[0:ig[-1]/2.]) * 1.05
;
;
;;; find the index where the background-limited case begins to dominate
;;iturn = ig[-1]+1
;;; take the right tail average and find the background-limited flux-limit slope
;if (n_elements(bins)-iturn ge 3) then dy2 = deriv(bins[iturn:-1],mins[iturn:-1]) else $
;if (n_elements(bins)-iturn eq 2) then dy2 = diff(mins[iturn:-1])/diff(bins[iturn:-1]) else $
;                                      dy2 = [-0.2,-0.2]
;resistant_mean,dy2,1.,mn,sig_mn,nrej,goodvec=igback
;k2 = median(dy2[igback]); > (-0.2)
;;; convert back to linear
;photlim = (expt/10.^3)^k1
;backlim = (expt/10.^3)^k2
;;; match photon- and background-limited case
;loc1 = value_locate(expt,10.^bins[iturn])
;c1 = photlim[loc1]/backlim[loc1]
;xlim = photlim + backlim*c1
;
;;; normalize full flux-limit to data
;loc2 = value_locate(expt,10.^bins[iturn])
;c2 = 10.^mins[iturn]/xlim[loc2]
;xlim *= c2
;
;;p = plot(expt,xflx,'.',col='dodger blue',/xlog,/ylog)
;;p = plot(expt,xlim,'--r',/ov)
;
;return, xlim
;
;
;;; data does not have photon- and background-limited case
;SINGLE_LIMIT: 
;bins = [mm[0]:mm[1]:diff(mm)/15.]
;mins = dblarr(n_elements(bins))
;
;;; loop through bins
;for i = 0,n_elements(bins)-2 do begin
;    ;; find sources in bin
;    ibin = where(lexp ge bins[i] and lexp lt bins[i+1],binct)
;    ;print, binct
;    if (binct eq 0) then continue
;    mins[i] = min(lflx[ibin])
;endfor
;;; add 1/2 step to bins
;bins += width(bins)/2.
;;; remove any bins that contain no sources
;izero = where(mins eq 0.,complement=ival)
;bins = bins[ival]
;mins = mins[ival]
;;; take the derivative of this flux limit line
;dy = deriv(bins,mins)
;
;
;;; iteratre through bins, starting from the least exposure time 
;;; and adding sources until one falls out
;;; assume this is a turning point in the slope and where 
;;; the background-limited case starts to impact
;resistant_mean,dy,1.,mn,sig_mn,nrej,goodvec=ig
;for i = 0,n_elements(ig)-1 do begin
;    resistant_mean,dy[ig[0:i]],1.,mn,sig_mn,nrej,goodvec=ig2
;    if (n_elements(ig2) gt 2 and nrej gt 0) then break
;endfor
;;; flux-limit slope
;k1 = median(dy[ig[ig2]])
;k2 = !Values.f_nan
;
;;; find the index to normalize flux-limit
;iturn = ig[ig2[-1]]+1
;;; convert back to linear
;xlim = (expt/10.^3)^k1
;;; normalize flux-limit to data
;loc = value_locate(expt,10.^bins[iturn])
;c2 = 10.^mins[iturn]/xlim[loc]
;xlim *= c2
;
;;p = plot(expt,xflx,'.',col='dodger blue',/xlog,/ylog)
;;p = plot(expt,xlim,'--r',/ov)
;
;return, xlim


;isort = sort(expt)
;lexp = alog10(expt[isort])
;lflx = alog10(xflx[isort])
;lerr = sqrt((xerr[isort]/(xflx[isort]*al
;;; number of sources in catalog
;ncat = n_elements(lexp)
;;binsz = floor(ncat*0.05)
;mmlexp = minmax(lexp)
;if (ncat lt 1000.) then nbins = 5. else 
;                        nbins = 20.
;binsz = diff(mmlexp)/nbins
;;; PHOTON-LIMITED CASE
;bins1 = []
;mins1 = []
;minerr1 = []
;;; loop counter
;;counter = 0
;nobj = 0.
;;; loop through bins
;for i = 0,nbins-2 do begin
;    binra = [mmlexp[0]+binsz*i,mmlexp[0]
;    ibin = where(lexp ge binra[0] and le
;    ;if (binct lt nobj/2.) then break
;    nobj = binct
;    if (binct eq 0) then continue
;    resistant_mean,lexp[ibin],2.,mn,sigm
;    binexp = lexp[ibin[ig]]
;    binflx = lflx[ibin[ig]]
;    binerr = lerr[ibin[ig]]
;    ilo = [0:floor(n_elements(binflx)*0.
;    clipexp = binexp[ilo]
;    clipflx = binflx[ilo]
;    cliperr = binerr[ilo]
;    loflx = min(clipflx-cliperr,imin)
;    loexp = clipexp[imin]
;    loerr = cliperr[imin]
;    ;; ihi = [-(floor(n_elements(binflx)
;    clipexp = binexp[-(floor(n_elements(
;    clipflx = binflx[-(floor(n_elements(
;    cliperr = binerr[-(floor(n_elements(
;    hiflx = min(clipflx-cliperr,imax)
;    hiexp = clipexp[imax]
;    hierr = cliperr[imax]
;    bins1 = [bins1,mn]
;    mins1 = [mins1,interpol([loflx,hiflx
;    minerr1 = [minerr1,interpol([loerr,h
;endfor
;if (n_elements(bins1) lt 3.) then stop
;dy1 = deriv(bins1,mins1)
;resistant_mean,dy1,2.0,mn1,sigm1,nrej1,g
;re1 = poly_fit(bins1[ig1],mins1[ig1],1)
;
;if (i lt nbins-2) then begin
;    bins2 = []
;    mins2 = []
;    minerr2 = []
;
;    for i = i,nbins-2 do begin
;        binra = [mmlexp[0]+binsz*i,mmlex
;        ibin = where(lexp ge binra[0] an
;        ;if (binct lt nobj/2.) then brea
;        ;nobj = binct
;        if (binct eq 0) then continue
;        resistant_mean,lexp[ibin],2.,mn,
;        binexp = lexp[ibin[ig]]
;        binflx = lflx[ibin[ig]]
;        binerr = lerr[ibin[ig]]
;        ilo = [0:floor(n_elements(binflx
;        clipexp = binexp[ilo]
;        clipflx = binflx[ilo]
;        cliperr = binerr[ilo]
;        loflx = min(clipflx-cliperr,imin
;        loexp = clipexp[imin]
;        loerr = cliperr[imin]
;        ;; ihi = [-(floor(n_elements(bin
;        clipexp = binexp[-(floor(n_eleme
;        clipflx = binflx[-(floor(n_eleme
;        cliperr = binerr[-(floor(n_eleme
;        hiflx = min(clipflx-cliperr,imax
;        hiexp = clipexp[imax]
;        hierr = cliperr[imax]
;        bins2 = [bins2,mn]
;        mins2 = [mins2,interpol([loflx,h
;        minerr2 = [minerr2,interpol([loe
;    endfor
;    dy2 = deriv(bins2,mins2)
;    resistant_mean,dy2,2.0,mn2,sigm2,nre
;    re2 = poly_fit(bins2[ig2],mins2[ig2]
;endif else if (i eq nbins) then begin
;    bins2 = mmlexp[0]+binsz##(i+[1,2,3])
;    mins2 = mins1[ig1[-1]]*[1,1,1]
;    minerr2 = minerr1[ig1[-1]]*[1,1,1]
;    re2 = poly_fit(bins2,mins2,1)
;endif 
;
;;; flux limit
;photlim = re1[0]+re1[1]*lexp
;backlim = re2[0]+re2[1]*lexp
;;!NULL = min(abs(photlim-backlim),icomb)
;;xlim = [photlim[0:icomb],backlim[icomb+
;xlim = 10.^photlim+10.^backlim
;;loc1 = value_locate(lexp,bins1[-1])
;;xlim = [photlim[0:loc1],backlim[loc1+1:
;;return, 10.^xlim
;
;p = plot(lexp,lflx,'.',col='light grey')
;o = errorplot(bins1,mins1,minerr1,'o',li
;o = errorplot(bins1[ig1],mins1[ig1],mine
;o = plot(lexp,photlim,'-r',/ov,yra=p.yra
;o = plot(bins2[ig2],mins2[ig2],'ob',/ov,
;o = plot(lexp,backlim,'-b',/ov,yra=p.yra
;o = plot(lexp,alog10(10.^photlim+10.^bac
;
;


