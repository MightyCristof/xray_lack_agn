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
    if (binct lt floor(nobj*0.75)) then break
    nobj = binct
    if (binct eq 0) then continue
    resistant_mean,lexp[ibin],2.,midpt,sigmn,nrej,goodvec=ig
    binexp = lexp[ibin[ig]]
    binflx = lflx[ibin[ig]]
    binerr = lerr[ibin[ig]]
    ilo = [0:floor(n_elements(binflx)*0.1)<n_elements(binflx)-1]
    clipexp = binexp[ilo]
    clipflx = binflx[ilo]
    cliperr = binerr[ilo]
    resistant_mean,clipflx-cliperr,5.,mn,sigmn,nrej,goodvec=ig
    loflx = min(clipflx[ig]-cliperr[ig],imin)
    loexp = clipexp[ig[imin]]
    loerr = cliperr[ig[imin]]
    ;; ihi = [-(floor(n_elements(binflx)*0.1)<n_elements(binflx)):-1]
    clipexp = binexp[-(floor(n_elements(binflx)*0.1)<n_elements(binflx)):-1]
    clipflx = binflx[-(floor(n_elements(binflx)*0.1)<n_elements(binflx)):-1]
    cliperr = binerr[-(floor(n_elements(binflx)*0.1)<n_elements(binflx)):-1]
    resistant_mean,clipflx-cliperr,5.,mn,sigmn,nrej,goodvec=ig
    hiflx = min(clipflx[ig]-cliperr[ig],imax)
    hiexp = clipexp[ig[imax]]
    hierr = cliperr[ig[imax]]
    bins1 = [bins1,midpt]
    if (hiflx-loflx eq 0.) then begin
        mins1 = [mins1,loflx]
        minerr1 = [minerr1,loerr]
    endif else begin
        mins1 = [mins1,interpol([loflx,hiflx],[loexp,hiexp],midpt)]
        minerr1 = [minerr1,interpol([loerr,hierr],[loexp,hiexp],midpt)]
    endelse
endfor

if (n_elements(bins1) lt 3.) then stop
dy1 = deriv(bins1,mins1)
resistant_mean,dy1,1.,mn1,sigm1,nrej1,goodvec=ig1
re1 = poly_fit(bins1[ig1],mins1[ig1],1)

if ((nbins gt 6) and (i lt nbins-2)) then begin
    bins2 = []
    mins2 = []
    minerr2 = []
    ;bintr = mmlexp[0]+binsz*i
    ;nbins = nbins-i
    if (i gt nbins/2) then i -= 1
    for i = i,nbins-2 do begin
        binra = [mmlexp[0]+binsz*i,mmlexp[0]+binsz*(i+1)]
        ibin = where(lexp ge binra[0] and lexp le binra[1],binct)
        if (binct lt floor(nobj*0.5)) then break
        nobj = binct
        if (binct eq 0) then continue
        resistant_mean,lexp[ibin],2.,midpt,sigmn,nrej,goodvec=ig
        binexp = lexp[ibin[ig]]
        binflx = lflx[ibin[ig]]
        binerr = lerr[ibin[ig]]
        ilo = [0:floor(n_elements(binflx)*0.1)<n_elements(binflx)-1]
        clipexp = binexp[ilo]
        clipflx = binflx[ilo]
        cliperr = binerr[ilo]
        resistant_mean,clipflx-cliperr,5.,mn,sigmn,nrej,goodvec=ig
        loflx = min(clipflx[ig]-cliperr[ig],imin)
        loexp = clipexp[ig[imin]]
        loerr = cliperr[ig[imin]]
        ;; ihi = [-(floor(n_elements(binflx)*0.1)<n_elements(binflx)):-1]
        clipexp = binexp[-(floor(n_elements(binflx)*0.1)<n_elements(binflx)):-1]
        clipflx = binflx[-(floor(n_elements(binflx)*0.1)<n_elements(binflx)):-1]
        cliperr = binerr[-(floor(n_elements(binflx)*0.1)<n_elements(binflx)):-1]
        resistant_mean,clipflx-cliperr,5.,mn,sigmn,nrej,goodvec=ig
        hiflx = min(clipflx[ig]-cliperr[ig],imax)
        hiexp = clipexp[ig[imax]]
        hierr = cliperr[ig[imax]]
        bins2 = [bins2,midpt]
        mins2 = [mins2,interpol([loflx,hiflx],[loexp,hiexp],midpt)]
        minerr2 = [minerr2,interpol([loerr,hierr],[loexp,hiexp],midpt)]
    endfor
    dy2 = deriv(bins2,mins2)
    resistant_mean,dy2,1.,mn2,sigm2,nrej2,goodvec=ig2
    re2 = poly_fit(bins2[ig2],mins2[ig2],1)
endif

;; calculate flux limit in log space
if (n_elements(re2) gt 0) then begin
    photlim = re1[0]+re1[1]*lexp
    backlim = re2[0]+re2[1]*lexp
    xlim = 10.^photlim+10.^backlim
endif else begin
    photlim = re1[0]+re1[1]*lexp
    xlim = 10.^photlim
endelse

return, xlim



END







;p = plot(lexp,lflx,'.',col='light grey',yra=[-30,0])
;p = errorplot(bins1,mins1,minerr1,'o',linestyle='',/ov)
;p = errorplot(bins1[ig1],mins1[ig1],minerr1[ig1],'or',sym_filled=1,linestyle='',/ov)
;p = errorplot(bins2,mins2,minerr2,'o',linestyle='',/ov)
;p = errorplot(bins2[ig2],mins2[ig2],minerr2[ig2],'ob',sym_filled=1,linestyle='',/ov)
;p = plot(lexp,photlim,'-r',/ov)
;p = plot(lexp,backlim,'-b',/ov)
;p = plot(lexp,alog10(10.^photlim+10.^backlim),'-',col='purple',/ov)
