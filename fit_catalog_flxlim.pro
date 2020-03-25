;; FOV 70% Effective Area
;; CHA - https://cxc.harvard.edu/proposer/POG/html/chap4.html#tth_sEc4.2.2
;; XMM - https://heasarc.gsfc.nasa.gov/docs/xmm/uhb/effareaoffaxis.html
;; NST - https://heasarc.gsfc.nasa.gov/docs/nustar/nustar_obsguide.pdf
PRO fit_catalog_flxlim, MULTI_SN = multi_sn, $
                        PLT = plt


common _inf_cha
common _inf_xmm
common _inf_nst
common _det_cha
common _det_xmm
common _det_nst
common _softx


;; Energy band 2-10 keV
cat_exp = 'CAT_EXP'+xfield
cat_flx = 'CAT_FLX'+xfield
cat_err = 'CAT_ERR'+xfield
;; position
cat_ra = 'CAT_RA'+xfield
cat_dec = 'CAT_DEC'+xfield
;; used
iimast = 'IIMAST'+xfield
cat_fov = 'FOV'+xfield
;; FOV 70% effective area
eff_fov = 'EFF'+xfield

for f = 0,nfield-1 do begin
    ;; choose instrument
    cat_fld = (strsplit(xfield[f],'_',/extract))[0]
    re = execute('ncat = n_elements('+cat_fld+')')
    ;; initialize instrument arrays
    re = execute(cat_exp[f]+' = dblarr(ncat)')
    re = execute(cat_flx[f]+' = dblarr(ncat)')
    re = execute(cat_err[f]+' = dblarr(ncat)')
    re = execute(cat_ra[f]+' = dblarr(ncat)')
    re = execute(cat_dec[f]+' = dblarr(ncat)')
    ;; energy conversion factors per instrument
    conv = xconv.(f)
    ;; sort energy bands by minimum conversion factor
    iorder = sort(abs((conv-conv[ixband[f]])/conv[ixband[f]]))
    ;; use only one filter for each source; boolean flag
    iiused = bytarr(ncat)
    ;; loop through energy bands
    for b = 0,nxband[f]-1 do begin
        ;; detection exists
        re = execute('iifill = '+cat_fld+'.'+tt.(f)[iorder[b]]+' gt 0. and '+cat_fld+'.'+ff.(f)[iorder[b]]+' gt 0. and '+cat_exp[f]+' eq 0.')
        ;; source is not already accounted for
        ifill = where(iifill eq 1 and iiused eq 0,nfill)
        iiused[ifill] = 1
        if (nfill gt 0.) then begin
            re = execute(cat_exp[f]+'[ifill] = '+cat_fld+'[ifill].'+tt.(f)[iorder[b]])
            re = execute(cat_flx[f]+'[ifill] = '+cat_fld+'[ifill].'+ff.(f)[iorder[b]]+' * conv[iorder[b]]')
            re = execute(cat_err[f]+'[ifill] = '+cat_fld+'[ifill].'+ee.(f)[iorder[b]]+' * conv[iorder[b]]')
        endif
    endfor
    ;; grab positions for used sources
    iused = where(iiused,nused)
    if (nused gt 0.) then begin
        re = execute(cat_ra[f]+'[iused] = '+cat_fld+'[iused].ra')
        re = execute(cat_dec[f]+'[iused] = '+cat_fld+'[iused].dec')
    endif
    ;; sort by exposure time
    re = execute('isort = sort('+cat_exp[f]+')')
    re = execute(cat_exp[f]+' = '+cat_exp[f]+'[isort]')
    re = execute(cat_flx[f]+' = '+cat_flx[f]+'[isort]')
    re = execute(cat_err[f]+' = '+cat_err[f]+'[isort]')
    re = execute(cat_ra[f]+' = '+cat_ra[f]+'[isort]')
    re = execute(cat_dec[f]+' = '+cat_dec[f]+'[isort]')
    ;; S/N cut
    re = execute('cat_sn = '+cat_flx[f]+'/'+cat_err[f])
    re = execute('icat_sn = where(cat_sn ge 3. and '+cat_exp[f]+' gt 0.,ng)')
    if (ng eq 0.) then stop
    re = execute(cat_exp[f]+' = '+cat_exp[f]+'[icat_sn]')
    re = execute(cat_flx[f]+' = '+cat_flx[f]+'[icat_sn]')
    re = execute(cat_err[f]+' = '+cat_err[f]+'[icat_sn]')
    re = execute(cat_ra[f]+' = '+cat_ra[f]+'[icat_sn]')
    re = execute(cat_dec[f]+' = '+cat_dec[f]+'[icat_sn]')
    ;; distance from field center (on-axis)
    re = execute('iimast'+xfield[f]+' = obs_cntr'+xfield[f]+'('+cat_ra[f]+','+cat_dec[f]+','+cat_exp[f]+','+cat_fov[f]+','+eff_fov[f]+')')
    re = execute('imast = where(iimast'+xfield[f]+',nmast)')
    if (nmast eq 0.) then stop
    re = execute(cat_exp[f]+' = '+cat_exp[f]+'[imast]')
    re = execute(cat_flx[f]+' = '+cat_flx[f]+'[imast]')
    re = execute(cat_err[f]+' = '+cat_err[f]+'[imast]')
endfor

;; flux limit
cat_lim = 'CAT_LIM'+xfield
for f = 0,nfield-1 do re = execute(cat_lim[f]+' = xray_flux_limit('+cat_exp[f]+','+cat_flx[f]+','+cat_err[f]+',6)')

sav_vars = [cat_exp,cat_flx,cat_err,cat_lim]
sav_inds = []

if keyword_set(multi_sn) then begin
    ;;; S/N ³ 3
    ;cat_exp3 = 'CAT_EXP3'+xfield
    ;cat_flx3 = 'CAT_FLX3'+xfield
    ;cat_err3 = 'CAT_ERR3'+xfield
    ;;; S/N cut
    ;for f = 0,nfield-1 do begin
    ;    re = execute('cat_sn = '+cat_flx[f]+'/'+cat_err[f])
    ;    re = execute('icat_sn = where(cat_sn ge 3.,ng)')
    ;    if (ng eq 0.) then stop
    ;    re = execute(cat_exp3[f]+' = '+cat_exp[f]+'[icat_sn]')
    ;    re = execute(cat_flx3[f]+' = '+cat_flx[f]+'[icat_sn]')
    ;    re = execute(cat_err3[f]+' = '+cat_err[f]+'[icat_sn]')
    ;endfor
    ;;; flux limit
    ;cat_lim3 = 'CAT_LIM3'+xfield
    ;for f = 0,nfield-1 do re = execute(cat_lim3[f]+' = xray_flux_limit('+cat_exp3[f]+','+cat_flx3[f]+','+cat_err3[f]+',6)')
    ;;; save vars
    ;sav_vars = [sav_vars,cat_exp3,cat_flx3,cat_err3,cat_lim3]
    ;sav_inds = [sav_inds]

    ;; S/N ³ 5
    cat_exp5 = 'CAT_EXP5'+xfield
    cat_flx5 = 'CAT_FLX5'+xfield
    cat_err5 = 'CAT_ERR5'+xfield
    ;; S/N cut
    for f = 0,nfield-1 do begin
        re = execute('cat_sn = '+cat_flx[f]+'/'+cat_err[f])
        re = execute('icat_sn = where(cat_sn ge 5.,ng)')
        if (ng eq 0.) then stop
        re = execute(cat_exp5[f]+' = '+cat_exp[f]+'[icat_sn]')
        re = execute(cat_flx5[f]+' = '+cat_flx[f]+'[icat_sn]')
        re = execute(cat_err5[f]+' = '+cat_err[f]+'[icat_sn]')
    endfor
    ;; flux limit
    cat_lim5 = 'CAT_LIM5'+xfield
    for f = 0,nfield-1 do re = execute(cat_lim5[f]+' = xray_flux_limit('+cat_exp5[f]+','+cat_flx5[f]+','+cat_err5[f]+',6)')
    ;; save vars
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




