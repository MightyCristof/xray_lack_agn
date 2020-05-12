;; FOV 70% Effective Area
;; CHA - https://cxc.harvard.edu/proposer/POG/html/chap4.html#tth_sEc4.2.2
;; XMM - https://heasarc.gsfc.nasa.gov/docs/xmm/uhb/effareaoffaxis.html
;; NST - https://heasarc.gsfc.nasa.gov/docs/nustar/nustar_obsguide.pdf
PRO catalog_fxlim, MULTI_SN = multi_sn


common _inf_cha
common _inf_xmm
common _inf_nst
common _det_cha
common _det_xmm
common _det_nst
common _xconv


;; Energy band 2-10 keV
cat_exp = 'CAT_EXP'+xfield
cat_flx = 'CAT_FLX'+xfield
cat_err = 'CAT_ERR'+xfield
;; position
cat_ra = 'CAT_RA'+xfield
cat_dec = 'CAT_DEC'+xfield
;; FOV where effective area ³ 70%
fov_eff = 'FOV_EFF'+xfield
fov_eff_cha = 7.*60.
fov_eff_xmm = 7.*60.
fov_eff_nst = 5.*60.


for i = 0,nfield-1 do begin
    ;; choose instrument
    cat_fld = (strsplit(xfield[i],'_',/extract))[0]
    re = execute('ncat = n_elements('+cat_fld+')')
    ;; initialize instrument arrays
    re = execute(cat_exp[i]+' = dblarr(ncat)')
    re = execute(cat_flx[i]+' = dblarr(ncat)')
    re = execute(cat_err[i]+' = dblarr(ncat)')
    re = execute(cat_ra[i]+' = dblarr(ncat)')
    re = execute(cat_dec[i]+' = dblarr(ncat)')
    ;; energy conversion factors per instrument
    conv = xconv.(i)
    ;; source assigned flag
    iifill = bytarr(ncat)
    ;; loop through energy bands
    for b = 0,nxband[i]-1 do begin
        ;; valid detection exists for specified energy band
        re = execute('iivalid = '+cat_fld+'.'+tt.(i)[ixband.(i)[b]]+' gt 0. and '+cat_fld+'.'+ff.(i)[ixband.(i)[b]]+' gt 0. and '+cat_exp[i]+' eq 0.')
        ;; detection has not already been accounted for
        ifill = where(iivalid and iifill eq 0,nfill)
        if (nfill gt 0.) then begin
            ;; account for catalog source
            iifill[ifill] = 1
            re = execute(cat_exp[i]+'[ifill] = '+cat_fld+'[ifill].'+tt.(i)[ixband.(i)[b]])
            re = execute(cat_flx[i]+'[ifill] = '+cat_fld+'[ifill].'+ff.(i)[ixband.(i)[b]]+' * conv[ixband.(i)[b]]')
            re = execute(cat_err[i]+'[ifill] = '+cat_fld+'[ifill].'+ee.(i)[ixband.(i)[b]]+' * conv[ixband.(i)[b]]')
        endif
    endfor
    ;; grab positions for used sources
    ifill = where(iifill,nfill)
    if (nfill gt 0.) then begin
        re = execute(cat_ra[i]+'[ifill] = '+cat_fld+'[ifill].ra')
        re = execute(cat_dec[i]+'[ifill] = '+cat_fld+'[ifill].dec')
    endif
    ;; sort by exposure time
    re = execute('isort = sort('+cat_exp[i]+')')
    re = execute(cat_exp[i]+' = '+cat_exp[i]+'[isort]')
    re = execute(cat_flx[i]+' = '+cat_flx[i]+'[isort]')
    re = execute(cat_err[i]+' = '+cat_err[i]+'[isort]')
    re = execute(cat_ra[i]+' = '+cat_ra[i]+'[isort]')
    re = execute(cat_dec[i]+' = '+cat_dec[i]+'[isort]')
    ;; S/N cut
    re = execute('cat_sn = '+cat_flx[i]+'/'+cat_err[i])
    re = execute('icat_sn = where(cat_sn ge 3. and '+cat_exp[i]+' gt 0.,ncat_sn)')
    if (ncat_sn eq 0.) then stop
    re = execute(cat_exp[i]+' = '+cat_exp[i]+'[icat_sn]')
    re = execute(cat_flx[i]+' = '+cat_flx[i]+'[icat_sn]')
    re = execute(cat_err[i]+' = '+cat_err[i]+'[icat_sn]')
    re = execute(cat_ra[i]+' = '+cat_ra[i]+'[icat_sn]')
    re = execute(cat_dec[i]+' = '+cat_dec[i]+'[icat_sn]')
    ;; distance from field center (on-axis)
    re = execute('iimast'+xfield[i]+' = obs_cntr'+xfield[i]+'('+cat_ra[i]+','+cat_dec[i]+','+cat_exp[i]+','+fov_eff[i]+')')
    re = execute('imast = where(iimast'+xfield[i]+',nmast)')
    if (nmast eq 0.) then stop
    re = execute(cat_exp[i]+' = '+cat_exp[i]+'[imast]')
    re = execute(cat_flx[i]+' = '+cat_flx[i]+'[imast]')
    re = execute(cat_err[i]+' = '+cat_err[i]+'[imast]')
endfor

;; flux limit
cat_lim = 'CAT_LIM'+xfield
for i = 0,nfield-1 do re = execute(cat_lim[i]+' = xray_flux_limit('+cat_exp[i]+','+cat_flx[i]+','+cat_err[i]+',6)')

sav_vars = [cat_exp,cat_flx,cat_err,cat_lim,fov_eff]
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


END




