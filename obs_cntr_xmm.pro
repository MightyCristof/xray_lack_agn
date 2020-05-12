FUNCTION obs_cntr_xmm, cat_ra, $
                       cat_dec, $
                       cat_exp, $
                       cat_fov


common _inf_xmm

;;; number of catalog sources
nxsrc = n_elements(cat_ra)
;; match catalog to master observations
spherematch,cat_ra,cat_dec,mast_xmm.ra,mast_xmm.dec,cat_fov/3600.,icat,ifield,sep_cntr,maxmatch=0
sep_cntr *= 3600.       ;; convert to arcsec
;; output matched observation data
iimast = bytarr(nxsrc)
;; exposure time and field center distance
mast_exp = dblarr(nxsrc)
mast_sep = dblarr(nxsrc)

;; loop over observations and choose closest field
uind = icat[uniq(icat,sort(icat))]
for i = 0,n_elements(uind)-1 do begin
    imatch = where(icat eq uind[i],mlen)
    if (mlen eq 0) then stop
    ;; closest exposure time CATALOG == MASTER
    ;min_diff = min(abs(cat_exp[icat[imatch]]-arch[ifield[imatch]].pn_time),ipick)
    ;; largest exposure time MASTER
    max_diff = max(mast_xmm[ifield[imatch]].pn_time,ipick)
    mast_exp[icat[imatch[ipick]]] = mast_xmm[ifield[imatch[ipick]]].pn_time
    mast_sep[icat[imatch[ipick]]] = sep_cntr[imatch[ipick]]
endfor

;; field exposure exists
iimast[where(mast_exp gt 0.,/null)] = 1

return, iimast


END





