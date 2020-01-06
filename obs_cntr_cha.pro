FUNCTION obs_cntr_cha, cat_ra, $
                       cat_dec, $
                       cat_exp, $
                       cat_fov, $
                       eff_fov


;; Chandra Master
mast_path = '/Users/ccarroll/Research/surveys/Chandra/*master*.fits'
;; 3XMM-DR8 Serendip Catalog Per-Observation Source Table
cat_path = '/Users/ccarroll/Research/surveys/Chandra/observation-source-2.fits'

;; OBSERVATIONS
;; number of elements in main sample
nsrc = n_elements(cat_ra)

;; Chandra Master Archive
arch = mrdfits(mast_path,1)
;; Master Catalog is updated more frequently than CSC2! 
;; avoid spurious non-detections!
cat = mrdfits(cat_path,1)
;; use only OBSID that are in catalots
mast_id = arch.obsid
cat_id = cat[where(cat.instrument eq 'ACIS',/null)].obsid
cat_id = cat_id[uniq(cat_id,sort(cat_id))]
match,mast_id,cat_id,imast,icat
iiarch = bytarr(n_elements(arch))
iiarch[imast] = 1
arch = arch[where(iiarch,/null)]
;; use only archived sources
arch = arch[where(arch.status eq 'ARCHIVED' or arch.status eq 'OBSERVED',/null)]
arch = arch[where(arch.detector eq 'ACIS-I',/null)]
;; ACIS-I FOV is 16'x16'
;; https://heasarc.gsfc.nasa.gov/docs/chandra/chandra.html
;; 50% Area
spherematch,cat_ra,cat_dec,arch.ra,arch.dec,cat_fov/3600.,icat,ifield,sep_cntr,maxmatch=0
sep_cntr *= 3600.       ;; convert to arcsec
;; output matched observation data
iimast = bytarr(nsrc)
;; exposure time and field center distance
mast_exp = dblarr(nsrc)
mast_sep = dblarr(nsrc)

;; loop over observations and choose closest field
uind = icat[uniq(icat,sort(icat))]
for i = 0,n_elements(uind)-1 do begin
    imatch = where(icat eq uind[i],mlen)
    if (mlen eq 0) then stop
    ;; closest exposure time CATALOG == MASTER
    ;min_diff = min(abs(cat_exp[icat[imatch]]-arch[ifield[imatch]].exposure),ipick)
    ;; largest exposure time MASTER
    max_diff = max(arch[ifield[imatch]].exposure,ipick)
    mast_exp[icat[imatch[ipick]]] = arch[ifield[imatch[ipick]]].exposure
    mast_sep[icat[imatch[ipick]]] = sep_cntr[imatch[ipick]]
endfor

;; field exposure exists and separation is inner 50%
iimast[where(mast_exp gt 0. and mast_sep le eff_fov,/null)] = 1

return, iimast


END




