FUNCTION obs_cntr_nst, cat_ra, $
                       cat_dec, $
                       cat_exp, $
                       cat_fov, $
                       eff_fov


;; Combined NuSTAR Fields
mast_path = '/Users/ccarroll/Research/surveys/NuSTAR/*master*.fits'
;; NuSTAR Catalogs
;cat_path = '/Users/ccarroll/Research/surveys/NuSTAR/combined_nustar_fields.fits'

;; OBSERVATIONS
nsrc = n_elements(cat_ra) ;; How many srcs
;;;; ========= X-ray Database Information ==============
;; Read in the NuSTAR observation information (HEASARC);
arch = mrdfits(mast_path,1)
;; Master Catalog is updated more frequently than individual catalogs! 
;; avoid spurious non-detections!
;cat = mrdfits(cat_path,1)
cat_id = create_nustar_master_obsid_list()
;; use only OBSID that are in catalots
mast_id = arch.obsid
match,mast_id,cat_id,imast,icat
iiarch = bytarr(n_elements(arch))
iiarch[imast] = 1
arch = arch[where(iiarch,/null)]
;; use only archived sources
arch = arch[where(arch.status eq 'ARCHIVED' or arch.status eq 'OBSERVED',/null)]
;; select SCIENCE mode
arch = arch[where(arch.observation_mode eq 'SCIENCE',/null)]
;; NuSTAR FOV is 13'x13'
;; https://heasarc.gsfc.nasa.gov/docs/nustar/nustar.html
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
    ;min_diff = min(abs(cat_exp[icat[imatch]]-arch[ifield[imatch]].ontime_a),ipick)
    ;; largest exposure time MASTER
    max_diff = max(arch[ifield[imatch]].ontime_a,ipick)
    mast_exp[icat[imatch[ipick]]] = arch[ifield[imatch[ipick]]].ontime_a
    mast_sep[icat[imatch[ipick]]] = sep_cntr[imatch[ipick]]
endfor

;; field exposure exists and separation is inner 50%
iimast[where(mast_exp gt 0. and mast_sep le eff_fov,/null)] = 1

return, iimast


END






