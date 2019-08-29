PRO infield_nustar


common _fits


;; Combined NuSTAR Fields
mast_path = '/Users/ccarroll/Research/surveys/NuSTAR/*master*.fits'
;; NuSTAR Catalogs
;cat_path = '/Users/ccarroll/Research/surveys/NuSTAR/combined_nustar_fields.fits'

;; OBSERVATIONS
nsrc = n_elements(ra) ;; How many srcs
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
fov_nst = 13.*60./2.
spherematch,ra,dec,arch.ra,arch.dec,fov_nst/3600.,isamp,ifield,sep_cntr,maxmatch=0
sep_cntr *= 3600.       ;; convert to arcsec
;; output matched observation data
iiinf_nst = bytarr(nsrc)
texp_nst = dblarr(nsrc)
sdst_nst = dblarr(nsrc)
;; tag main sample sources as "in field"
iiinf_nst[isamp] = 1
;; loop over observations and choose closest field
uind = isamp[uniq(isamp,sort(isamp))]
for i = 0,n_elements(uind)-1 do begin
    imatch = where(isamp eq uind[i],mlen)
    if (mlen eq 0) then stop
    min_sep = min(sep_cntr[imatch],imin)
    texp_nst[isamp[imatch[imin]]] = arch[ifield[imatch[imin]]].ontime_a
    sdst_nst[isamp[imatch[imin]]] = min_sep
endfor
save,iiinf_nst,texp_nst,sdst_nst,fov_nst,file='infield_nst.sav'


END





