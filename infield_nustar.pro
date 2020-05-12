PRO infield_nustar


common _fits


;; Combined NuSTAR Fields
mast_path = file_search('/Users/ccarroll/Research/surveys/NuSTAR/*master*.fits')
;; NuSTAR Catalogs
;cat_path = '/Users/ccarroll/Research/surveys/NuSTAR/combined_nustar_fields.fits'

;; OBSERVATIONS
nsrc = n_elements(ra) ;; How many srcs
;;;; ========= X-ray Database Information ==============
;; Read in the NuSTAR observation information (HEASARC);
mast_nst = mrdfits(mast_path[(sort(mast_path))[-1]],1)
;; Master Catalog is updated more frequently than individual catalogs! 
;; avoid spurious non-detections!
;cat = mrdfits(cat_path,1)
cat_id = create_nustar_master_obsid_list()
;; use only OBSID that are in catalots
mast_id = mast_nst.obsid
match,mast_id,cat_id,imast,icat
iimast_nst = bytarr(n_elements(mast_nst))
iimast_nst[imast] = 1
mast_nst = mast_nst[where(iimast_nst,/null)]
;; use only mast_nstived sources
mast_nst = mast_nst[where(mast_nst.status eq 'ARCHIVED' or mast_nst.status eq 'OBSERVED',/null)]
;; select SCIENCE mode
mast_nst = mast_nst[where(mast_nst.observation_mode eq 'SCIENCE',/null)]
;; output matched observation data
iiinf_nst = bytarr(nsrc)            ;; in field
texp_nst = dblarr(nsrc)-9999.       ;; exposure time (ontime)
sdst_nst = dblarr(nsrc)-9999.       ;; separation distance from field center

;; NuSTAR FOV is 13'x13'
;; https://heasarc.gsfc.nasa.gov/docs/nustar/nustar.html
fov_nst = 13.*60./2.
;; match to master catalog
spherematch,ra,dec,mast_nst.ra,mast_nst.dec,fov_nst/3600.,isamp,ifield,sep_cntr,maxmatch=0
sep_cntr *= 3600.       ;; convert to arcsec
;; tag main sample sources as "in field"
iiinf_nst[isamp] = 1

;; loop over observations and choose closest field
uind = isamp[uniq(isamp,sort(isamp))]
for i = 0,n_elements(uind)-1 do begin
    imatch = where(isamp eq uind[i],mlen)
    if (mlen eq 0) then stop
    min_sep = min(sep_cntr[imatch],imin)
    texp_nst[isamp[imatch[imin]]] = total(mast_nst[ifield[imatch]].ontime_a)
    sdst_nst[isamp[imatch[imin]]] = min_sep
endfor
save,mast_nst,iiinf_nst,texp_nst,sdst_nst,fov_nst,file='infield_nst.sav'


END





