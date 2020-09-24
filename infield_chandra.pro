PRO infield_chandra, INF_ONLY = inf_only


common _fits

;; Chandra Master
mast_path = file_search('/Users/ccarroll/Research/surveys/Chandra/*master*.fits')
;; 3XMM-DR8 Serendip Catalog Per-Observation Source Table
cat_path = '/Users/ccarroll/Research/surveys/Chandra/observation-source-2.fits'

;; OBSERVATIONS
;; number of elements in main sample
nsrc = n_elements(ra)

;; Chandra Master mast_chaive
mast_cha = mrdfits(mast_path[(sort(mast_path))[-1]],1)
;; Master Catalog is updated more frequently than CSC2! 
;; avoid spurious non-detections!
cat = mrdfits(cat_path,1)
;; use only OBSID that are in catalots
mast_id = mast_cha.obsid
cat_id = cat[where(cat.instrument eq 'ACIS',/null)].obsid
cat_id = cat_id[uniq(cat_id,sort(cat_id))]
match,mast_id,cat_id,imast,icat
iimast_cha = bytarr(n_elements(mast_cha))
iimast_cha[imast] = 1
mast_cha = mast_cha[where(iimast_cha,/null)]
;; use only mast_chaived sources
mast_cha = mast_cha[where(mast_cha.status eq 'ARCHIVED' or mast_cha.status eq 'OBSERVED',/null)]
mast_cha = mast_cha[where(mast_cha.detector eq 'ACIS-I',/null)]
;; output matched observation data
iiinf_cha = bytarr(nsrc)            ;; in field
texp_cha = dblarr(nsrc)-9999.       ;; exposure time (ontime)
sdst_cha = dblarr(nsrc)-9999.       ;; separation distance from field center

;; ACIS-I FOV is 16'x16'
;; https://heasarc.gsfc.nasa.gov/docs/chandra/chandra.html
fov_cha = 16.9*60./2.
;; match to master catalog
spherematch,ra,dec,mast_cha.ra,mast_cha.dec,fov_cha/3600.,isamp,ifield,sep_cntr,maxmatch=0
sep_cntr *= 3600.       ;; convert to arcsec
;; tag main sample sources as "in field"
iiinf_cha[isamp] = 1

if keyword_set(inf_only) then begin
    save,iiinf_cha,file='infield_cha.sav'
    return
endif

;; loop over observations and choose closest field
uind = isamp[uniq(isamp,sort(isamp))]
for i = 0,n_elements(uind)-1 do begin
    imatch = where(isamp eq uind[i],mlen)
    if (mlen eq 0) then stop
    min_sep = min(sep_cntr[imatch],imin)
    texp_cha[isamp[imatch[imin]]] = total(mast_cha[ifield[imatch]].exposure)
    sdst_cha[isamp[imatch[imin]]] = min_sep
endfor
save,mast_cha,iiinf_cha,texp_cha,sdst_cha,fov_cha,file='infield_cha.sav'


END




