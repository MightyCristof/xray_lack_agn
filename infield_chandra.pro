PRO infield_chandra


common _fits

;; Chandra Master
mast_path = '/Users/ccarroll/Research/surveys/Chandra/*master*.fits'
;; 3XMM-DR8 Serendip Catalog Per-Observation Source Table
cat_path = '/Users/ccarroll/Research/surveys/Chandra/observation-source-2.fits'

;; OBSERVATIONS
;; number of elements in main sample
nsrc = n_elements(ra)

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
fov_cha = 16.*60./2.
spherematch,ra,dec,arch.ra,arch.dec,fov_cha/3600.,isamp,ifield,sep_cntr,maxmatch=0
sep_cntr *= 3600.       ;; convert to arcsec
;; output matched observation data
iiinf_cha = bytarr(nsrc)
texp_cha = dblarr(nsrc)
sdst_cha = dblarr(nsrc)
;; tag main sample sources as "in field"
iiinf_cha[isamp] = 1
;; loop over observations and choose closest field
uind = isamp[uniq(isamp,sort(isamp))]
for i = 0,n_elements(uind)-1 do begin
    imatch = where(isamp eq uind[i],mlen)
    if (mlen eq 0) then stop
    min_sep = min(sep_cntr[imatch],imin)
    texp_cha[isamp[imatch[imin]]] = arch[ifield[imatch[imin]]].exposure
    sdst_cha[isamp[imatch[imin]]] = min_sep
endfor
save,iiinf_cha,texp_cha,sdst_cha,fov_cha,file='infield_cha.sav'


END




