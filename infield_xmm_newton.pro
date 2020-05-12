PRO infield_xmm_newton


common _fits   

;; XMM Master
mast_path = file_search('/Users/ccarroll/Research/surveys/XMM/*master*.fits')
;; 3XMM-DR8 Serendip Catalog
cat_path = '/Users/ccarroll/Research/surveys/XMM/3XMM_DR8cat_v1.0.fits'

;; OBSERVATIONS
;; number of elements in main sample
nsrc = n_elements(ra)

;; XMM Master mast_xmmive
mast_xmm = mrdfits(mast_path[(sort(mast_path))[-1]],1)
;; Master Catalog is updated more frequently than 3XMM-DR8! 
;; avoid spurious non-detections!
cat = mrdfits(cat_path,1)
;; use only OBSID that are in catalots
mast_id = mast_xmm.obsid
cat_id = cat.obs_id
cat_id = cat_id[uniq(cat_id,sort(cat_id))]
match,mast_id,cat_id,imast,icat
iimast_xmm = bytarr(n_elements(mast_xmm))
iimast_xmm[imast] = 1
mast_xmm = mast_xmm[where(iimast_xmm,/null)]
;; use only mast_xmmived sources (possibly use )
mast_xmm = mast_xmm[where(mast_xmm.status eq 'ARCHIVED' or mast_xmm.status eq 'OBSERVED',/null)]  ;; observed sources
mast_xmm = mast_xmm[where(mast_xmm.pn_time gt 0.,/null)]                                      ;; ensure PN observation
mast_xmm = mast_xmm[where(mast_xmm.duration gt 0.,/null)]                                     ;; sanity check
iimode = strmatch(mast_xmm.pn_mode,'*FLG*',/fold) or $                          ;; ensure Large-Window or Full-Frame mode
         strmatch(mast_xmm.pn_mode,'*FF*',/fold) or $
         strmatch(mast_xmm.pn_mode,'*EFF*',/fold)
mast_xmm = mast_xmm[where(iimode,/null)]
;; output matched observation data
iiinf_xmm = bytarr(nsrc)            ;; in field
texp_xmm = dblarr(nsrc)-9999.       ;; exposure time (ontime)
sdst_xmm = dblarr(nsrc)-9999.       ;; separation distance from field center

;; XMM PN MOS FOV is ~27.5'x27.5'; use FOV inscribed circle--being conservative
;; https://heasarc.gsfc.nasa.gov/docs/xmm/xmm.html
fov_xmm = 33.*60./2.
;; match to master catalog
spherematch,ra,dec,mast_xmm.ra,mast_xmm.dec,fov_xmm/3600.,isamp,ifield,sep_cntr,maxmatch=0
sep_cntr *= 3600.       ;; convert to arcsec
;; tag main sample sources as "in field"
iiinf_xmm[isamp] = 1

;; loop over observations and choose closest field
uind = isamp[uniq(isamp,sort(isamp))]
for i = 0,n_elements(uind)-1 do begin
    imatch = where(isamp eq uind[i],mlen)
    if (mlen eq 0) then stop
    min_sep = min(sep_cntr[imatch],imin)
    texp_xmm[isamp[imatch[imin]]] = total(mast_xmm[ifield[imatch]].pn_time)
    sdst_xmm[isamp[imatch[imin]]] = min_sep
endfor
save,mast_xmm,iiinf_xmm,texp_xmm,sdst_xmm,fov_xmm,file='infield_xmm.sav'


END




