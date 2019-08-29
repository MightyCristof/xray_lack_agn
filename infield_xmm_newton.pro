PRO infield_xmm_newton


common _fits   

;; XMM Master
mast_path = '/Users/ccarroll/Research/surveys/XMM/*master*.fits'
;; 3XMM-DR8 Serendip Catalog
cat_path = '/Users/ccarroll/Research/surveys/XMM/3XMM_DR8cat_v1.0.fits'

;; OBSERVATIONS
;; number of elements in main sample
nsrc = n_elements(ra)

;; XMM Master Archive
arch = mrdfits(mast_path,1)
;; Master Catalog is updated more frequently than 3XMM-DR8! 
;; avoid spurious non-detections!
cat = mrdfits(cat_path,1)
;; use only OBSID that are in catalots
mast_id = arch.obsid
cat_id = cat.obs_id
cat_id = cat_id[uniq(cat_id,sort(cat_id))]
match,mast_id,cat_id,imast,icat
iiarch = bytarr(n_elements(arch))
iiarch[imast] = 1
arch = arch[where(iiarch,/null)]
;; use only archived sources (possibly use )
arch = arch[where(arch.status eq 'ARCHIVED' or arch.status eq 'OBSERVED',/null)]  ;; observed sources
arch = arch[where(arch.pn_time gt 0.,/null)]                                      ;; ensure PN observation
arch = arch[where(arch.duration gt 0.,/null)]                                     ;; sanity check
iimode = strmatch(arch.pn_mode,'*FLG*',/fold) or $                          ;; ensure Large-Window or Full-Frame mode
         strmatch(arch.pn_mode,'*FF*',/fold) or $
         strmatch(arch.pn_mode,'*EFF*',/fold)
arch = arch[where(iimode,/null)]
;; XMM PN MOS FOV is ~27.5'x27.5'; use FOV inscribed circle--being conservative
;; https://heasarc.gsfc.nasa.gov/docs/xmm/xmm.html
fov_xmm = 27.5*60./2.
spherematch,ra,dec,arch.ra,arch.dec,fov_xmm/3600.,isamp,ifield,sep_cntr,maxmatch=0
sep_cntr *= 3600.       ;; convert to arcsec
;; output matched observation data
iiinf_xmm = bytarr(nsrc)        ;; in field
texp_xmm = dblarr(nsrc)         ;; exposure time (ontime)
sdst_xmm = dblarr(nsrc)         ;; separation distance from field center
;; tag main sample sources as "in field"
iiinf_xmm[isamp] = 1
;; loop over observations and choose closest field
uind = isamp[uniq(isamp,sort(isamp))]
for i = 0,n_elements(uind)-1 do begin
    imatch = where(isamp eq uind[i],mlen)
    if (mlen eq 0) then stop
    min_sep = min(sep_cntr[imatch],imin)
    texp_xmm[isamp[imatch[imin]]] = arch[ifield[imatch[imin]]].pn_time
    sdst_xmm[isamp[imatch[imin]]] = min_sep
endfor
save,iiinf_xmm,texp_xmm,sdst_xmm,fov_xmm,file='infield_xmm.sav'


END




