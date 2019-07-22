PRO infield_xmm_newton


common _fits   

;; XMM Master
pth1 = '/Users/ccarroll/Research/surveys/XMM/*master*.fits'

;; OBSERVATIONS
;; number of elements in main sample
nsrc = n_elements(ra)

;; use only archived sources (possibly use )
arch = mrdfits(pth1,1)
arch = arch[where(arch.status eq 'ARCHIVED' or arch.status eq 'OBSERVED')]  ;; observed sources
arch = arch[where(arch.pn_time gt 0.)]                                      ;; ensure PN observation
arch = arch[where(arch.duration gt 0.)]                                     ;; sanity check
iimode = strmatch(arch.pn_mode,'*FLG*',/fold) or $                          ;; ensure Large-Window or Full-Frame mode
         strmatch(arch.pn_mode,'*FF*',/fold) or $
         strmatch(arch.pn_mode,'*EFF*',/fold)
arch = arch[where(iimode)]

;; XMM PN MOS FOV is ~27.5'x27.5'; use FOV inscribed circle--being conservative
;; https://heasarc.gsfc.nasa.gov/docs/xmm/xmm.html
fov_xmm = 27.5*60./2.

spherematch,ra,dec,arch.ra,arch.dec,fov_xmm/3600.,isamp,ifield,sep_cntr,maxmatch=0
;; tag main sample sources as "in field"

iiinf_xmm = bytarr(nsrc)            ;; in field
texp_xmm = dblarr(nsrc)         ;; exposure time (ontime)
sdst_xmm = dblarr(nsrc)         ;; separation distance from field center
iiinf_xmm[isamp] = 1
texp_xmm[isamp] = arch[ifield].duration
sdst_xmm[isamp] = sep_cntr*3600.            ;; convert to arcsec
save,iiinf_xmm,texp_xmm,sdst_xmm,fov_xmm,file='infield_xmm.sav'


END




