PRO xray_field_xmm


common _fits   

;; XMM Master
pth1 = '/Users/ccarroll/Research/surveys/XMM/*master*.fits'

;; OBSERVATIONS
;; number of elements in main sample
nsrc = n_elements(ra)

;; use only archived sources (possibly use )
arch = mrdfits(pth1,1)
arch = arch[where(arch.status eq 'ARCHIVED' or arch.status eq 'OBSERVED')]

;; XMM MOS FOV is ~33'x33'; use PN FOV inscribed circle--being conservative
fov = 26.4/2./60.
spherematch,ra,dec,arch.ra,arch.dec,fov,isamp,ifield,sep_cntr,maxmatch=0
;; tag main sample sources as "in field"

iiinf_xmm_field = bytarr(nsrc)            ;; in field
texp_xmm_field = dblarr(nsrc)         ;; exposure time (ontime)
sdst_xmm_field = dblarr(nsrc)         ;; separation distance from field center
iiinf_xmm_field[isamp] = 1
texp_xmm_field[isamp] = arch[ifield].duration
sdst_xmm_field[isamp] = sep_cntr
save,iiinf_xmm_field,texp_xmm_field,sdst_xmm_field,file='xfield_xmm.sav'


END




