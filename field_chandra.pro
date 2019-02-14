PRO field_chandra


common _fits

;; Chandra Master
pth1 = '/Users/ccarroll/Research/surveys/Chandra/*master*.fits'

;; OBSERVATIONS
;; number of elements in main sample
nsrc = n_elements(ra)

;; use only archived sources
arch = mrdfits(pth1,1)
arch = arch[where(arch.status eq 'ARCHIVED' or arch.status eq 'OBSERVED')]
arch = arch[where(arch.detector eq 'ACIS-I')]; or arch.detector eq 'HRC-I')]

;; ACIS-I FOV is 16'x16'; use inscribed circle--being conservative
fov = 16./2./60.
spherematch,ra,dec,arch.ra,arch.dec,fov,isamp,ifield,sep_cntr,maxmatch=0
;; tag main sample sources as "in field"
iiinf_cha_field = bytarr(nsrc)
texp_cha_field = dblarr(nsrc)
sdst_cha_field = dblarr(nsrc)
iiinf_cha_field[isamp] = 1
texp_cha_field[isamp] = arch[ifield].exposure
sdst_cha_field[isamp] = sep_cntr
save,iiinf_cha_field,texp_cha_field,sdst_cha_field,file='xfield_cha.sav'


END




