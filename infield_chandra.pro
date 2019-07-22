PRO infield_chandra


common _fits

;; Chandra Master
pth1 = '/Users/ccarroll/Research/surveys/Chandra/*master*.fits'

;; OBSERVATIONS
;; number of elements in main sample
nsrc = n_elements(ra)

;; use only archived sources
arch = mrdfits(pth1,1)
arch = arch[where(arch.status eq 'ARCHIVED' or arch.status eq 'OBSERVED')]
arch = arch[where(arch.detector eq 'ACIS-I')]

;; ACIS-I FOV is 16'x16'
;; https://heasarc.gsfc.nasa.gov/docs/chandra/chandra.html
fov_cha = 16.*60./2.

spherematch,ra,dec,arch.ra,arch.dec,fov_cha/3600.,isamp,ifield,sep_cntr,maxmatch=0
;; tag main sample sources as "in field"
iiinf_cha = bytarr(nsrc)
texp_cha = dblarr(nsrc)
sdst_cha = dblarr(nsrc)
iiinf_cha[isamp] = 1
texp_cha[isamp] = arch[ifield].exposure
sdst_cha[isamp] = sep_cntr*3600.            ;; convert to arcsec
save,iiinf_cha,texp_cha,sdst_cha,fov_cha,file='infield_cha.sav'


END




