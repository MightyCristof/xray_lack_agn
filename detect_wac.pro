PRO detect_wac


common _fits
nsrc = n_elements(ra)

cat_dir = '/Users/ccarroll/Research/surveys/WISE/AGN Catalog 2017/table1-R90.fits'
;cat = '/Users/ccarroll/Research/surveys/WISE/AGN Catalog 2017/table2-C75.fits'
wac = mrdfits(cat_dir,1)

;; sources within 6 arcsec separation
sepdist = 6.
wac_sep = sepdist/3600.
spherematch,wac.ra,wac.dec,ra,dec,wac_sep,icat,isamp,sep_wac

;; byte flag for detected/matched in WISE AGN Catalog
iidet_wac = bytarr(nsrc)
iidet_wac[isamp] = 1

;; save data
save,iidet_wac,wac,/compress,file='detections_wac.sav'


END




