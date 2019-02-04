PRO xray_detect_wac


common _inf_fits

cat_dir = '/Users/ccarroll/Research/surveys/WISE/AGN Catalog 2017/table1-R90.fits'
;cat = '/Users/ccarroll/Research/surveys/WISE/AGN Catalog 2017/table2-C75.fits'
wac = mrdfits(cat_dir,1)

;; sources within 6 arcsec separation
sepdist = 6.
wac_sep = sepdist/3600.
spherematch,wac.ra,wac.dec,ra_inf,dec_inf,wac_sep,icat,isamp,sep_wac

;; byte flag for detected/matched in WISE AGN Catalog
iiwac = bytarr(nsrc)
iiwac[isamp] = 1
;; WISE AGN Catalog subset
wac = wac[icat]
struct_add_field,wac,'SEP_WAC',sep_wac

save,iiwac,wac,/compress,file='xdetect_wac.sav'


END




