PRO detect_nustar


common _fits
common _inf_nst
nsrc = n_elements(ra)

;; Combined NuSTAR Fields
pth1 = '/Users/ccarroll/Research/surveys/NuSTAR/nustar_combined_catalogs.fits'

;; DETECTIONS
;; load the combined NuSTAR sources
nst = mrdfits(pth1,1)

;; match to sample sources
;; NUSTAR PSF up to 70"
;; https://heasarc.gsfc.nasa.gov/docs/nustar/NuSTAR_observatory_guide-v1.0.pdf
;; we will use XMM 6.25" separation
psf_nst = 6.25
spherematch,ra,dec,nst.ra,nst.dec,psf_nst/3600.,isamp,inustar,sepx
sepx *= 3600.

tags = tag_names(nst)
nvars = n_elements(tags)

;; conflicting definition RA/Dec
ipos = where(strmatch(tags,'RA',/fold) or strmatch(tags,'DEC',/fold))
nst_vars = tags
nst_vars[ipos] = tags[ipos]+'_NST'

;; declare and initialize NuSTAR variables, matched to main sample
for i = 0,nvars-1 do begin
    re = execute('type = typename(nst.'+tags[i]+')')
    re = execute(nst_vars[i]+' = make_array(nsrc,/'+type+')')
    re = execute(nst_vars[i]+'[isamp] = nst[inustar].'+tags[i])
endfor
;; sample separation
sep_nst = dblarr(nsrc)
sep_nst[isamp] = sepx

;; boolean flag for cross-matched with NuSTAR combined catalog
iix_nst = bytarr(nsrc)
iix_nst[isamp] = 1
;; ensure valid photometry
phot = tags[where(strmatch(nst_vars,'?BF'),nphot)]
photerr = tags[where(strmatch(nst_vars,'E_?BF'),nphoterr)]
;; ensure valid exposure time
time = tags[where(strmatch(nst_vars,'?EXP'),ntime)]
;; boolean flag for valid detections in 3XMM
re = execute('iidet_nst = '+strjoin("(finite("+phot+") and "+phot+" gt 0. and finite("+photerr+") and "+photerr+" gt 0. and finite("+time+") and "+time+" gt 0.)"," or "))

;; save detection data
nst_str = 'NST,SEP_NST,IIX_NST,IIDET_NST,'+strjoin(nst_vars,',')
re = execute('save,'+nst_str+',/compress,file="detections_nst.sav"')

;; update in-field data
inew = where(iix_nst eq 1 and iiinf_nst eq 0,ctnew)
if (ctnew gt 0) then begin
    iiinf_nst[inew] = 1
    texp_nst[inew] = -1.
    sdst_nst[inew] = -1.
    inf_str = strjoin(scope_varname(common='_INF_NST'),',')
    re = execute('save,'+inf_str+',file="infield_nst.sav"')
endif


END





