PRO detect_nustar


common _fits
common _inf_nst
nsrc = n_elements(ra)

;; Combined NuSTAR Fields
pth1 = '/Users/ccarroll/Research/surveys/NuSTAR/combined_nustar_fields.fits'

;; DETECTIONS
;; load the combined NuSTAR sources
nst = mrdfits(pth1,1)

spherematch,ra,dec,nst.ra,nst.dec,6./3600.,isamp,inustar,sep_nu
sep_nu *= 3600.

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

;; boolean flag for valid detection in any band
iidet_nst = bytarr(nsrc)
iidet_nst[isamp] = 1
idet_nst = where(iidet_nst)

;; save detection data
nst_str = 'NST,IIDET_NST,IDET_NST,'+strjoin(nst_vars,',')
re = execute('save,'+nst_str+',/compress,file="detections_nst.sav"')

;; update in-field data
inew = where(iidet_nst eq 1 and iiinf_nst eq 0,ctnew)
if (ctnew gt 0) then begin
    iiinf_nst[inew] = 1
    texp_nst[inew] = -9999.
    sdst_nst[inew] = -9999.
    inf_str = strjoin(scope_varname(common='_INF_NST'),',')
    re = execute('save,'+inf_str+',file="infield_nst.sav"')
endif


END





