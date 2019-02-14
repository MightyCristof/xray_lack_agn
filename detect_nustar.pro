PRO detect_nustar


common _inf_fits
common _inf_nst

;; Combined NuSTAR Fields
pth1 = '/Users/ccarroll/Research/surveys/NuSTAR/nustar_fields.fits'

;; running for "Luminous AGN Lacking X-rays"

;; DETECTIONS
;; load the combined NuSTAR sources
nst = mrdfits(pth1,1)
nst = clean_detect_nustar(nst)

spherematch,ra_inf,dec_inf,nst.xra,nst.xdec,6./3600.,isamp,inustar,sep_nu
sep_nu *= 3600.

tags = tag_names(nst)
nvars = n_elements(tags)

;; conflicting definition RA/Dec
ipos = where(strmatch(tags,'XRA') or strmatch(tags,'XDEC'))
pos = (strsplit(tags[ipos],'X',/extract)).ToArray()
vars_nst = tags
vars_nst[ipos] = pos+'_NST'
vars_nst[where(strmatch(tags,'FIELD'))] += '_NST'

;; declare and initialize NuSTAR variables, matched to main sample
for i = 0,nvars-1 do begin
    re = execute('type = typename(nst.'+tags[i]+')')
    re = execute(vars_nst[i]+' = make_array(nsrc,/'+type+')')
    re = execute(vars_nst[i]+'[isamp] = nst[inustar].'+tags[i])
endfor

;; determine source detections (exposure time, flux, flux error)
xband = (strsplit(tags[where(strmatch(tags,'?BF'))],'BF',/extract,/regex)).ToArray()
iix = 'II'+xband
xband = tags[where(strmatch(tags,'?BF'))]
xberr = tags[where(strmatch(tags,'E_?BF'))]
xbexp = tags[where(strmatch(tags,'?EXP'))]
for i = 0,n_elements(xband)-1 do begin
    xstr = xband[i]+' gt 0. and '+xberr[i]+' gt 0. and '+xbexp[i]+' gt 0.'
    re = execute(iix[i]+' = '+xstr)
endfor
iixstr = strjoin(iix,' or ')
;; boolean flag for valid detection in any band
re = execute('iidet_nst = '+iixstr)
idet_nst = where(iidet_nst)

nst_str = 'NST,IIDET_NST,IDET_NST,'+strjoin(vars_nst,',')
re = execute('save,'+nst_str+',/compress,file="xdetect_nst.sav"')


END





