PRO detect_chandra
               

common _fits
common _inf_cha
nsrc = n_elements(ra)

;; Chandra Source Catalog 2
pth1 = '/Users/ccarroll/Research/surveys/Chandra/chandra-source-catalog-2.fits'

;; DETECTIONS
;; load Chandra Source Catalog 2
cha = mrdfits(pth1,1)

;; add error column
cat_tags = tag_names(cha)
xband = ['B','H','M','S','U','W']
xbflx = cat_tags[where(strmatch(cat_tags,'FLUX_POWLAW_APER90_?'))]
xhilim = cat_tags[where(strmatch(cat_tags,'FLUX_POWLAW_APER90_HILIM_?'))]
xlolim = cat_tags[where(strmatch(cat_tags,'FLUX_POWLAW_APER90_LOLIM_?'))]
xberr = xbflx+'_ERR'
for i = 0,n_elements(xband)-1 do begin
    re = execute(xberr[i]+' = median([[cha.'+xbflx[i]+'-cha.'+xlolim[i]+'],[cha.'+xhilim[i]+'-cha.'+xbflx[i]+']],/even,dim=2)')
    re = execute('ifin = where(finite('+xberr[i]+'),complement=inan,ncomplement=nanct)')
    if (nanct gt 0.) then re = execute(xberr[i]+'[inan] = 0.')
    re = execute('struct_add_field,cha,xberr[i],'+xberr[i])
endfor

;; match to sample sources
spherematch,ra,dec,cha.ra,cha.dec,6./3600.,isamp,imatch,sep_cha

tags = ['NAME','RA','DEC','ACIS_TIME', $
        'FLUX_POWLAW_APER90_B','FLUX_POWLAW_APER90_B_ERR', $
        'FLUX_POWLAW_APER90_H','FLUX_POWLAW_APER90_H_ERR', $
        'FLUX_POWLAW_APER90_M','FLUX_POWLAW_APER90_M_ERR', $
        'FLUX_POWLAW_APER90_S','FLUX_POWLAW_APER90_S_ERR', $
        'FLUX_POWLAW_APER90_U','FLUX_POWLAW_APER90_U_ERR', $
        'FLUX_POWLAW_APER90_W','FLUX_POWLAW_APER90_W_ERR', $
        'HARD_HS','HARD_HS_LOLIM','HARD_HS_HILIM', $
        'DITHER_WARNING_FLAG','PILEUP_FLAG','SAT_SRC_FLAG','VAR_FLAG', $
        'STREAK_SRC_FLAG','VAR_INTER_HARD_FLAG','MAN_ADD_FLAG' $
        ]
nvars = n_elements(tags)

;; conflicting definition RA/Dec
ipos = where(strmatch(tags,'RA') or strmatch(tags,'DEC'))
cha_vars = tags
cha_vars[ipos] = tags[ipos]+'_CHA'

;; declare and initialize cha variables, matched to main sample
for i = 0,nvars-1 do begin
    re = execute('type = typename(cha.'+tags[i]+')')
    if (type eq 'LONG64') then type = 'L64'
    re = execute(cha_vars[i]+' = make_array(nsrc,/'+type+')')
    re = execute(cha_vars[i]+'[isamp] = cha[imatch].'+tags[i])
endfor

;; convert 90% aperture flux to 100%
aper90 = tags[where(strmatch(cha_vars,'*APER90*'),naper)]
if (naper gt 0.) then for i = 0,naper-1 do re = execute(aper90[i]+' *= 1.1')


;; boolean flag for detection in CSC2
iidet_cha = bytarr(nsrc)
iidet_cha[isamp] = 1
idet_cha = where(iidet_cha)

;; save data
cha_str = 'CHA,IIDET_CHA,IDET_CHA,'+strjoin(cha_vars,',')
re = execute('save,'+cha_str+',/compress,file="detections_cha.sav"')


END




