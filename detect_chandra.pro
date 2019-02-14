PRO detect_chandra
               

common _inf_fits
common _inf_cha

;; Chandra Source Catalog 2
pth1 = '/Users/ccarroll/Research/surveys/Chandra/chandra-source-catalog-2.fits'

;; running for "Luminous AGN Lacking X-rays"

;; DETECTIONS
;; load Chandra Source Catalog 2
cha = mrdfits(pth1,1)
;; use reliable source detections
cha = clean_detect_chandra(cha)

spherematch,ra_inf,dec_inf,cha.ra,cha.dec,6./3600.,isamp,imatch,sep_cha
tags = ['NAME','RA','DEC', $
        ;'CONF_FLAG','PILEUP_FLAG','SAT_SRC_FLAG', $
        'FLUX_POWLAW_APER90_B','FLUX_POWLAW_APER90_LOLIM_B','FLUX_POWLAW_APER90_HILIM_B', $
        'FLUX_POWLAW_APER90_H','FLUX_POWLAW_APER90_LOLIM_H','FLUX_POWLAW_APER90_HILIM_H', $
        'FLUX_POWLAW_APER90_M','FLUX_POWLAW_APER90_LOLIM_M','FLUX_POWLAW_APER90_HILIM_M', $
        'FLUX_POWLAW_APER90_S','FLUX_POWLAW_APER90_LOLIM_S','FLUX_POWLAW_APER90_HILIM_S', $
        'FLUX_POWLAW_APER90_U','FLUX_POWLAW_APER90_LOLIM_U','FLUX_POWLAW_APER90_HILIM_U', $
        'FLUX_POWLAW_APER90_W','FLUX_POWLAW_APER90_LOLIM_W','FLUX_POWLAW_APER90_HILIM_W', $
        'HARD_HS','HARD_HS_LOLIM','HARD_HS_HILIM', $
        'ACIS_TIME' $
        ]
nvars = n_elements(tags)

;; conflicting definition RA/Dec
ipos = where(strmatch(tags,'RA') or strmatch(tags,'DEC'))
cha_vars = tags
cha_vars[ipos] = tags[ipos]+'_CHA'

;; declare and initialize cha variables, matched to main sample
for i = 0,nvars-1 do begin
    re = execute('type = typename(cha.'+tags[i]+')')
    re = execute(cha_vars[i]+' = make_array(nsrc,/'+type+')')
    re = execute(cha_vars[i]+'[isamp] = cha[imatch].'+tags[i])
endfor

;; determine source detections (exposure time, flux, flux error)
xband = (strsplit(tags[where(strmatch(tags,'FLUX_POWLAW_APER90_?'))],'FLUX_POWLAW_APER90_?',/extract,/regex)).ToArray()
iix = 'II'+xband
xband = tags[where(strmatch(tags,'FLUX_POWLAW_APER90_?'))]
xhilim = tags[where(strmatch(tags,'FLUX_POWLAW_APER90_HILIM_?'))]
xlolim = tags[where(strmatch(tags,'FLUX_POWLAW_APER90_LOLIM_?'))]
for i = 0,n_elements(xband)-1 do begin
    xstr = xband[i]+' gt 0. and '+xlolim[i]+' gt 0. and '+xhilim[i]+' gt 0.'
    re = execute(iix[i]+' = '+xstr)
endfor
iixstr = '('+strjoin(iix," or ")+') and ACIS_TIME'
;; boolean flag for valid detection in any band
re = execute('iidet_cha = '+iixstr)
idet_cha = where(iidet_cha)

cha_str = 'CHA,IIDET_CHA,IDET_CHA,'+strjoin(cha_vars,',')
re = execute('save,'+cha_str+',/compress,file="xdetect_cha.sav"')


END




