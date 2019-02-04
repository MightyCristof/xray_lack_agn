PRO xray_detect_xmm


common _inf_fits
common _inf_xmm

;; XMM DR8 Source Catalog
pth1 = '/Users/ccarroll/Research/surveys/XMM/3XMM_DR8cat_v1.0.fits'

;; running for "Luminous AGN Lacking X-rays"

;; DETECTIONS
;; load XMM DR8 source catalog
xmm = mrdfits(pth1,1)
;; use only sources without the possibility of spurious detection (SUM_FLAG == 0 or 1)
;; use only unique sources
xmm = clean_detect_xmm(xmm)

spherematch,ra_inf,dec_inf,xmm.ra,xmm.dec,6./3600.,isamp,imatch,sep_xmm
tags = ['OBS_ID','OBS_CLASS', $
        'PN_FILTER','M1_FILTER','M2_FILTER','PN_SUBMODE','M1_SUBMODE','M2_SUBMODE', $
        'SC_RA','SC_DEC', $
        'EP_ONTIME','PN_ONTIME','M1_ONTIME','M2_ONTIME', $
        'SUM_FLAG','EP_FLAG','PN_FLAG','M1_FLAG','M2_FLAG','TSERIES','SPECTRA','HIGH_BACKGROUND', $
        'SC_EP_1_FLUX','SC_EP_1_FLUX_ERR', $
        'SC_EP_2_FLUX','SC_EP_2_FLUX_ERR', $
        'SC_EP_3_FLUX','SC_EP_3_FLUX_ERR', $
        'SC_EP_4_FLUX','SC_EP_4_FLUX_ERR', $
        'SC_EP_5_FLUX','SC_EP_5_FLUX_ERR', $
        'SC_EP_8_FLUX','SC_EP_8_FLUX_ERR', $
        'SC_EP_9_FLUX','SC_EP_9_FLUX_ERR', $
        'SC_HR1','SC_HR1_ERR','SC_HR2','SC_HR2_ERR','SC_HR3','SC_HR3_ERR','SC_HR4','SC_HR4_ERR' $
        ]
nvars = n_elements(tags)

;; conflicting definition RA/Dec
ipos = where(strmatch(tags,'SC_RA') or strmatch(tags,'SC_DEC'))
pos = (strsplit(tags[ipos],'SC_',/extract,/regex)).ToArray()
xmm_vars = tags
xmm_vars[ipos] = pos+'_XMM'

;; declare and initialize XMM variables, matched to main sample
for i = 0,nvars-1 do begin
    re = execute('type = typename(xmm.'+tags[i]+')')
    re = execute(xmm_vars[i]+' = make_array(nsrc,/'+type+')')
    re = execute(xmm_vars[i]+'[isamp] = xmm[imatch].'+tags[i])
endfor

;; determine source detections (exposure time, flux, flux error)
xband = (strsplit(tags[where(strmatch(tags,'SC_EP_?_FLUX'))],'SC_EP_?_FLUX',/extract)).ToArray()
iix = 'II'+xband
xband = tags[where(strmatch(tags,'SC_EP_?_FLUX'))]
xberr = tags[where(strmatch(tags,'SC_EP_?_FLUX_ERR'))]
for i = 0,n_elements(xband)-1 do begin
    xstr = xband[i]+' gt 0. and '+xberr[i]+' gt 0.'
    re = execute(iix[i]+' = '+xstr)
endfor
iixstr = '('+strjoin(iix,' or ')+') and EP_ONTIME'
;; boolean flag for valid detection in any band
re = execute('iidet_xmm = '+iixstr)
idet_xmm = where(iidet_xmm)


xmm_str = 'XMM,IIDET_XMM,IDET_XMM,'+strjoin(xmm_vars,',')
re = execute('save,'+xmm_str+',/compress,file="xdetect_xmm.sav"')


END




;xmmra = xmm.ra
;xmmde = xmm.dec
;;; ID of observation
;obs_id = long(xmm.obs_id)
;uniqid = obs_id[uniq(obs_id,sort(obs_id))]
;numid = n_elements(uniqid)
;;; array for unique observation field centers
;obsra = dblarr(numid)
;obsde = dblarr(numid)
;
;;; average RA/Dec of all sources in field to estimate field center position
;for i = 0,numid-1 do begin
;    ii = where(uniqid[i] eq obs_id,idlen)
;    thisra = xmmra[ii]
;    thisde = xmmde[ii]
;    obsde[i] = mean(thisde[ii])
;    delmm = max(thisra)-min(thisra)
;    ;; some observations wrap in RA; correct and average
;    if (delmm*cos(obsde[i]*!const.pi/180.) gt 10.) then begin
;        if (median(thisra) lt 180.) then thisra[where(thisra gt 180.)] -= 360. else $
;                                         thisra[where(thisra lt 180.)] += 360.
;    endif
;    obsra[i] = mean(thisra)
;endfor


