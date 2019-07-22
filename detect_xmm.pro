PRO detect_xmm


common _fits
nsrc = n_elements(ra)

;; XMM DR8 Source Catalog
pth1 = '/Users/ccarroll/Research/surveys/XMM/3XMM_DR8cat_v1.0.fits'

;; DETECTIONS
;; load XMM DR8 source catalog
xmm = mrdfits(pth1,1)

spherematch,ra,dec,xmm.ra,xmm.dec,6./3600.,isamp,imatch,sep_xmm
tags = ['OBS_ID','OBS_CLASS', $
        'PN_FILTER','PN_SUBMODE', $
        'SC_RA','SC_DEC', $
        'EP_ONTIME','PN_ONTIME', $
        'SUM_FLAG','EP_FLAG','PN_FLAG','TSERIES','SPECTRA','HIGH_BACKGROUND', $
        'PN_1_FLUX','PN_1_FLUX_ERR', $
        'PN_2_FLUX','PN_2_FLUX_ERR', $
        'PN_3_FLUX','PN_3_FLUX_ERR', $
        'PN_4_FLUX','PN_4_FLUX_ERR', $
        'PN_5_FLUX','PN_5_FLUX_ERR', $
        'PN_8_FLUX','PN_8_FLUX_ERR', $
        'PN_9_FLUX','PN_9_FLUX_ERR', $
        'SC_HR1','SC_HR1_ERR','SC_HR2','SC_HR2_ERR','SC_HR3','SC_HR3_ERR','SC_HR4','SC_HR4_ERR', $
        'SUM_FLAG','SRCID']
nvars = n_elements(tags)

;; conflicting definition RA/Dec
ipos = where(strmatch(tags,'SC_RA') or strmatch(tags,'SC_DEC'))
pos = (strsplit(tags[ipos],'SC_',/extract,/regex)).ToArray()
xmm_vars = tags
xmm_vars[ipos] = pos+'_XMM'

;; declare and initialize XMM variables, matched to main sample
for i = 0,nvars-1 do begin
    re = execute('type = typename(xmm.'+tags[i]+')')
    if (type eq 'LONG64') then type = 'L64'
    re = execute(xmm_vars[i]+' = make_array(nsrc,/'+type+')')
    re = execute(xmm_vars[i]+'[isamp] = xmm[imatch].'+tags[i])
endfor

;; boolean flag for valid detection in any band
iidet_xmm = bytarr(nsrc)
iidet_xmm[isamp] = 1
idet_xmm = where(iidet_xmm)

;; save data
xmm_str = 'XMM,IIDET_XMM,IDET_XMM,'+strjoin(xmm_vars,',')
re = execute('save,'+xmm_str+',/compress,file="detections_xmm.sav"')


END




