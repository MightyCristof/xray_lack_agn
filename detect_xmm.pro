PRO detect_xmm


common _fits
common _inf_xmm

nsrc = n_elements(ra)

;; XMM DR8 Source Catalog
pth1 = '/Users/ccarroll/Research/surveys/XMM/3XMM_DR8cat_v1.0.fits'

;; DETECTIONS
;; load XMM DR8 source catalog
xmm = mrdfits(pth1,1)

;; XMM PSF FWHM = 12.5"
;; https://heasarc.nasa.gov/docs/xmm/uhb/onaxisxraypsf.html#uhb:fig:onaxis_psf
psf_xmm = 6.25
spherematch,ra,dec,xmm.ra,xmm.dec,psf_xmm/3600.,isamp,imatch,sepx
sepx *= 3600.

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
;; sample separation
sep_xmm = dblarr(nsrc)
sep_xmm[isamp] = sepx

;; boolean flag for cross-matched with 3XMM
iix_xmm = bytarr(nsrc)
iix_xmm[isamp] = 1
;; ensure valid photometry
phot = tags[where(strmatch(xmm_vars,'PN_?_FLUX'),nphot)]
photerr = tags[where(strmatch(xmm_vars,'PN_?_FLUX_ERR'),nphoterr)]
re = execute('iiphot = '+strjoin("(finite("+phot+") and "+phot+" gt 0. and finite("+photerr+") and "+photerr+" gt 0.)"," or "))
;; ensure valid exposure time
iitime = finite(pn_ontime) and pn_ontime gt 0.
;; boolean flag for valid detections in 3XMM
iidet_xmm = iiphot and iitime
;; boolean flag for infield non-detections
iinon_xmm = iiinf_xmm and ~iix_xmm

;; "clean" X-ray observations
;; passes all quality flags
iinoflag_xmm = (SUM_FLAG EQ 0 OR SUM_FLAG EQ 1) and $      ;; no spurious fields
                strmatch(PN_SUBMODE,'*Full*')              ;; full CCD chip readout "SCIENCE MODE"
;; and fail-safe is in XMM catalog
iiclean_xmm = iix_xmm and iinoflag_xmm
;; removed sources
iidirty_xmm = iix_xmm and ~iiclean_xmm
;; unflag detections where X-ray observations are not clean
iidet_xmm[where(iidirty_xmm,/null)] = 0

;; save detection data
xmm_str = 'XMM,SEP_XMM,IIX_XMM,IIDET_XMM,IINON_XMM,IICLEAN_XMM,IIDIRTY_XMM,'+strjoin(xmm_vars,',')
re = execute('save,'+xmm_str+',/compress,file="detections_xmm.sav"')

;; update in-field data
;inew = where(iix_xmm eq 1 and iiinf_xmm eq 0,ctnew)
;if (ctnew gt 0) then begin
;    iiinf_xmm[inew] = 1
;    texp_xmm[inew] = -1.
;    sdst_xmm[inew] = -1.
;    inf_str = strjoin(scope_varname(common='_INF_XMM'),',')
;    re = execute('save,'+inf_str+',file="infield_xmm.sav"')
;endif


END




