FUNCTION clean_detect_xmm, in_data


data = in_data

;; no spurious fields
ig = where(data.sum_flag eq 0 or data.sum_flag eq 1,/null)
data = data[ig]
;; unique source
ig = uniq(data.srcid,sort(data.srcid))
data = data[ig]
;; exposure > 0
ig = where(data.ep_ontime gt 0.)
data = data[ig]
;; full CCD chip readout - SCIENCE mode
tags = tag_names(data)
imode = where(strmatch(tags,'*SUBMODE*'))
mode_str = strjoin('strmatch(DATA.'+tags[imode]+',"*Full*")',' and ')
re = execute('ig = where('+mode_str+')')
data = data[ig]
;; flux > 0
iflux = where(strmatch(tags,'EP_?_FLUX'))
flux_str = strjoin('DATA.'+tags[iflux]+' gt 0.',' or ')
re = execute('ig = where('+flux_str+')')
data = data[ig]
;; err > 0
ierr = where(strmatch(tags,'EP_?_FLUX_ERR'))
err_str = strjoin('DATA.'+tags[ierr]+' gt 0.',' or ')
re = execute('ig = where('+err_str+')')
data = data[ig]
;; SC flux > 0
iscflux = where(strmatch(tags,'SC_EP_?_FLUX'))
scflux_str = strjoin('DATA.'+tags[iscflux]+' gt 0.',' or ')
re = execute('ig = where('+scflux_str+')')
data = data[ig]
;; SC err > 0
iscerr = where(strmatch(tags,'SC_EP_?_FLUX_ERR'))
scerr_str = strjoin('DATA.'+tags[iscerr]+' gt 0.',' or ')
re = execute('ig = where('+scerr_str+')')
data = data[ig]


;ig = where()
;; exposure > 0
;ig = where(data.ep_ontime gt 0.,/null)
;data = data[ig]
;ig = where((data.sc_ep_4_flux gt 0. and data.sc_ep_4_flux_err gt 0.) or (data.sc_ep_5_flux gt 0. and data.sc_ep_5_flux_err gt 0.),/null)
;data = data[ig]
;ig = where((data.sc_ep_4_flux/data.sc_ep_4_flux_err gt 2.) or (data.sc_ep_5_flux/data.sc_ep_5_flux_err gt 2.),/null)
;data = data[ig]

;; sort by exposure time
isort = sort(data.ep_ontime)
data = data[isort]

;; create EP_4+EP_5 channel
;inon = where(data.sc_ep_4_flux eq 0. or data.sc_ep_5_flux eq 0.)
;new_epf = data.sc_ep_4_flux + data.sc_ep_5_flux
;new_epe = sqrt(data.sc_ep_4_flux_err^2 + data.sc_ep_5_flux_err^2)
;new_epf[inon] = 0.
;new_epe[inon] = 0.
;struct_add_field,data,'SC_EP_45_FLUX',new_epf
;struct_add_field,data,'SC_EP_45_FLUX_ERR',new_epe

return, data


END





