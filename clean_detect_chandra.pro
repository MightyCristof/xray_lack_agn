FUNCTION clean_detect_chandra, in_data


data = in_data

;; passes all quality flags or was manually added after review
;data.CONF_FLAG eq 'F' and $
ig = where((data.DITHER_WARNING_FLAG eq 'F' and $
            data.PILEUP_FLAG eq 'F' and $
            data.SAT_SRC_FLAG eq 'F' and $
            data.VAR_FLAG eq 'F' and $
            data.STREAK_SRC_FLAG eq 'F' and $
            data.VAR_INTER_HARD_FLAG eq 'F' ) $
            OR data.MAN_ADD_FLAG eq 'T')
data = data[ig]

;; exposure > 0
ig = where(data.acis_time gt 0.,/null)
data = data[ig]
;; flux > 0
tags = tag_names(data)
iflux = where(strmatch(tags,'FLUX_POWLAW_APER90_?'))
flux_str = strjoin('DATA.'+tags[iflux]+' gt 0.',' or ')
re = execute('ig = where('+flux_str+')')
data = data[ig]
;; err > 0
ihilim = where(strmatch(tags,'FLUX_POWLAW_APER90_HILIM_?'))
ilolim = where(strmatch(tags,'FLUX_POWLAW_APER90_LOLIM_?'))
islim = '(DATA.'+tags[ihilim]+' gt 0. and DATA.'+tags[ilolim]+' gt 0.)'
lim_str = strjoin(islim,' or ')
re = execute('ig = where('+lim_str+')')
data = data[ig]
;ig = where(data.flux_powlaw_aper90_h/(data.flux_powlaw_aper90_hilim_h-data.flux_powlaw_aper90_h) gt 2.)
;data = data[ig]

;; sort by exposure time
isort = sort(data.acis_time)
data = data[isort]

;; create flux error
;new_err = data.flux_powlaw_aper90_hilim_h - data.flux_powlaw_aper90_h
;struct_add_field,data,'FLUX_POWLAW_APER90_H_ERR',new_err

return, data


END


