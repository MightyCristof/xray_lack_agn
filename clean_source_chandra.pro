FUNCTION clean_source_chandra, in_data


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
;; flux > 0, err > 0
tags = tag_names(data)
iflux = where(strmatch(tags,'FLUX_POWLAW_APER90_?'))
ihilim = where(strmatch(tags,'FLUX_POWLAW_APER90_HILIM_?'))
ilolim = where(strmatch(tags,'FLUX_POWLAW_APER90_LOLIM_?'))
flux_str = strjoin('(DATA.'+tags[iflux]+' gt 0. and DATA.'+tags[ihilim]+' gt 0. and DATA.'+tags[ilolim]+' gt 0.)',' or ')
re = execute('ig = where('+flux_str+')')
data = data[ig]

;; convert 90% aperture flux to 100%
iaper90 = where(strmatch(tags,'*APER90*'),naper)
for i = 0,naper-1 do re = execute('DATA.'+tags[iaper90[i]]+' *= 1.1')

;; sort by exposure time
isort = sort(data.acis_time)
data = data[isort]

return, data


END


