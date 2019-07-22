PRO clean_source_chandra


common _det_cha

vars = scope_varname(common = '_det_cha')


;; exposure time > 0 
iitime = ACIS_TIME gt 0.

;; photometry exists in at least one band
xbflx = vars[where(strmatch(vars,'FLUX_POWLAW_APER90_?'))]
xberr = vars[where(strmatch(vars,'FLUX_POWLAW_APER90_?_ERR'))]
phot_str = strjoin('('+xbflx+' gt 0. and '+xberr+' gt 0.)',' or ')
re = execute('iiphot = '+phot_str)

;; passes all quality flags or was manually added after review
iipass = (DITHER_WARNING_FLAG eq 'F' and PILEUP_FLAG eq 'F' and SAT_SRC_FLAG eq 'F' and  $
          VAR_FLAG eq 'F' and STREAK_SRC_FLAG eq 'F' and VAR_INTER_HARD_FLAG eq 'F' ) $
          OR MAN_ADD_FLAG eq 'T'


;; clean photometry
iiclean_cha = iitime and iiphot and iipass
;; removed sources
iiflag_cha = iidet_cha and  ~iiclean_cha
save,iiclean_cha,iiflag_cha,file='cleaned_cha.sav'


END




