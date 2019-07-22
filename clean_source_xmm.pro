PRO clean_source_xmm, in_data


common _det_xmm

vars = scope_varname(common = '_det_xmm')


;; exposure > 0
iitime = PN_ONTIME gt 0.

;; photometry exists in at least one band
xbflx = vars[where(strmatch(vars,'PN_?_FLUX'))]
xberr = vars[where(strmatch(vars,'PN_?_FLUX_ERR'))]
phot_str = strjoin('('+xbflx+' gt 0. and '+xberr+' gt 0.)',' or ')
re = execute('iiphot = '+phot_str)

;; passes all quality flags
iipass = (SUM_FLAG EQ 0 OR SUM_FLAG EQ 1) and $      ;; no spurious fields
          strmatch(PN_SUBMODE,'*Full*')              ;; full CCD chip readout "SCIENCE MODE"

;; unique source (only 17 sources)
;iu = uniq(SRCID,sort(SRCID))


;; clean photometry
iiclean_xmm = iitime and iiphot and iipass
;; removed sources
iiflag_xmm = iidet_xmm and ~iiclean_xmm
save,iiclean_xmm,iiflag_xmm,file='cleaned_xmm.sav'


END



