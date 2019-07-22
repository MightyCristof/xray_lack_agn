PRO clean_source_nustar


common _det_nst

vars = scope_varname(common='_det_nst')

;; exposure > 0
xbexp = vars[where(strmatch(vars,'?EXP'))]
time_str = strjoin(+xbexp+' gt 0.',' or ')
re = execute('iitime = '+time_str)

;; photometry exists in at least one band
xbflx = vars[where(strmatch(vars,'?BF'))]
xberr = vars[where(strmatch(vars,'E_?BF'))]
phot_str = strjoin('('+xbflx+' gt 0. and '+xberr+' gt 0.)',' or ')
re = execute('iiphot = '+phot_str)

;; passes all quality flags
;iipass = 

;; clean photometry
iiclean_nst = iitime and iiphot; and iipass
;; removed sources
iiflag_nst = iidet_nst and ~iiclean_nst
save,iiclean_nst,iiflag_nst,file='cleaned_nst.sav'


END





