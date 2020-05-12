PRO clean_source_nustar


common _det_nst

;; clean photometry
iiclean_nst = iidet_nst
;; removed sources
iidirty_nst = iix_nst and ~iiclean_nst
save,iiclean_nst,iidirty_nst,file='cleaned_nst.sav'


END





