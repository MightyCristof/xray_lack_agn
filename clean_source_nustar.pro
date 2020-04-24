PRO clean_source_nustar


common _det_nst

;; clean photometry
iiclean_nst = iidet_nst
;; removed sources
iiflag_nst = iidet_nst and ~iiclean_nst
save,iiclean_nst,iiflag_nst,file='cleaned_nst.sav'


END





