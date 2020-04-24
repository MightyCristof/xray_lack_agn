PRO clean_source_xmm, in_data


common _det_xmm

;; clean photometry
;; passes all quality flags
iiclean_xmm = (SUM_FLAG EQ 0 OR SUM_FLAG EQ 1) and $      ;; no spurious fields
               strmatch(PN_SUBMODE,'*Full*')              ;; full CCD chip readout "SCIENCE MODE"

;; removed sources
iiflag_xmm = iidet_xmm and ~iiclean_xmm
save,iiclean_xmm,iiflag_xmm,file='cleaned_xmm.sav'


END



