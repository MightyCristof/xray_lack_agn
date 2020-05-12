PRO clean_source_xmm, in_data


common _det_xmm

;; clean photometry
;; passes all quality flags
iinoflag_xmm = (SUM_FLAG EQ 0 OR SUM_FLAG EQ 1) and $      ;; no spurious fields
                strmatch(PN_SUBMODE,'*Full*')              ;; full CCD chip readout "SCIENCE MODE"
;; and fail-safe is in XMM catalog
iiclean_xmm = iix_xmm and iinoflag_xmm

;; removed sources
iidirty_xmm = iix_xmm and ~iiclean_xmm
save,iiclean_xmm,iidirty_xmm,file='cleaned_xmm.sav'


END



