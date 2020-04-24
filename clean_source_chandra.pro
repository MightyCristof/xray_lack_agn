PRO clean_source_chandra


common _det_cha

;; clean photometry
;; passes all quality flags or was manually added after review
iiclean_cha = (DITHER_WARNING_FLAG eq 'F' and PILEUP_FLAG eq 'F' and SAT_SRC_FLAG eq 'F' and  $
               VAR_FLAG eq 'F' and STREAK_SRC_FLAG eq 'F' and VAR_INTER_HARD_FLAG eq 'F' ) $
               OR MAN_ADD_FLAG eq 'T'

;; removed sources
iiflag_cha = iidet_cha and ~iiclean_cha
save,iiclean_cha,iiflag_cha,file='cleaned_cha.sav'


END




