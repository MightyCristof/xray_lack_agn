FUNCTION lxir_c17, log_lir, $
                   SCATTER = scat


;; LX-LIR relationship of Chen+17
;; log L6um < 44.79 :
;; log LX = 0.84±0.03 * log(L6um/10^45) + 44.60±0.01
;; log L6um ³ 44.79 :
;; log LX = 0.40±0.03 * log(L6um/10^45) + 44.51±0.01
nobj = n_elements(log_lir)
sig = 0.3

log_lxir = dblarr(nobj)
scat = dblarr(nobj)
iilo = log_lir lt 44.79
ilo = where(iilo,loct)
ihi = where(~iilo,hict)
if (loct gt 0.) then begin
    log_lxir[ilo] = (0.84)*(log_lir[ilo]-45.) + (44.60)
    scat[ilo] = randomn(seed,loct)*sig
endif
if (hict gt 0.) then begin
    log_lxir[ihi] = (0.40)*(log_lir[ihi]-45.) + (44.51)
    scat[ihi] = randomn(seed,hict)*sig
endif

return, log_lxir


END






