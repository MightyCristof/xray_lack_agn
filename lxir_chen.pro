FUNCTION lxir_chen, log_lir, $
                    SCATTER = scatter


;; LX-LIR relationship of Chen+17
;; log L6um < 44.79 :
;; log LX = 0.84±0.03 * log(L6um/10^45) + 44.60±0.01
;; log L6um ³ 44.79 :
;; log LX = 0.40±0.03 * log(L6um/10^45) + 44.51±0.01
nobj = n_elements(log_lir)
sig = 0.3

log_lxir = dblarr(nobj)
iilo = log_lir lt 44.79
ilo = where(iilo,loct)
ihi = where(~iilo,hict)
if (loct gt 0.) then log_lxir[ilo] = (0.84)*(log_lir[ilo]-45.) + (44.60)
if keyword_set(scatter) then log_lxir[ilo] += randomn(seed,loct)*sig
;if keyword_set(scatter) then log_lxir[ilo] += 2.*sig*(randomu(seed,loct)-0.5)
if (hict gt 0.) then log_lxir[ihi] = (0.40)*(log_lir[ihi]-45.) + (44.51)
if keyword_set(scatter) then log_lxir[ihi] += randomn(seed,hict)*sig
;if keyword_set(scatter) then log_lxir[ihi] += 2.*sig*(randomu(seed,hict)-0.5)

return, log_lxir


END



