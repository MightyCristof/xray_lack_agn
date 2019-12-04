FUNCTION lxir_fiore, log_lir, $
                     SCATTER = scatter


;; LX-LIR relationship of Fiore+09
;; log L(5.8 um) < 43.04
;; log L(2 - 10 keV) = log L(5.8 um) - 0.3
;; log L(5.8 um) > 43.04 : 
;; log L(2 - 10 keV) = 43.574 + 0.72(log L(5.8 um) - 44.2)
nobj = n_elements(log_lir)
sig = 0.3

log_lxir = dblarr(nobj)
iilo = log_lir lt 43.04
ilo = where(iilo,loct)
ihi = where(~iilo,hict)
if (loct gt 0.) then log_lxir[ilo] = log_lir[ilo] - 0.3
if keyword_set(scatter) then log_lxir[ilo] += randomn(seed,loct)*sig
if (hict gt 0.) then log_lxir[ihi] = 43.574 + 0.72*(log_lir[ihi] - 44.2)
if keyword_set(scatter) then log_lxir[ihi] += randomn(seed,hict)*sig

return, log_lxir


END



