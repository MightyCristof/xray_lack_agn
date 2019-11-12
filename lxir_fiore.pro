FUNCTION lxir_fiore, log_lir


;; LX-LIR relationship of Fiore+09
;; log L(5.8 um) < 43.04
;; log L(2 - 10 keV) = log L(5.8 um) - 0.3
;; log L(5.8 um) > 43.04 : 
;; log L(2 - 10 keV) = 43.574 + 0.72(log L(5.8 um) - 44.2)
nobj = n_elements(log_lir)

log_lxir = dblarr(nobj)
iilow = log_lir lt 43.04
log_lxir[where(iilow)] = log_lir[where(iilow)] - 0.3
log_lxir[where(~iilow)] = 43.574 + 0.72*(log_lir[where(~iilow)] - 44.2)

return, log_lxir


END



