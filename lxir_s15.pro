FUNCTION lxir_s15, log_lir, $
                   SCATTER = scat


;; LX-LIR relationship of Stern+15

;; log (2 10 keV) = 40.981 + 1.024x - 0.047x^2
;; where x = log(vLv 6um) - 41.

nobj = n_elements(log_lir)
sig = 0.3

log_lxir = dblarr(nobj)
scat = dblarr(nobj)
log_lxir[*] = 40.981 + 1.024*(log_lir-41.) - 0.047*(log_lir-41.)^2.
scat[*] = randomn(seed,nobj)*sig

return, log_lxir


END



