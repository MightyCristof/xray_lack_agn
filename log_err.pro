FUNCTION log_err, lin_flx, $
                  lin_err
                  

;; output array
log_err = dblarr(n_elements(lin_err))
;; where data exists
ival = where(lin_flx gt 0. and lin_err gt 0.,ctval)         
if (ctval eq 0.) then begin
    print, 'NO VALID INPUT'
    return, 0
endif
log_err[ival] = sqrt((lin_err[ival]/(lin_flx[ival]*alog(10.)))^2.)

return, log_err


END





