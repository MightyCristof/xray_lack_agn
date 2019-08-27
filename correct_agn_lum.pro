FUNCTION correct_agn_lum, lum_in, $
                          wave_in, $
                          flux_in, $
                          e_flux_in, $
                          param_in, $
                          z_in, $
                          fagn_in, $
                          OVER = over, $
                          UNDER = under, $
                          IICORR = iilum


lum = lum_in
wave = wave_in
flux = flux_in
e_flux = e_flux_in
param = param_in
ebv = reform(param[0,*])
z = z_in
coeff = param[2:5,*]
fagn = fagn_in

;; calculate template flux at WISE 
temp_phot = interp_ir_temp(wave,flux,ebv,z,coeff)
match,wave,[3.4,4.6,12.,22.],iwise,i2
;; count the number of corrected sources
iilum = bytarr(n_elements(z))


if keyword_set(over) then begin
    ;; is the template constrained by the source WISE photometry - overestimate
    data_upper = flux[iwise,*]+e_flux[iwise,*]
    temp_over = temp_phot gt data_upper and data_upper gt 0.

    iiover = total(temp_over[1:2,*],1) ge 2. or total(temp_over[2:3,*],1) ge 2.
    iover = where(iiover and lum gt 0.,ctover)
    if (ctover gt 0) then begin
        flux_over = interp_6micr_agn(wave,flux[*,iover])
        lum_over = interp_6micr_lum(flux_over,z[iover],/log);+alog10(fagn[iover])
        icorr = where(lum_over lt lum[iover],ctcorr)
        if (ctcorr gt 0) then begin
            lum[iover[icorr]] = lum_over[icorr]
            iilum[iover[icorr]] = 1
        endif
    endif
endif

if keyword_set(under) then begin
    ;; is the template constrained by the source WISE photometry - underestimate
    data_lower = flux[iwise,*]-e_flux[iwise,*]
    temp_under = temp_phot lt data_lower and data_lower gt 0.

    iiunder = total(temp_under[1:2,*],1) ge 2. or total(temp_under[2:3,*],1) ge 2.
    iunder = where(iiunder and lum gt 0.,ctunder)
    if (ctunder gt 0) then begin
        flux_under = interp_6micr_agn(wave,flux[*,iunder])
        lum_under = interp_6micr_lum(flux_under,z[iunder],/log);+alog10(fagn[iunder])
        icorr = where(lum_under gt lum[iunder],ctcorr)
        if (ctcorr gt 0) then begin
            lum[iunder[icorr]] = lum_under[icorr]
            iilum[iunder[icorr]] = 1
        endif
    endif
endif

  
return, lum


END

