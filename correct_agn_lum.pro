FUNCTION correct_agn_lum, lum_in, $
                          wave_in, $
                          flux_in, $
                          e_flux_in, $
                          param_in, $
                          z_in, $
                          fagn_in, $
                          OVER = over, $
                          UNDER = under, $
                          NCORR = ctcorr


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
ctcorr = 0L

if keyword_set(over) then begin
;; is the template constrained by the source WISE photometry - overestimate
data_upper = flux[iwise,*]+e_flux[iwise,*]
delta_over = temp_phot gt data_upper
;iiover = (total(data_upper[2:3,*] gt 0. and delta_over[2:3,*],1) eq 2) OR $         ;; W3+W4 exist and template overshoots W3+W4
;         (reform(data_upper[3,*]) eq 0. and total(delta_over[1:2,*],1) eq 2) OR $   ;; W4 missing and template overshoots W2+W3
;         (reform(data_upper[2,*]) eq 0. and total(delta_over[[1,3],*],1) eq 2)      ;; W3 missing and template overshoots W2+W4 (found none)
;iiover = total(data_upper gt 0. and delta_over,1) ge 2.
iiover = total(data_upper[2:3,*] gt 0. and delta_over[2:3,*],1) ge 2.
iover = where(iiover)
flux_over = interp_6micr_agn(wave,flux[*,iover])
lum_over = interp_6micr_lum(flux_over,z[iover],/log)+alog10(fagn[iover])
icorr = where(lum_over lt lum[iover],ctover)
if (ctover gt 0) then lum[iover[icorr]] = lum_over[icorr]
;lum[iover] = lum_over
;print, ctover
ctcorr += ctover
endif

if keyword_set(under) then begin
;; is the template constrained by the source WISE photometry - underestimate
data_lower = flux[iwise,*]-e_flux[iwise,*]
delta_under = temp_phot lt data_lower
iiunder = total(data_lower[2:3,*] gt 0. and delta_under[2:3,*],1) ge 2.
iunder = where(iiunder)

flux_under = interp_6micr_agn(wave,flux[*,iunder])
lum_under = interp_6micr_lum(flux_under,z[iunder],/log)+alog10(fagn[iunder])
icorr = where(lum_under gt lum[iunder],ctunder)
if (ctunder gt 0) then lum[iunder[icorr]] = lum_under[icorr]
;lum[iunder] = lum_under
;print, ctunder
ctcorr += ctunder
endif

  
return, lum


END

