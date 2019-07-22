FUNCTION interp_6micr_agn, wave_in, $
                           flux_in
                           ;ebv_in, $
                           ;z_in, $
                           ;c_agn_in, $
                           ;DERED = dered


isort = sort(wave_in)
wave = wave_in[isort]
flux = flux_in[isort,*]

nsrc = n_elements(flux[0,*])
sixm_agn = dblarr(nsrc)

;nu6 = (!const.c*1e6)/6.
;for i = 0,nsrc-1 do sixm_agn[i] = 1e-29*nu6*interpol(flux[*,i],wave,6.)
nu = (!const.c*1e6)/wave
for i = 0,nsrc-1 do sixm_agn[i] = interpol(1e-29*nu*flux[*,i],wave,6.)

return, sixm_agn


;common _comp
;
;isort = sort(wave_in)
;wave = wave_in[isort]
;flux = flux_in[isort,*]
;c_agn = reform(c_agn_in)
;if keyword_set(dered) then ebv = reform(ebv_in) else ebv = dblarr(n_elements(ebv_in))
;z = reform(z_in)
;fagn = reform(fagn_in)
;
;nsrc = n_elements(z)
;
;wave_phot = wave
;wave_temp = comp.wav#(1+z)
;nu_phot = (!const.c*1e6)/wave_phot
;nu_temp = (!const.c*1e6)/wave_temp
;flux_phot = 1e-29*rebin(nu_phot,n_elements(wave_in),n_elements(z))*flux
;flux_agn = 1e-29*nu_temp*(comp.agn#c_agn*10.^(-0.4*comp.kap#ebv))
;
;sixm_agn = dblarr(nsrc)
;
;for i = 0,nsrc-1 do sixm_agn[i] = 1e-29*nu6*interpol(flux[*,i],wave,6.)*fagn[0]
;
;for i = 0,nsrc-1 do begin
;    w1_temp = interpol(flux_temp[*,i],wave_temp[*,i],3.4)
;    w2_temp = interpol(flux_temp[*,i],wave_temp[*,i],4.6)
;    w3_temp = interpol(flux_temp[*,i],wave_temp[*,i],12.)
;    w4_temp = interpol(flux_temp[*,i],wave_temp[*,i],22.)
;    wise_temp[*,i] = [w1_temp,w2_temp,w3_temp,w4_temp]
;
;    wise_temp[*,i] *= 1e29/nu_phot[11:14]
;endfor
;
;;return, wise_temp

END


