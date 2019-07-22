FUNCTION interp_ir_temp, wave_in, $
                         flux_in, $
                         ebv_in, $
                         z_in, $
                         coeff_in


common _comp


isort = sort(wave_in)
wave = wave_in[isort]
flux = flux_in[isort,*]
c_agn = reform(coeff_in[0,*])
c_ell = reform(coeff_in[1,*])
c_sfg = reform(coeff_in[2,*])
c_irr = reform(coeff_in[3,*])
ebv = reform(ebv_in)
z = reform(z_in)

nsrc = n_elements(z)

wave_phot = wave
wave_temp = comp.wav#(1+z)
nu_phot = (!const.c*1e6)/wave_phot
nu_temp = (!const.c*1e6)/wave_temp
flux_phot = 1e-29*rebin(nu_phot,n_elements(wave_in),n_elements(z))*flux
flux_agn = 1e-29*nu_temp*(comp.agn#c_agn*10.^(-0.4*comp.kap#ebv))
flux_temp = 1e-29*nu_temp*(comp.ell#c_ell + comp.sfg#c_sfg + comp.irr#c_irr) + flux_agn

wise_temp = dblarr(4,nsrc)
for i = 0,nsrc-1 do begin
    w1_temp = interpol(flux_temp[*,i],wave_temp[*,i],3.4)
    w2_temp = interpol(flux_temp[*,i],wave_temp[*,i],4.6)
    w3_temp = interpol(flux_temp[*,i],wave_temp[*,i],12.)
    w4_temp = interpol(flux_temp[*,i],wave_temp[*,i],22.)
    wise_temp[*,i] = [w1_temp,w2_temp,w3_temp,w4_temp]

    wise_temp[*,i] *= 1e29/nu_phot[11:14]
endfor

return, wise_temp

END




