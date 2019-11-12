FUNCTION interp_agn6micr, ebv_in, $
                          z_in, $
                          coeff_in, $
                          DERED = dered


common _comp


nsrc = n_elements(z_in)

if keyword_set(dered) then ebv = dblarr(nsrc) else $
                           ebv = reform(ebv_in)
z = reform(z_in)
c_agn = reform(coeff_in)
;c_ell = reform(coeff_in[1,*])
;c_sfg = reform(coeff_in[2,*])
;c_irr = reform(coeff_in[3,*])


wave_temp = comp.wav#(1+z)
nu_temp = !const.c/(wave_temp/1e6)
flux_agn = 1e-29*c_agn##comp.agn * 10.^(-0.4*comp.kap#ebv)

wise6agn = dblarr(nsrc)
for i = 0,nsrc-1 do wise6agn[i] = interpol(flux_agn[*,i]*nu_temp[*,i],wave_temp[*,i],6.)

return, wise6agn

END




