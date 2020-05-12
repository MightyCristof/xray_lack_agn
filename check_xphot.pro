PRO check_xphot, ind, $
                 XCHA = xcha, $
                 XXMM = xxmm, $
                 XNST = xnst
                      
                      
common _det_cha
common _det_xmm
common _det_nst

for i = 0,n_elements(ind)-1 do begin
    if keyword_set(xcha) then begin
        print, 'EXP TIME:   '+strtrim(acis_time[ind[i]],2)
        print, '  H BAND:   '+strtrim(flux_powlaw_aper90_h[ind[i]],2)+', '+strtrim(flux_powlaw_aper90_h_err[ind[i]],2)
        print, '  B BAND:   '+strtrim(flux_powlaw_aper90_b[ind[i]],2)+', '+strtrim(flux_powlaw_aper90_b_err[ind[i]],2)
        print, '  W BAND:   '+strtrim(flux_powlaw_aper90_w[ind[i]],2)+', '+strtrim(flux_powlaw_aper90_w_err[ind[i]],2)
        print, '  S BAND:   '+strtrim(flux_powlaw_aper90_s[ind[i]],2)+', '+strtrim(flux_powlaw_aper90_s_err[ind[i]],2)
        print, '  M BAND:   '+strtrim(flux_powlaw_aper90_m[ind[i]],2)+', '+strtrim(flux_powlaw_aper90_m_err[ind[i]],2)
        print, '  U BAND:   '+strtrim(flux_powlaw_aper90_u[ind[i]],2)+', '+strtrim(flux_powlaw_aper90_u_err[ind[i]],2)
    endif

    if keyword_set(xxmm) then begin
        print, 'EXP TIME:   '+strtrim(pn_ontime[ind[i]],2)
        print, '  4 BAND:   '+strtrim(pn_4_flux[ind[i]],2)+', '+strtrim(pn_4_flux_err[ind[i]],2)
        print, '  5 BAND:   '+strtrim(pn_5_flux[ind[i]],2)+', '+strtrim(pn_5_flux_err[ind[i]],2)
        print, '  3 BAND:   '+strtrim(pn_3_flux[ind[i]],2)+', '+strtrim(pn_3_flux_err[ind[i]],2)
        print, '  9 BAND:   '+strtrim(pn_9_flux[ind[i]],2)+', '+strtrim(pn_9_flux_err[ind[i]],2)
        print, '  2 BAND:   '+strtrim(pn_2_flux[ind[i]],2)+', '+strtrim(pn_2_flux_err[ind[i]],2)
        print, '  8 BAND:   '+strtrim(pn_8_flux[ind[i]],2)+', '+strtrim(pn_8_flux_err[ind[i]],2)
    endif

    if keyword_set(xnst) then begin
        print, 'S EXP:   '+strtrim(sexp[ind[i]],2)+'      S BAND:   '+strtrim(sbf[ind[i]],2)+', '+strtrim(e_sbf[ind[i]],2)
        print, 'H EXP:   '+strtrim(hexp[ind[i]],2)+'      H BAND:   '+strtrim(hbf[ind[i]],2)+', '+strtrim(e_hbf[ind[i]],2)
        print, 'F EXP:   '+strtrim(fexp[ind[i]],2)+'      F BAND:   '+strtrim(fbf[ind[i]],2)+', '+strtrim(e_fbf[ind[i]],2)
    endif
endfor


END