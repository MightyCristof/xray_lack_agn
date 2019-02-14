PRO convert_xray05_2kev, PLT = plt


;; load variables
common _inf_nst
common _inf_xmm
common _inf_cha
common _det_nst
common _det_xmm
common _det_cha
common _soft210


;; Instrument variables: exposure time, flux, error
ff_052 = ff+'_052'
ee_052 = ee+'_052'

;; Flux conversions
;; WebPIMMS parameters: Galactic NH=2E20, power law photon index=1.8
;; https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3pimms/w3pimms.pl
case phot_ind of
    2.0: begin
        ;; NuSTAR                   ;; energy band (keV)
        out052_nst_s = 1.328E+00    ;; 3-8
        out052_nst_h = 0.    ;; 8-24
        out052_nst_f = 0.    ;; 3-24
        ;; XMM
        out052_xmm_1 = 0.    ;; 0.2-0.5
        out052_xmm_2 = 2.087E+00    ;; 0.5-1
        out052_xmm_3 = 0.    ;; 1-2
        out052_xmm_4 = 0.    ;; 2-4.5
        out052_xmm_5 = 0.    ;; 4.5-12
        out052_xmm_8 = 0.    ;; 0.2-12
        out052_xmm_9 = 0.    ;; 0.5-4.5
        ;; Chandra
        out052_cha_b = 0.    ;; 0.5-7
        out052_cha_h = 0.    ;; 2-7
        out052_cha_m = 2.590E+00    ;; 1.2-2
        out052_cha_s = 1.629E+00    ;; 0.5-1.2
        out052_cha_u = 0.    ;; 0.2-0.5
        out052_cha_w = 0.    ;; 0.1-10
        end
    1.8: begin
        ;; NuSTAR                   ;; energy band (keV)
        out052_nst_s = 9.716E-01    ;; 3-8
        out052_nst_h = 0.    ;; 8-24
        out052_nst_f = 0.    ;; 3-24
        ;; XMM
        out052_xmm_1 = 0.    ;; 0.2-0.5
        out052_xmm_2 = 2.248E+00    ;; 0.5-1
        out052_xmm_3 = 0.    ;; 1-2
        out052_xmm_4 = 0.    ;; 2-4.5
        out052_xmm_5 = 0.    ;; 4.5-12
        out052_xmm_8 = 0.    ;; 0.2-12
        out052_xmm_9 = 0.    ;; 0.5-4.5
        ;; Chandra
        out052_cha_b = 0.    ;; 0.5-7
        out052_cha_h = 0.    ;; 2-7
        out052_cha_m = 2.388E+00    ;; 1.2-2
        out052_cha_s = 1.720E+00    ;; 0.5-1.2
        out052_cha_u = 0.    ;; 0.2-0.5
        out052_cha_w = 0.    ;; 0.1-10
        end
    1.4: begin
        ;; NuSTAR                   ;; energy band (keV)
        out052_nst_s = 5.249E-01    ;; 3-8
        out052_nst_h = 0.    ;; 8-24
        out052_nst_f = 0.    ;; 3-24
        ;; XMM
        out052_xmm_1 = 0.    ;; 0.2-0.5
        out052_xmm_2 = 2.643E+00    ;; 0.5-1
        out052_xmm_3 = 0.    ;; 1-2
        out052_xmm_4 = 0.    ;; 2-4.5
        out052_xmm_5 = 0.    ;; 4.5-12
        out052_xmm_8 = 0.    ;; 0.2-12
        out052_xmm_9 = 0.    ;; 0.5-4.5
        ;; Chandra
        out052_cha_b = 0.    ;; 0.5-7
        out052_cha_h = 0.    ;; 2-7
        out052_cha_m = 2.064E+00    ;; 1.2-2
        out052_cha_s = 1.940E+00    ;; 0.5-1.2
        out052_cha_u = 0.    ;; 0.2-0.5
        out052_cha_w = 0.    ;; 0.1-10
        ;; 2-10 keV to 0.5-2 keV
        ;out052_in210 = 3.300E-01
        end
    else: print, 'GAMMA NOT IN RANGE'
endcase
;; energy conversions
xconv = ['OUT052_NST_S','OUT052_NST_H','OUT052_NST_F', $
         'OUT052_XMM_1','OUT052_XMM_2','OUT052_XMM_3','OUT052_XMM_4','OUT052_XMM_5','OUT052_XMM_8','OUT052_XMM_9', $
         'OUT052_CHA_B','OUT052_CHA_H','OUT052_CHA_M','OUT052_CHA_S','OUT052_CHA_U','OUT052_CHA_W']

;; structure containing flux conversions for all instruments/energy bands
gmma = 'GAMMA:'+string(phot_ind,"(d3.1)")
ff_str = '"'+ff+'",0d'
re = execute('cnv_052 = create_struct(name=gmma,'+strjoin(ff_str,",")+')')
;; fill conversion structure by photon index
;; apply conversion to fluxes and errors
for i = 0,nxband-1 do begin
    re = execute('cnv_052.(i) = '+xconv[i])
    re = execute(ff_052[i]+' = '+ff[i]+' * cnv_052.(i)')
    re = execute(ee_052[i]+' = '+ee[i]+' * cnv_052.(i)')
endfor

sav_vars = ['FF_052','EE_052',ff_052,ee_052, $
            'CNV_052']
sav_inds = []


;; commit fluxes and errors to source for use
;; currently choosing one band per telescope
xray_exp_052 = 'EXP'+xfield+'_052'
xray_flx_052 = 'FLX'+xfield+'_052'
xray_err_052 = 'ERR'+xfield+'_052'
xray_cnv_052 = 'CNV'+xfield+'_052'

ixband_052 = [0,4,13]
used_exp = tt[ixband_052]
used_flx = ff_052[ixband_052]
used_err = ee_052[ixband_052]
used_cnv = 'CNV_052.'+(tag_names(cnv_052))[ixband_052]

;; X-ray detections with chosen band
iidet_052 = iidet+'_052'
idet_052 = idet+'_052'

for i = 0,nfield-1 do begin
    re = execute(xray_exp_052[i]+' = '+used_exp[i])
    re = execute(xray_flx_052[i]+' = '+used_flx[i])
    re = execute(xray_err_052[i]+' = '+used_err[i])
    re = execute(xray_cnv_052[i]+' = '+used_cnv[i])
    re = execute(iidet_052[i]+' = '+xray_exp_052[i]+' gt 0. and '+xray_flx_052[i]+' gt 0. and '+xray_err_052[i]+' gt 0.')
    re = execute(idet_052[i]+' = where('+iidet_052[i]+')')
endfor

sav_vars = [sav_vars,'XRAY_EXP_052','XRAY_FLX_052','XRAY_ERR_052','XRAY_CNV_052', $
                     xray_exp_052,xray_flx_052,xray_err_052,xray_cnv_052]
sav_inds = [sav_inds,'IXBAND_052','IIDET_052','IDET_052', $
                                   iidet_052, idet_052]

;; SAVE all variables
sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="xray_soft052.sav"')


;; check flux-by-flux
if keyword_set(plt) then begin
    ;; plot detections
    iplt = where(iidet_nst_052 or iidet_xmm_052 or iidet_cha_052)
    
    e = {sym_size:0.5,sym_filled:1,color:'dodger blue', $
         xr:[1e-17,1e-11],yr:[1e-17,1e-11],xlog:1,ylog:1, $
         aspect_ratio:1,dimension:[1200,400], $
         buffer:1}
    ;; Chandra vs. XMM
    p = plot(flx_xmm_052[iplt],flx_cha_052[iplt],'S',_extra=e,layout=[3,1,1])
    p = plot(e.xr,e.yr,'--',/ov)
    p.xtitle = '$F_{XMM,0.5-2keV}(2-4.5 keV)$'
    p.ytitle = '$F_{Chandra,0.5-2keV}(2-7 keV)$'
    ;; NuSTAR vs. XMM
    p = plot(flx_xmm_052[iplt],flx_nst_052[iplt],'S',_extra=e,layout=[3,1,2],/current)
    p = plot(e.xr,e.yr,'--',/ov)
    p.xtitle = '$F_{XMM,0.5-2keV}(2-4.5 keV)$'
    p.ytitle = '$F_{NuSTAR,0.5-2keV}(3-8 keV)$'
    ;; NuSTAR vs. Chandra
    p = plot(flx_cha_052[iplt],flx_nst_052[iplt],'S',_extra=e,layout=[3,1,3],/current)
    p = plot(e.xr,e.yr,'--',/ov)
    p.xtitle = '$F_{Chanrda,0.5-2keV}(2-7 keV)$'
    p.ytitle = '$F_{NuSTAR,0.5-2keV}(3-8 keV)$'
    ;; title lists photon index
    t = text(0.5,0.9,'$\Gamma = $'+string(phot_ind,'(d3.1)'),alignment=0.5,/normal)
    ;; save!
    p.save,'compare_soft052_gamma_'+string(phot_ind,'(d3.1)')+'.png'
endif

END





