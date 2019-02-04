PRO xray_lacking_agn_soft210, phot_ind, $
                              PLT = plt


;; load variables
;common _inf_fits
common _inf_nst
common _inf_xmm
common _inf_cha
common _det_nst
common _det_xmm
common _det_cha

;; X-ray fields
xfield = '_'+['NST','XMM','CHA']
nfield = n_elements(xfield)
texp = 'TEXP'+xfield
sdst = 'SDST'+xfield

;; number of sources in X-ray fields
;nsrc = n_elements(iiinf_nst)
xinf = 'IINF'+xfield
xdet = 'IDET'+xfield
xxinf = 'I'+xinf
xxdet = 'I'+xdet
;for i = 0,nfield-1 do begin
;    re = execute(xinf[i]+' = where(IIINF'+xfield[i]+',NINF'+xfield[i]+')')  ;; in field
;    re = execute(xdet[i]+' = where(IIDET'+xfield[i]+',NDET'+xfield[i]+')')  ;; detections
;endfor
iidet = iidet_nst or iidet_xmm or iidet_cha
idet = where(iidet,ndet)

sav_vars = ['XFIELD','NFIELD','TEXP','SDST','PHOT_IND']
sav_inds = ['XXINF','XXDET','XINF','XDET','IIDET','IDET']

;; Instrument variables: exposure time, flux, error
tt = [['S','H','F']+'EXP', $
      strarr(7)+'EP_ONTIME', $
      strarr(6)+'ACIS_TIME']
ff = [['S','H','F']+'BF', $
      'SC_EP_'+['1','2','3','4','5','8','9']+'_FLUX', $
      'FLUX_POWLAW_APER90_'+['B','H','M','S','U','W']]
ee = ['E_'+['S','H','F']+'BF', $
      'SC_EP_'+['1','2','3','4','5','8','9']+'_FLUX_ERR', $
      'FLUX_POWLAW_APER90_LOLIM_'+['B','H','M','S','U','W']]
ff_210 = ff+'_210'
ee_210 = ee+'_210'
nxband = n_elements(ff)

;; Flux conversions
;; WebPIMMS parameters: Galactic NH=2E20, power law photon index=1.8
;; https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3pimms/w3pimms.pl
case phot_ind of
    2.0: begin
        ;; NuSTAR                   ;; energy band (keV)
        out210_nst_s = 1.639E+00    ;; 3-8
        out210_nst_h = 1.462E+00    ;; 8-24
        out210_nst_f = 7.729E-01    ;; 3-24
        ;; XMM
        out210_xmm_1 = 3.077E+00    ;; 0.2-0.5
        out210_xmm_2 = 2.577E+00    ;; 0.5-1
        out210_xmm_3 = 2.370E+00    ;; 1-2
        out210_xmm_4 = 1.988E+00    ;; 2-4.5
        out210_xmm_5 = 1.638E+00    ;; 4.5-12
        out210_xmm_8 = 4.448E-01    ;; 0.2-12
        out210_xmm_9 = 7.616E-01    ;; 0.5-4.5
        ;; Chandra
        out210_cha_b = 6.297E-01    ;; 0.5-7
        out210_cha_h = 1.285E+00    ;; 2-7
        out210_cha_m = 3.198E+00    ;; 1.2-2
        out210_cha_s = 2.011E+00    ;; 0.5-1.2
        out210_cha_u = 3.077E+00    ;; 0.2-0.5
        out210_cha_w = 4.648E-01    ;; 0.1-10
        end
    1.8: begin
        ;; NuSTAR                   ;; energy band (keV)
        out210_nst_s = 1.614E+00    ;; 3-8
        out210_nst_h = 1.169E+00    ;; 8-24
        out210_nst_f = 6.781E-01    ;; 3-24
        ;; XMM
        out210_xmm_1 = 5.147E+00    ;; 0.2-0.5
        out210_xmm_2 = 3.734E+00    ;; 0.5-1
        out210_xmm_3 = 2.993E+00    ;; 1-2
        out210_xmm_4 = 2.160E+00    ;; 2-4.5
        out210_xmm_5 = 1.488E+00    ;; 4.5-12
        out210_xmm_8 = 5.178E-01    ;; 0.2-12
        out210_xmm_9 = 9.391E-01    ;; 0.5-4.5
        ;; Chandra
        out210_cha_b = 7.400E-01    ;; 0.5-7
        out210_cha_h = 1.334E+00    ;; 2-7
        out210_cha_m = 3.969E+00    ;; 1.2-2
        out210_cha_s = 2.858E+00    ;; 0.5-1.2
        out210_cha_u = 5.147E+00    ;; 0.2-0.5
        out210_cha_w = 5.541E-01    ;; 0.1-10
        end
    1.4: begin
        ;; NuSTAR                   ;; energy band (keV)
        out210_nst_s = 1.591E+00    ;; 3-8
        out210_nst_h = 7.576E-01    ;; 8-24
        out210_nst_f = 5.132E-01    ;; 3-24
        ;; XMM
        out210_xmm_1 = 1.467E+01    ;; 0.2-0.5
        out210_xmm_2 = 8.008E+00    ;; 0.5-1
        out210_xmm_3 = 4.875E+00    ;; 1-2
        out210_xmm_4 = 2.600E+00    ;; 2-4.5
        out210_xmm_5 = 1.246E+00    ;; 4.5-12
        out210_xmm_8 = 6.309E-01    ;; 0.2-12
        out210_xmm_9 = 1.399E+00    ;; 0.5-4.5
        ;; Chandra
        out210_cha_b = 9.818E-01    ;; 0.5-7
        out210_cha_h = 1.452E+00    ;; 2-7
        out210_cha_m = 6.254E+00    ;; 1.2-2
        out210_cha_s = 5.880E+00    ;; 0.5-1.2
        out210_cha_u = 1.467E+01    ;; 0.2-0.5
        out210_cha_w = 7.141E-01    ;; 0.1-10
        end
    else: print, 'PHOTON INDEX NOT IN RANGE'
endcase
;; energy conversions
xconv = ['OUT210_NST_S','OUT210_NST_H','OUT210_NST_F', $
         'OUT210_XMM_1','OUT210_XMM_2','OUT210_XMM_3','OUT210_XMM_4','OUT210_XMM_5','OUT210_XMM_8','OUT210_XMM_9', $
         'OUT210_CHA_B','OUT210_CHA_H','OUT210_CHA_M','OUT210_CHA_S','OUT210_CHA_U','OUT210_CHA_W']

;; structure containing flux conversions for all instruments/energy bands
gmma = 'GAMMA:'+string(phot_ind,"(d3.1)")
ff_str = '"'+ff+'",0d'
re = execute('cnv_210 = create_struct(name=gmma,'+strjoin(ff_str,",")+')')
;; fill conversion structure by photon index
;; apply conversion to fluxes and errors
for i = 0,nxband-1 do begin
    re = execute('cnv_210.(i) = '+xconv[i])
    re = execute(ff_210[i]+' = '+ff[i]+' * cnv_210.(i)')
    re = execute(ee_210[i]+' = '+ee[i]+' * cnv_210.(i)')
endfor

sav_vars = [sav_vars,'TT','FF','EE', $
                     'FF_210','EE_210',ff_210,ee_210, $
                     'NXBAND','CNV_210']
sav_inds = [sav_inds]


;; commit fluxes and errors to source for use
;; currently choosing one band per telescope
xray_exp_210 = 'EXP'+xfield+'_210'
xray_flx_210 = 'FLX'+xfield+'_210'
xray_err_210 = 'ERR'+xfield+'_210'
xray_cnv_210 = 'CNV'+xfield+'_210'

ixband_210 = [0,6,11]
used_exp = tt[ixband_210]
used_flx = ff_210[ixband_210]
used_err = ee_210[ixband_210]
used_cnv = 'CNV_210.'+(tag_names(cnv_210))[ixband_210]
for i = 0,nfield-1 do begin
    re = execute(xray_exp_210[i]+' = '+used_exp[i])
    re = execute(xray_flx_210[i]+' = '+used_flx[i])
    re = execute(xray_err_210[i]+' = '+used_err[i])
    re = execute(xray_cnv_210[i]+' = '+used_cnv[i])
endfor

sav_vars = [sav_vars,'XRAY_EXP_210','XRAY_FLX_210','XRAY_ERR_210','XRAY_CNV_210', $
                     xray_exp_210,xray_flx_210,xray_err_210,xray_cnv_210]
sav_inds = [sav_inds,'IXBAND_210']

;; SAVE all variables
sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="xray_soft210.sav"')


;; check flux-by-flux
if keyword_set(plt) then begin
    e = {sym_size:0.5,sym_filled:1,color:'dodger blue', $
         xr:[1e-17,1e-11],yr:[1e-17,1e-11],xlog:1,ylog:1, $
         aspect_ratio:1,dimension:[1200,400], $
         buffer:0}
    ;; Chandra vs. XMM
    p = plot(flx_xmm_210[idet],flx_cha_210[idet],'S',_extra=e,layout=[3,1,1])
    p = plot(e.xr,e.yr,'--',/ov)
    p.xtitle = '$F_{XMM,2-10keV}(2-4.5 keV)$'
    p.ytitle = '$F_{Chandra,2-10keV}(2-7 keV)$'
    ;; NuSTAR vs. XMM
    p = plot(flx_xmm_210[idet],flx_nst_210[idet],'S',_extra=e,layout=[3,1,2],/current)
    p = plot(e.xr,e.yr,'--',/ov)
    p.xtitle = '$F_{XMM,2-10keV}(2-4.5 keV)$'
    p.ytitle = '$F_{NuSTAR,2-10keV}(3-8 keV)$'
    ;; NuSTAR vs. Chandra
    p = plot(flx_cha_210[idet],flx_nst_210[idet],'S',_extra=e,layout=[3,1,3],/current)
    p = plot(e.xr,e.yr,'--',/ov)
    p.xtitle = '$F_{Chanrda,2-10keV}(2-7 keV)$'
    p.ytitle = '$F_{NuSTAR,2-10keV}(3-8 keV)$'
    ;; title lists photon index
    t = text(0.5,0.9,'$\Gamma = $'+string(phot_ind,'(d3.1)'),alignment=0.5,/normal)
    ;; save!
    p.save,'compare_soft210_gamma_'+string(phot_ind,'(d3.1)')+'.png'
    ;; gamma = 1.4, 
    ;; exposure time vs exposure time
    ;; flux limit vs flux limit
    ;; hickox+06 (+07 ??) photon index function of x-ray flux
endif


END





