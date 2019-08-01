PRO convert_softx, phot_ind, $
                   PLT = plt


;; load variables
common _fits
common _inf_cha
common _inf_xmm
common _inf_nst
common _det_cha
common _det_xmm
common _det_nst

;; X-ray fields
xfield = '_'+['CHA','XMM','NST']
nfield = n_elements(xfield)

;; number of sources in X-ray fields
nsrc = n_elements(ra)
iiinf = strjoin('IIINF'+xfield,' or ')
re = execute('iiinf = '+iiinf)
if (n_elements(where(iiinf)) ne nsrc) then begin
    print, 'NUMBER OF IN-FIELD SOURCES DOES NOT MATCH'
    stop
endif

;; Instrument variables: exposure time, flux, error
tt = {CHA:strarr(6)+'ACIS_TIME', $
      XMM:strarr(7)+'PN_ONTIME', $
      NST:['S','H','F']+'EXP'}
ff = {CHA:'FLUX_POWLAW_APER90_'+['B','H','M','S','U','W'], $
      XMM:'PN_'+['1','2','3','4','5','8','9']+'_FLUX', $
      NST:['S','H','F']+'BF'}
ee = {CHA:'FLUX_POWLAW_APER90_'+['B','H','M','S','U','W']+'_ERR', $
      XMM:'PN_'+['1','2','3','4','5','8','9']+'_FLUX_ERR', $
      NST:'E_'+['S','H','F']+'BF'}
;; total number of energy bands
nxband = intarr(nfield)
for i = 0,nfield-1 do nxband[i] = n_elements(ff.(i))

sav_vars = ['XFIELD','NFIELD','PHOT_IND','NSRC','TT','FF','EE','NXBAND']
sav_inds = ['IIINF']


;; Flux conversions
;; WebPIMMS parameters: Galactic NH=2E20, power law photon index=1.8
;; https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3pimms/w3pimms.pl
case phot_ind of
    2.0: begin
        xconv = {CHA:[6.297E-01,1.285E+00,3.198E+00,2.011E+00,3.077E+00,4.648E-01], $
                 XMM:[3.077E+00,2.577E+00,2.370E+00,1.988E+00,1.638E+00,4.448E-01,7.616E-01], $
                 NST:[1.639E+00,1.462E+00,7.729E-01]}
        ;; Chandra                  ;; energy band (keV)
        ;out210_cha_b = 6.297E-01    ;; 0.5-7
        ;out210_cha_h = 1.285E+00    ;; 2-7
        ;out210_cha_m = 3.198E+00    ;; 1.2-2
        ;out210_cha_s = 2.011E+00    ;; 0.5-1.2
        ;out210_cha_u = 3.077E+00    ;; 0.2-0.5
        ;out210_cha_w = 4.648E-01    ;; 0.1-10
        ;; XMM
        ;out210_xmm_1 = 3.077E+00    ;; 0.2-0.5
        ;out210_xmm_2 = 2.577E+00    ;; 0.5-1
        ;out210_xmm_3 = 2.370E+00    ;; 1-2
        ;out210_xmm_4 = 1.988E+00    ;; 2-4.5
        ;out210_xmm_5 = 1.638E+00    ;; 4.5-12
        ;out210_xmm_8 = 4.448E-01    ;; 0.2-12
        ;out210_xmm_9 = 7.616E-01    ;; 0.5-4.5
        ;; NuSTAR
        ;out210_nst_s = 1.639E+00    ;; 3-8
        ;out210_nst_h = 1.462E+00    ;; 8-24
        ;out210_nst_f = 7.729E-01    ;; 3-24
        end
    1.8: begin
        xconv = {CHA:[7.400E-01,1.334E+00,3.969E+00,2.858E+00,5.147E+00,5.541E-01], $
                 XMM:[5.147E+00,3.734E+00,2.993E+00,2.160E+00,1.488E+00,5.178E-01,9.391E-01], $
                 NST:[1.614E+00,1.169E+00,6.781E-01]}
        ;; Chandra                  ;; energy band (keV)
        ;; b = 7.400E-01    ;; 0.5-7
        ;; h = 1.334E+00    ;; 2-7
        ;; m = 3.969E+00    ;; 1.2-2
        ;; s = 2.858E+00    ;; 0.5-1.2
        ;; u = 5.147E+00    ;; 0.2-0.5
        ;; w = 5.541E-01    ;; 0.1-10
        ;; XMM
        ;; 1 = 5.147E+00    ;; 0.2-0.5
        ;; 2 = 3.734E+00    ;; 0.5-1
        ;; 3 = 2.993E+00    ;; 1-2
        ;; 4 = 2.160E+00    ;; 2-4.5
        ;; 5 = 1.488E+00    ;; 4.5-12
        ;; 8 = 5.178E-01    ;; 0.2-12
        ;; 9 = 9.391E-01    ;; 0.5-4.5
        ;; NuSTAR
        ;; s = 1.614E+00    ;; 3-8
        ;; h = 1.169E+00    ;; 8-24
        ;; f = 6.781E-01    ;; 3-24
        end
    1.4: begin
        xconv = {CHA:[9.818E-01,1.452E+00,6.254E+00,5.880E+00,1.467E+01,7.141E-01], $
                 XMM:[1.467E+01,8.008E+00,4.875E+00,2.600E+00,1.246E+00,6.309E-01,1.399E+00], $
                 NST:[1.591E+00,7.576E-01,5.132E-01]}
        ;; Chandra                  ;; energy band (keV)
        ;; b = 9.818E-01    ;; 0.5-7
        ;; h = 1.452E+00    ;; 2-7
        ;; m = 6.254E+00    ;; 1.2-2
        ;; s = 5.880E+00    ;; 0.5-1.2
        ;; u = 1.467E+01    ;; 0.2-0.5
        ;; w = 7.141E-01    ;; 0.1-10
        ;; XMM
        ;; 1 = 1.467E+01    ;; 0.2-0.5
        ;; 2 = 8.008E+00    ;; 0.5-1
        ;; 3 = 4.875E+00    ;; 1-2
        ;; 4 = 2.600E+00    ;; 2-4.5
        ;; 5 = 1.246E+00    ;; 4.5-12
        ;; 8 = 6.309E-01    ;; 0.2-12
        ;; 9 = 1.399E+00    ;; 0.5-4.5
        ;; NuSTAR
        ;; s = 1.591E+00    ;; 3-8
        ;; h = 7.576E-01    ;; 8-24
        ;; f = 5.132E-01    ;; 3-24
        end
    else: print, 'PHOTON INDEX NOT IN RANGE'
endcase

;; flux and flux error for converted 2-10 keV
ff_210 = {CHA:ff.cha+'_210',XMM:ff.xmm+'_210',NST:ff.nst+'_210'}
ee_210 = {CHA:ee.cha+'_210',XMM:ee.xmm+'_210',NST:ee.nst+'_210'}

;; apply conversion factor and fill converted 2-10keV flux structure
for f = 0,nfield-1 do begin
    for b = 0,nxband[f]-1 do begin
        re = execute(ff_210.(f)[b]+' = '+ff.(f)[b]+' * xconv.(f)[b]')
        re = execute(ee_210.(f)[b]+' = '+ee.(f)[b]+' * xconv.(f)[b]')
    endfor
endfor

sav_vars = [sav_vars,'XCONV','FF_210','EE_210']
sav_inds = [sav_inds]


;; create output arrays for X-ray detections
xray_exp = 'EXP'+xfield
xray_flx = 'FLX'+xfield
xray_err = 'ERR'+xfield
xray_cnv = 'CNV'+xfield
for i = 0,nfield-1 do begin
    re = execute(xray_exp[i]+' = dblarr(nsrc)')
    re = execute(xray_flx[i]+' = dblarr(nsrc)')
    re = execute(xray_err[i]+' = dblarr(nsrc)')
    re = execute(xray_cnv[i]+' = dblarr(nsrc)')    
endfor

;; closest energy band to 2-10 keV from each instrument
;; Chandra==2-7keV, XMM==2-4.5keV, NuSTAR=3-8keV
ixband = [1,3,0]   
;; for each instrument sort by least fractional conversion factor from original chosen band (i.e., ixband_210)
for f = 0,nfield-1 do begin
    conv = xconv.(f)
    iorder = sort(abs((conv-conv[ixband[f]])/conv[ixband[f]]))
    for b = 0,nxband[f]-1 do begin
        re = execute('iuse = where('+tt.(f)[iorder[b]]+' gt 0. and '+ff_210.(f)[iorder[b]]+' gt 0. and '+xray_exp[f]+' eq 0.,nuse)')
        if (nuse gt 0.) then begin
            re = execute(xray_exp[f]+'[iuse] = '+tt.(f)[iorder[b]]+'[iuse]')
            re = execute(xray_flx[f]+'[iuse] = '+ff_210.(f)[iorder[b]]+'[iuse]')
            re = execute(xray_err[f]+'[iuse] = '+ee_210.(f)[iorder[b]]+'[iuse]')
            re = execute(xray_cnv[f]+'[iuse] = xconv.(f)[iorder[b]]')
        endif
    endfor
endfor

sav_vars = [sav_vars,xray_exp,xray_flx,xray_err,xray_cnv]
sav_inds = [sav_inds,'IXBAND']

;; SAVE all variables
sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="softx_conversions.sav"')


;; check flux-by-flux
if keyword_set(plt) then begin
    ;; plot detections
    iplt = where(iidet_nst or iidet_xmm or iidet_cha)
    
    e = {sym_size:0.5,sym_filled:1,color:'dodger blue', $
         xr:[1e-17,1e-11],yr:[1e-17,1e-11],xlog:1,ylog:1, $
         aspect_ratio:1,dimension:[1200,400], $
         buffer:1}
    ;; Chandra vs. XMM
    p = plot(flx_cha[iplt],flx_xmm[iplt],'S',_extra=e,layout=[3,1,1])
    p = plot(e.xr,e.yr,'--',/ov)
    p.xtitle = '$F_{Chandra,2-10keV}(2-7 keV)$'
    p.ytitle = '$F_{XMM,2-10keV}(2-4.5 keV)$'
    ;; NuSTAR vs. XMM
    p = plot(flx_xmm[iplt],flx_nst[iplt],'S',_extra=e,layout=[3,1,2],/current)
    p = plot(e.xr,e.yr,'--',/ov)
    p.xtitle = '$F_{XMM,2-10keV}(2-4.5 keV)$'
    p.ytitle = '$F_{NuSTAR,2-10keV}(3-8 keV)$'
    ;; NuSTAR vs. Chandra
    p = plot(flx_cha[iplt],flx_nst[iplt],'S',_extra=e,layout=[3,1,3],/current)
    p = plot(e.xr,e.yr,'--',/ov)
    p.xtitle = '$F_{Chanrda,2-10keV}(2-7 keV)$'
    p.ytitle = '$F_{NuSTAR,2-10keV}(3-8 keV)$'
    ;; title lists photon index
    t = text(0.5,0.9,'$\Gamma = $'+string(phot_ind,'(d3.1)'),alignment=0.5,/normal)
    ;; save!
    p.save,'plot_softx_gamma_'+string(phot_ind,'(d3.1)')+'.png'
endif


END





