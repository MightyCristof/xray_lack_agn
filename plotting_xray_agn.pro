



;; all X-ray detections
iidet_all = (iidet_nst_210 or iidet_xmm_210 or iidet_cha_210) and iiebv and iichi
idet_all = where(iidet_all)

ix_nst = where(iidet_nst_210 and iiebv and iichi)
ix_xmm = where(iidet_xmm_210 and iiebv and iichi)
ix_cha = where(iidet_cha_210 and iiebv and iichi)
;; corrected fluxes
iitype1 = iilir and iiebv and iichi and ebv lt 0.2 and agnf.obs gt 0.9
itype1 = where(iitype1)

wave = wave_inf
flux = flux_inf[*,itype1]
e_flux = e_flux_inf[*,itype1]
bin = bin_inf[*,itype1]
objid = objid_inf[itype1]
param = param_inf[*,itype1]
z = z_inf[itype1]
ebv_inf = ebv
ebv = ebv[itype1]
lum = lir
fagn = agnf[itype1].obs

;; calculate template flux at WISE 
temp_phot = interp_ir_temp(wave,flux,ebv,z,param[2:5,*])

;; is the template constrained by the source WISE photometry - overestimate
data_upper = flux[5:8,*]+e_flux[5:8,*]
delta_over = temp_phot gt data_upper
iiover = (total(data_upper[2:3,*] gt 0. and delta_over[2:3,*],1) eq 2) OR $         ;; W3+W4 exist and template overshoots
         (reform(data_upper[3,*]) eq 0. and total(delta_over[1:2,*],1) eq 2) OR $   ;; W4 missing and template overshoots W2+W3
         (reform(data_upper[2,*]) eq 0. and total(delta_over[[1,3],*],1) eq 2)      ;; W3 missing and template overshoots W2+W4 (found none)
iover = where(iiover)

flux_over = interp_6micr_agn(wave,flux[*,iover],ebv[iover],z[iover],param[2,iover],fagn[iover],/dered)
lum_over = interp_6micr_lum(flux_over,z[iover],/log)

;; is the template constrained by the source WISE photometry - underestimate
data_lower = flux[5:8,*]-e_flux[5:8,*]
delta_under = temp_phot lt data_lower
iiunder = (total(data_lower[2:3,*] gt 0. and delta_under[2:3,*],1) eq 2)
iunder = where(iiunder)

flux_under = interp_6micr_agn(wave,flux[*,iunder],ebv[iunder],z[iunder],param[2,iunder],fagn[iunder],/dered)
lum_under = interp_6micr_lum(flux_under,z[iunder],/log)

lum[itype1[iover]] = lum_over
lum[itype1[iunder]] = lum_under

;; LX-LIR relation line
rel_lir = dindgen(60,start=20)/2.+20
rel_lx = dblarr(n_elements(rel_lir))
ilo = where(rel_lir lt 44.79)
ihi = where(rel_lir ge 44.79)
rel_lx[ilo] = 0.84*(rel_lir[ilo]-45.)+44.60
rel_lx[ihi] = 0.40*(rel_lir[ihi]-45.)+44.51

e = {symbol:'o',linestyle:'',sym_filled:1,sym_size:0.5,transparency:75, $
     xra:[42.,48.],yra:[42.,48.],aspect_ratio:1,margin:[0.2,0.2,0.2,0.2]}

p = plot(lir[ix_nst],lx_nst_210[ix_nst],col='purple',title='NuSTAR',_extra=e,layout=[3,2,1],dimension=[900,450],ytitle='$L_X$')
p = plot(rel_lir,rel_lx,/ov)
p = plot(lir[ix_xmm],lx_xmm_210[ix_xmm],col='orange',title='XMM',_extra=e,layout=[3,2,2],/current)
p = plot(rel_lir,rel_lx,/ov)
p = plot(lir[ix_cha],lx_cha_210[ix_cha],col='dodger blue',title='Chandra',_extra=e,layout=[3,2,3],/current)
p = plot(rel_lir,rel_lx,/ov)

p = plot(lir[iagn_det_nst_210],lx_nst_210[iagn_det_nst_210],col='purple',_extra=e,layout=[3,2,4],/current)
p = plot(rel_lir,rel_lx,/ov)
p = plot(lir[iagn_det_xmm_210],lx_xmm_210[iagn_det_xmm_210],col='orange',_extra=e,layout=[3,2,5],/current,xtitle='$L_{IR}$')
p = plot(rel_lir,rel_lx,/ov)
p = plot(lir[iagn_det_cha_210],lx_cha_210[iagn_det_cha_210],col='dodger blue',_extra=e,layout=[3,2,6],/current)
p = plot(rel_lir,rel_lx,/ov)



p.save,'lir_editted.png'






;; plot for NERQUAM 2019
lumx = lx_cha_210
lumx[where(lumx eq 0.)] = lx_xmm_210[where(lumx eq 0.)]
lumx[where(lumx eq 0.)] = lx_nst_210[where(lumx eq 0.)]

rel_lir = dindgen(60,start=20)/2.+20
rel_lx = dblarr(n_elements(rel_lir))
ilo = where(rel_lir lt 44.79)
ihi = where(rel_lir ge 44.79)
rel_lx[ilo] = 0.84*(rel_lir[ilo]-45.)+44.60
rel_lx[ihi] = 0.40*(rel_lir[ihi]-45.)+44.51

e = {symbol:'o',linestyle:'',sym_filled:1,sym_size:0.5,transparency:75, $
     xra:[41.,49.],yra:[40.,48.],aspect_ratio:1}

p = plot(lir[idet_all],lumx[idet_all],col='dodger blue',xtitle='$log  \itL\rm_{IR}$',ytitle='$log  \itL\rm_X$',_extra=e)
p = plot(rel_lir,rel_lx,'--',/ov)

;; WAC
p = plot(lir[ixdet],lumx[ixdet],col='dodger blue',xtitle='$log L_{IR}$',ytitle='$log L_X$',_extra=e)
p = plot(rel_lir,rel_lx,/ov)




p = plot(lir[ix_nst],lx_nst_210[ix_nst],col='purple',title='NuSTAR',_extra=e,layout=[3,1,1],dimension=[900,300],ytitle='$log L_X$')
p = plot(rel_lir,rel_lx,/ov)
p = plot(lir[ix_xmm],lx_xmm_210[ix_xmm],col='orange',title='XMM',_extra=e,layout=[3,1,2],/current,xtitle='$log L_{IR}$')
p = plot(rel_lir,rel_lx,/ov)
p = plot(lir[ix_cha],lx_cha_210[ix_cha],col='dodger blue',title='Chandra',_extra=e,layout=[3,1,3],/current)
p = plot(rel_lir,rel_lx,/ov)







