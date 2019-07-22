PRO lum_ratio 


common _fits       
common _inf_cha    
common _inf_xmm    
common _inf_nst    
common _det_cha    
common _det_xmm    
common _det_nst 
common _det_wac    
common _softx      
common _fluxlim    
common _comp       
common _agn_lum    
common _clean_cha
common _clean_xmm
common _clean_nst
common _quality













END





;;----------------------------------------------------------------------------------------
;; SAMPLE AGN
;;----------------------------------------------------------------------------------------


iiagn_det = 'IIAGN_DET'+xfield             ;; Candidate AGN w.  X-ray detections
iiagn_lim = 'IIAGN_LIM'+xfield             ;; Candidate AGN w.o X-ray detections
iagn_det = 'IAGN_DET'+xfield
iagn_lim = 'IAGN_LIM'+xfield

iiagn_det_210 = iiagn_det+'_210'
iiagn_lim_210 = iiagn_lim+'_210'
iagn_det_210 = iagn_det+'_210'
iagn_lim_210 = iagn_lim+'_210'

for i = 0,nfield-1 do begin
    re = execute(iiagn_det_210[i]+' = iiagn and '+iidet_210[i])
    re = execute(iiagn_lim_210[i]+' = iiagn and '+iiinf[i])
    re = execute(iagn_det_210[i]+' = where('+iiagn_det_210[i]+')')
    re = execute(iagn_lim_210[i]+' = where('+iiagn_lim_210[i]+')')
endfor

;iiagn_det_052 = iiagn_det+'_052'
;iiagn_lim_052 = iiagn_lim+'_052'
;iagn_det_052 = iagn_det+'_052'
;iagn_lim_052 = iagn_lim+'_052'
;
;for i = 0,nfield-1 do begin
;    re = execute(iiagn_det_052[i]+' = iiagn and '+iidet_052[i])
;    re = execute(iiagn_lim_052[i]+' = iiagn and '+iiinf[i])
;    re = execute(iagn_det_052[i]+' = where('+iiagn_det_052[i]+')')
;    re = execute(iagn_lim_052[i]+' = where('+iiagn_lim_052[i]+')')
;endfor

sav_vars = [sav_vars,'AGNF15']
sav_inds = [sav_inds,'IIIRB','IIIRW','IIEBV','IIAGN', $
                     'IIAGN_DET','IIAGN_LIM','IAGN_DET','IAGN_LIM', $
                     'IIAGN_DET_210','IIAGN_LIM_210','IAGN_DET_210','IAGN_LIM_210', $
                     iiagn_det_210,iiagn_lim_210,iagn_det_210,iagn_lim_210];, $
                     ;'IIAGN_DET_052','IIAGN_LIM_052','IAGN_DET_052','IAGN_LIM_052', $
                     ;iiagn_det_052,iiagn_lim_052,iagn_det_052,iagn_lim_052]

;;----------------------------------------------------------------------------------------
;; LUMINOSITY RATIOS -- PROXY FOR OBSCURATION
;;----------------------------------------------------------------------------------------
;; detections/limits luminosity ratios
lldet_210 = 'LLDET'+xfield+'_210'
lllim_210 = 'LLLIM'+xfield+'_210'
for i = 0,nfield-1 do begin
    re = execute(lldet_210[i]+' = dblarr(nsrc)-9999.')
    re = execute(lllim_210[i]+' = dblarr(nsrc)-9999.')
    re = execute(lldet_210[i]+'['+iagn_det_210[i]+'] = '+lx_210[i]+'['+iagn_det_210[i]+']-lxir_210['+iagn_det_210[i]+']')
    re = execute(lllim_210[i]+'['+iagn_lim_210[i]+'] = '+lxlim_210[i]+'['+iagn_lim_210[i]+']-lxir_210['+iagn_lim_210[i]+']')
endfor

;lldet_052 = 'LLDET'+xfield+'_052'
;lllim_052 = 'LLLIM'+xfield+'_052'
;for i = 0,nfield-1 do begin
;    re = execute(lldet_052[i]+' = dblarr(nsrc)-9999.')
;    re = execute(lllim_052[i]+' = dblarr(nsrc)-9999.')
;    re = execute(lldet_052[i]+'['+iagn_det_052[i]+'] = '+lx_052[i]+'['+iagn_det_052[i]+']-lxir_052['+iagn_det_052[i]+']')
;    re = execute(lllim_052[i]+'['+iagn_lim_052[i]+'] = '+lxlim_052[i]+'['+iagn_lim_052[i]+']-lxir_052['+iagn_lim_052[i]+']')
;endfor

;; 10-40 keV
;lldet_1040 = 'LLDET'+xfield+'_1040'
;lllim_1040 = 'LLLIM'+xfield+'_1040'
;for i = 0,nfield-1 do begin
;    re = execute(lldet_1040[i]+' = dblarr(nsrc)-9999.')
;    re = execute(lllim_1040[i]+' = dblarr(nsrc)-9999.')
;    re = execute(lldet_1040[i]+'['+iagn_det_210[i]+'] = ('+lx_210[i]+'['+iagn_det_210[i]+']+alog10(1.163E+00))-(lxir_210['+iagn_det_210[i]+']+alog10(1.163E+00))')
;    re = execute(lllim_1040[i]+'['+iagn_lim_210[i]+'] = '+lxlim_210[i]+'['+iagn_lim_210[i]+']-lxir_210['+iagn_lim_210[i]+']')
;endfor

sav_vars = [sav_vars,'LLDET_210','LLLIM_210',lldet_210,lllim_210]
                     ;'LLDET_052','LLLIM_052',lldet_052,lllim_052, $
                     ;'LLDET_1040','LLLIM_1040',lldet_1040,lllim_1040, $
sav_inds = [sav_inds]


;;----------------------------------------------------------------------------------------
;; Luminosity data subset by field
;;----------------------------------------------------------------------------------------
if keyword_set(combine) then begin

    ;; data value: XMM over Chandra over NuSTAR
    lld_210 = lldet_cha_210
    lld_210[where(lld_210 eq -9999.)] = lldet_xmm_210[where(lld_210 eq -9999.)]
    lld_210[where(lld_210 eq -9999.)] = lldet_nst_210[where(lld_210 eq -9999.)]
    lll_210 = lllim_cha_210
    lll_210[where(lll_210 eq -9999.)] = lllim_xmm_210[where(lll_210 eq -9999.)]
    lll_210[where(lll_210 eq -9999.)] = lllim_nst_210[where(lll_210 eq -9999.)]
    ;; convert to 10-40 keV
    lld_1040 = lld_210
    lld_1040[where(lld_210 gt -99.,/null)] += alog10(1.163E+00)
    lll_1040 = lll_210
    lll_1040[where(lll_210 gt -99.,/null)] += alog10(1.163E+00)
    ;; combined detections/limits indices
    ixdet = where(lld_210 gt -9999.)
    ixlim = where(lll_210 gt -9999.)
    ;; AGN Catalog matches
    iwdet = where(lld_210 gt -9999. and iiwac)
    iwlim = where(lll_210 gt -9999. and iiwac)
    ;; log E(B-V)    
    lebv = alog10(ebv)>(-2.5)
    lebv = lebv + randomu(seed,n_elements(ebv))*0.05-0.05/2.

    sav_vars = [sav_vars,'LLD_210','LLL_210', $
                         ;'lld_1040','LLL_1040', $
                         'LEBV']
    sav_inds = [sav_inds,'IXDET','IXLIM','IWDET','IWLIM']
 endif
