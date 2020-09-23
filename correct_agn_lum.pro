FUNCTION correct_agn_lum, obswav, $
			              in_flx, $
			              in_err, $
                          in_fits


common _comp

;; determine which template components
components = tag_names(comp)
;; all possible templates (SED modeling procedure can handle max=5 templates)
temps = ['AGN','ELL','SFG','IRR','DST']   
;; match input components (use MATCH2.PRO to keep named order of TEMPS; MATCH.PRO alphabetizes; important for plotting purposes)
match2,components,temps,icomp,itemp
;; ensure we contain at least one valid template and sort
if (total(itemp ne -1) le 0) then stop		           
temps = temps[where(itemp ne -1)]
ntemps = n_elements(temps)
;; extract template normalizations
coeff = in_fits[2:2+ntemps-1,*]
;; number of input sources and output correction factor
nobj = n_elements(in_fits[0,*])
lcorr = dblarr(nobj)-9999.
;; extract indices of sources to correct
iagn = where(strmatch(temps,'AGN'),nagn)
if (nagn eq 0) then return, 0.
ind = where(coeff[iagn,*] gt 0.,nind)
if (nind eq 0) then return, 0.

;; output correction array
;lcorr = dblarr(n_elements(in_lum))
;; subset sources to correct
coeff = coeff[*,ind]
flx = in_flx[*,ind]
err = in_err[*,ind]
fits = in_fits[*,ind]
;; extract model parameters
ebv = fits[0,*]
z = fits[1,*]

;; calculate wavelength and frequency for sources and templates
objwav = rebin(obswav,n_elements(obswav),nind)
objnu = (!const.c*1e6)/objwav
tempwav = comp.wav#reform(1+z)
tempnu = (!const.c*1e6)/tempwav
;; target 6-micron wavelength in observed frame
targ_wav = 6.*(1+z)

;; covert data from flux density [microjansky] to flux [erg/s/cm2]
err *= 1e-29 * objnu         
flx *= 1e-29 * objnu

;; reconstruct models
;; convert models from flux density [microjansky] to flux [erg/s/cm2]
for i = 0,ntemps-1 do re = execute(temps[i]+' = 1e-29 * tempnu * (coeff[i,*]##comp.'+temps[i]+')')  ;; construct models
agn *= 10.^(-0.4 * comp.kap # ebv)                                                                  ;; add AGN extinction
re = execute('model = '+strjoin(temps,"+"))                                                         ;; coadded models

;; convert to log scale
err = abs((err)/(flx*alog(10)))
flx = alog10(flx)
for i = 0,ntemps-1 do re = execute(temps[i]+' = alog10('+temps[i]+')')
model = alog10(model)

;;========================================================================================
;;
;;              ADD AGN LUM CORRECTION HERE
;;
;;========================================================================================
temp_phot = dblarr(4,nind)
match,obswav,[3.4,4.6,12.,22.],iwise,i2
wav_wise = obswav[iwise]
flx_wise = flx[iwise,*]
err_wise = err[iwise,*]              ;; all non-finite values changed to zero
flx_perr = flx_wise+3.*err_wise
flx_merr = flx_wise-3.*err_wise

temp_flx = dblarr(n_elements(wav_wise),nind)
for i = 0,nind-1 do temp_flx[*,i] = interpol(model[*,i],tempwav[*,i],wav_wise)

;; where SED fit is above/below both W3+W4
iiabove = err_wise gt 0. and temp_flx gt flx_perr       ;; non-finite values accounted for here
iibelow = err_wise gt 0. and temp_flx lt flx_merr
iabove = where(total(iiabove[2:3,*],1) eq 2,nabove)
ibelow = where(total(iibelow[2:3,*],1) eq 2,nbelow)
if (nabove gt 0) then $
    for i = 0,nabove-1 do lcorr[ind[iabove[i]]] = interpol(flx_wise[*,iabove[i]],wav_wise,6.)-interpol(temp_flx[*,iabove[i]],wav_wise,targ_wav[iabove[i]])
if (nbelow gt 0) then $
    for i = 0,nbelow-1 do lcorr[ind[ibelow[i]]] = interpol(flx_wise[*,ibelow[i]],wav_wise,6.)-interpol(temp_flx[*,ibelow[i]],wav_wise,targ_wav[ibelow[i]])

icorr = where(lcorr ne -9999.,corrct)
if (corrct gt 0) then lcorr[icorr] = 10.^lcorr[icorr]
lcorr[where(lcorr eq -9999.,/null)] = 0.

return, lcorr


END












;temp_hi = temp_flx gt flx_perr
;temp_lo = temp_flx lt flx_merr
;;; 1st round, correct on W4
;iihi = total(temp_hi[-1,*],1) eq 1. and coeff[0,*] gt 0. and err_wise[-1,*] gt 0.
;iilo = total(temp_lo[-1,*],1) eq 1. and coeff[0,*] gt 0. and err_wise[-1,*] gt 0.
;ihi = where(iihi,cthi)
;ilo = where(iilo,ctlo)
;
;if (cthi gt 0) then begin
;    for i = 0,cthi-1 do begin
;        phot6 = interpol(flx_wise[*,ihi[i]],wav_wise,6.)
;        temp6 = interpol(temp_flx[*,ihi[i]],wav_wise,6.)
;        del6 = phot6-temp6
;        lum[ihi[i]] += del6
;    endfor
;endif
;
;if (ctlo gt 0) then begin
;    for i = 0,ctlo-1 do begin
;        phot6 = interpol(flx_wise[*,ilo[i]],wav_wise,6.)
;        temp6 = interpol(temp_flx[*,ilo[i]],wav_wise,6.)
;        del6 = phot6-temp6
;        lum[ilo[i]] -= del6
;    endfor
;endif
;
;;; 2nd round, correct on W3, where W4 doesn't exist
;iihi = total(temp_hi[-2,*],1) eq 1. and coeff[0,*] gt 0. and err_wise[-1,*] eq 0. and err_wise[-2,*] gt 0.
;iilo = total(temp_lo[-2,*],1) eq 1. and coeff[0,*] gt 0. and err_wise[-1,*] eq 0. and err_wise[-2,*] gt 0.
;ihi = where(iihi,cthi)
;ilo = where(iilo,ctlo)
;
;if (cthi gt 0) then begin
;    for i = 0,cthi-1 do begin
;        phot6 = interpol(flx_wise[*,ihi[i]],wav_wise,6.)
;        temp6 = interpol(temp_flx[*,ihi[i]],wav_wise,6.)
;        del6 = phot6-temp6
;        lum[ihi[i]] += del6
;    endfor
;endif
;
;if (ctlo gt 0) then begin
;    for i = 0,ctlo-1 do begin
;        phot6 = interpol(flx_wise[*,ilo[i]],wav_wise,6.)
;        temp6 = interpol(temp_flx[*,ilo[i]],wav_wise,6.)
;        del6 = phot6-temp6
;        lum[ilo[i]] -= del6
;    endfor
;endif


