FUNCTION correct_agn_lum, in_lum, $
                          obswav, $
			              in_flx, $
			              in_err, $
                          in_fits

common _comp
;; determine which template components
components = tag_names(comp)
;; all possible templates (SED modeling procedure can handle max=5 templates)
temps = ['AGN','ELL','SFG','IRR','DST']   
;; colors for plotting
col = [[204,121,167],[213,94,0],[0,158,115],[0,114,178],[240,228,66]]
;col = ['purple','purple','red','dark green','dark green','medium blue','brown']
;; match input components (use MATCH2.PRO to keep named order of TEMPS; MATCH.PRO alphabetizes; important for plotting purposes)
match2,components,temps,icomp,itemp
;; ensure we contain at least one valid template and sort
if (total(itemp ne -1) le 0) then stop		           
temps = temps[where(itemp ne -1)]
col = col[*,where(itemp ne -1)]
ntemps = n_elements(temps)

;; extract indices of sources to plot
if (n_elements(ind) eq 0) then ind = lindgen(n_elements(in_fits[0,*]))
nobj = n_elements(ind)
flx = in_flx[*,ind]
err = in_err[*,ind]
fits = in_fits[*,ind]
lum = in_lum[ind]

;; extract model parameters
ebv = fits[0,*]
z = fits[1,*]
coeff = fits[2:2+ntemps-1,*]
chi = fits[-2:-1,*]

;; calculate wavelength and frequency for sources and templates
if keyword_set(restframe) then begin
    objwav = obswav#(1+z)^(-1)
    objnu = (!const.c*1e6)/objwav#(1.+z)^(-1)
    tempwav = rebin(comp.wav,n_elements(comp),nobj)
    tempnu = (!const.c*1e6)/tempwav#(1.+z)^(-1)
    xtitle = '$Rest wavelength [ \mum ]$'
endif else begin
    objwav = rebin(obswav,n_elements(obswav),nobj)
    objnu = (!const.c*1e6)/objwav
    tempwav = comp.wav#reform(1+z)
    tempnu = (!const.c*1e6)/tempwav
    xtitle = '$Observed wavelength [ \mum ]$'
endelse

;; covert data from flux density [microjansky] to flux [erg/s/cm2]
err *= 1e-29 * objnu         
flx *= 1e-29 * objnu
;; reconstruct models
;; convert models from flux density [microjansky] to flux [erg/s/cm2]
agn = 1e-29 * tempnu * (coeff[0,*]##comp.(where(strmatch(tag_names(comp),'AGN*')))) * 10.^(-0.4 * comp.kap # ebv)                         ;; AGN model
for i = 1,ntemps-1 do re = execute(temps[i]+' = 1e-29 * tempnu * (coeff[i,*]##comp.'+temps[i]+')')  ;; galaxy models
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

temp_phot = dblarr(4,nobj)
match,obswav,[3.4,4.6,12.,22.],iwise,i2
wav_wise = obswav[iwise]
flx_wise = flx[iwise,*]
err_wise = err[iwise,*]>0.
flx_perr = flx_wise+err_wise*3.
flx_merr = flx_wise-err_wise*3.

temp_flx = dblarr(n_elements(wav_wise),nobj)
for i = 0,n_elements(ind)-1 do temp_flx[*,i] = interpol(model[*,i],tempwav[*,i],wav_wise)

temp_hi = temp_flx gt flx_perr
temp_lo = temp_flx lt flx_merr
;; 1st round, correct on W4
iihi = total(temp_hi[-1,*],1) eq 1. and coeff[0,*] gt 0. and err_wise[-1,*] gt 0.
iilo = total(temp_lo[-1,*],1) eq 1. and coeff[0,*] gt 0. and err_wise[-1,*] gt 0.
ihi = where(iihi,cthi)
ilo = where(iilo,ctlo)

if (cthi gt 0) then begin
    for i = 0,cthi-1 do begin
        phot6 = interpol(flx_wise[*,ihi[i]],wav_wise,6.)
        temp6 = interpol(temp_flx[*,ihi[i]],wav_wise,6.)
        del6 = phot6-temp6
        lum[ihi[i]] += del6
    endfor
endif

if (ctlo gt 0) then begin
    for i = 0,ctlo-1 do begin
        phot6 = interpol(flx_wise[*,ilo[i]],wav_wise,6.)
        temp6 = interpol(temp_flx[*,ilo[i]],wav_wise,6.)
        del6 = phot6-temp6
        lum[ilo[i]] -= del6
    endfor
endif

;; 2nd round, correct on W3, where W4 doesn't exist
iihi = total(temp_hi[-2,*],1) eq 1. and coeff[0,*] gt 0. and err_wise[-1,*] eq 0. and err_wise[-2,*] gt 0.
iilo = total(temp_lo[-2,*],1) eq 1. and coeff[0,*] gt 0. and err_wise[-1,*] eq 0. and err_wise[-2,*] gt 0.
ihi = where(iihi,cthi)
ilo = where(iilo,ctlo)

if (cthi gt 0) then begin
    for i = 0,cthi-1 do begin
        phot6 = interpol(flx_wise[*,ihi[i]],wav_wise,6.)
        temp6 = interpol(temp_flx[*,ihi[i]],wav_wise,6.)
        del6 = phot6-temp6
        lum[ihi[i]] += del6
    endfor
endif

if (ctlo gt 0) then begin
    for i = 0,ctlo-1 do begin
        phot6 = interpol(flx_wise[*,ilo[i]],wav_wise,6.)
        temp6 = interpol(temp_flx[*,ilo[i]],wav_wise,6.)
        del6 = phot6-temp6
        lum[ilo[i]] -= del6
    endfor
endif

return, lum


END



