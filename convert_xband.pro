PRO convert_xband, phot_ind


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

nsrc = n_elements(ra)
;; combine sources in X-ray fields
iiinf = strjoin('IIINF'+xfield,' or ')
re = execute('iiinf = '+iiinf)
if (n_elements(where(iiinf)) ne nsrc) then begin
    print, 'NUMBER OF IN-FIELD SOURCES DOES NOT MATCH'
    stop
endif
;; source matched to an X-ray catalog
iix = strjoin('IIX'+xfield,' or ')
re = execute('iix = '+iix)
;; source with valid X-ray detection
iidet = strjoin('IIDET'+xfield,' or ')
re = execute('iidet = '+iidet)

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
sav_inds = ['IIINF','IIX','IIDET']


;; Flux conversions
;; WebPIMMS parameters: Galactic NH=2E20, power law photon index=1.8
;; https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3pimms/w3pimms.pl
;; Chandra *1.1 to convert 90% aperture flux to 100%
;; energy band (keV)
;; Chandra          XMM (PN)         NuSTAR    
;; b  0.5-7         1  0.2-0.5       s  3-8
;; h  2-7           2  0.5-1         h  8-24
;; m  1.2-2         3  1-2           f  3-24
;; s  0.5-1.2       4  2-4.5
;; u  0.2-0.5       5  4.5-12
;; w  0.1-10        8  0.2-12
;;                  9  0.5-4.5
case phot_ind of
    2.0: begin
        print, 'UPDATE VALUES'
        xconv = {CHA:[6.297E-01,1.285E+00,3.198E+00,2.011E+00,3.077E+00,4.648E-01] * 1.1, $
                 XMM:[3.077E+00,2.577E+00,2.370E+00,1.988E+00,1.638E+00,4.448E-01,7.616E-01], $
                 NST:[1.640E+00,1.463E+00,7.732E-01]}
        end
    1.8: begin
        xconv = {CHA:[7.350E-01,1.334E+00,3.951E+00,2.794E+00,4.933E+00,5.488E-01] * 1.1, $
                 XMM:[4.933E+00,3.638E+00,2.976E+00,2.159E+00 ,1.488E+00,5.132E-01,9.310E-01], $
                 NST:[1.615E+00,1.170E+00,6.783E-01]}
        end
    1.4: begin
        print, 'UPDATE VALUES'
        xconv = {CHA:[9.818E-01,1.452E+00,6.254E+00,5.880E+00,1.467E+01,7.141E-01] * 1.1, $
                 XMM:[1.467E+01,8.008E+00,4.875E+00,2.600E+00,1.246E+00,6.309E-01,1.399E+00], $
                 NST:[1.591E+00,7.580E-01,5.134E-01]}
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

sav_vars = [sav_vars,'XCONV',ff_210.cha,ff_210.xmm,ff_210.nst,ee_210.cha,ee_210.xmm,ee_210.nst]
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
;; Chandra H == 2-7 keV, XMM PN4 == 2-4.5 keV, NuSTAR S == 3-8 keV
ixband = {cha:1+lonarr(nxband[0]), $
          xmm:3+lonarr(nxband[1]), $
          nst:0+lonarr(nxband[2])}
for i = 0,nfield-1 do begin
    ;; energy conversion factors per instrument
    conv = xconv.(i)
    ;; sort by least fractional conversion factor from original chosen band
    ixband.(i) = sort(abs(conv-conv[ixband.(i)[0]]))
    ;; source assigned flag
    iifill = bytarr(nsrc)
    for b = 0,nxband[i]-1 do begin
        ;; valid detection exists for specified energy band
        re = execute('iivalid = '+tt.(i)[ixband.(i)[b]]+' gt 0. and '+ff_210.(i)[ixband.(i)[b]]+' gt 0. and '+ee_210.(i)[ixband.(i)[b]]+' gt 0. and '+xray_exp[i]+' eq 0.')
        ;; detection has not already been assigned
        ifill = where(iivalid and iifill eq 0,nfill)
        if (nfill gt 0.) then begin
            ;; account for source assignment
            iifill[ifill] = 1
            re = execute(xray_exp[i]+'[ifill] = '+tt.(i)[ixband.(i)[b]]+'[ifill]')
            re = execute(xray_flx[i]+'[ifill] = '+ff_210.(i)[ixband.(i)[b]]+'[ifill]')
            re = execute(xray_err[i]+'[ifill] = '+ee_210.(i)[ixband.(i)[b]]+'[ifill]')
            re = execute(xray_cnv[i]+'[ifill] = conv[ixband.(i)[b]]')
        endif
    endfor
    nfill = total(iifill)
    re = execute('ndet = total(IIDET'+xfield[i]+')')
    if (nfill ne ndet) then begin
        print, 'ERROR OCCURED! CONVERTED SOURCES DO NOT MATCH DETECTIONS!'
        stop
    endif
endfor

sav_vars = [sav_vars,xray_exp,xray_flx,xray_err,xray_cnv]
sav_inds = [sav_inds,'IXBAND']

;; SAVE all variables
sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="xband_conversions.sav"')


END





