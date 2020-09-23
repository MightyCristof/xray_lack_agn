FUNCTION stack_nondet_flux, logl6um, $
                            dlumsq, $
                            wbinc, $
                            wdc, $
                            det, $
                            rl


;inon = where(iifinal_non and iiinf_cha and iidet_wac,nobj)
;fx_dist_avg = dblarr(10000)
;fx_dist_med = dblarr(10000)
logfx_dist_avg = dblarr(10000)
logfx_dist_med = dblarr(10000)
for i = 0,9999 do begin
    ;; SINGLE DRAW
    rl_avg = mc_nondet_dist(logl6um,wbinc,wdc,det,rl)
    loglx_non = rl_avg + logl6um
    ;; LOG SPACE
    logfx_non = loglx_non - alog10(4.*!const.pi*dlumsq)
    logfx_dist_avg[i] = mean(logfx_non)
    logfx_dist_med[i] = median(logfx_non)
    ;; LINEAR SPACE
    ;lx_non = 10.^loglx_non
    ;fx_non = lx_non / (4.*!const.pi*dlumsq)
    ;fx_dist_avg[i] = mean(fx_non)
    ;fx_dist_med[i] = median(fx_non)
endfor

;yfxavg = histogram(fx_dist_avg,locations=xfxavg,bin=scott(fx_dist_avg))
yfxavg = histogram(logfx_dist_avg,locations=xfxavg,bin=scott(logfx_dist_avg))
yfxmed = histogram(logfx_dist_med,locations=xfxmed,bin=scott(logfx_dist_med))

fx_stack = {fxavg:mean(logfx_dist_avg),yavg:yfxavg,xavg:xfxavg, $
            fxmed:mean(logfx_dist_med),ymed:yfxmed,xmed:xfxmed}

return, fx_stack


END





;if keyword_set(oldstuff) then begin
;    readcol,'wagn1_10x05.out',wbinc,wdc,format='d,d',/silent,skipline=587
;    readcol,'wagn1.dat',det,rl,format='d,d'
;    ;inon = where(iifinal_non and iidet_wac,nobj)
;    rl_struct = {rl_dist:dblarr(nobj)}
;    rl_struct = replicate(rl_struct,10000)
;    rl_avg = dblarr(nobj)
;    ;for i = 0,9999 do rl_struct[i].rl_dist = mc_nondet_dist(lir[inon],wbinc,wdc,det,rl)
;    for i = 0,9999 do rl_struct[i].rl_dist = mc_nondet_dist(lir,wbinc,wdc,det,rl)
;
;    ;; bin size for histograms
;    ;bins = dblarr(10000)
;    ;for i = 0,10000-1 do bins[i] = scott(rl_struct[i].rl_dist)    
;    ;binsz = mean(bins)
;    binsz = 0.5
;    binmin = wbinc[0]-0.25
;    binmax = wbinc[-1]-0.25
;    hist_struct = {yhist:dblarr(10)}
;    hist_struct = replicate(hist_struct,10000)
;    ;; histogram each distribution
;    for i = 0,9999 do hist_struct[i].yhist = histogram(rl_struct[i].rl_dist,locations=tempx,bin=binsz,min=binmin,max=binmax)
;    xhist = tempx+0.25
;
;    rl_avg = mean(hist_struct.yhist,dim=2)
;    rl_med = median(hist_struct.yhist,dim=2)
;
;    loglx_non = rl_avg + loglir[inon]
;    lx_non = 10.^loglx_non
;    fx_non = lx_non / (4.*!const.pi*dl2[inon])
;endif


