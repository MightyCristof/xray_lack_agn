FUNCTION stack_nondet_flux, logl6um, $
                            dlumsq, $
                            binc, $
                            dc, $
                            det, $
                            ll


;inon = where(iifinal_non and iiinf_cha and iiwac,nobj)
;fx_dist_avg = dblarr(10000)
;fx_dist_med = dblarr(10000)
logfx_dist_avg = dblarr(10000)
logfx_dist_med = dblarr(10000)
for i = 0,9999 do begin
    ;; SINGLE DRAW
    ll_avg = mc_nondet_dist(logl6um,binc,dc,det,ll)
    loglx_non = ll_avg + logl6um
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
;    readcol,'wagn1_10x05.out',binc,dc,format='d,d',/silent,skipline=587
;    readcol,'wagn1.dat',det,ll,format='d,d'
;    ;inon = where(iifinal_non and iiwac,nobj)
;    ll_struct = {ll_dist:dblarr(nobj)}
;    ll_struct = replicate(ll_struct,10000)
;    ll_avg = dblarr(nobj)
;    ;for i = 0,9999 do ll_struct[i].ll_dist = mc_nondet_dist(lir[inon],binc,dc,det,ll)
;    for i = 0,9999 do ll_struct[i].ll_dist = mc_nondet_dist(lir,binc,dc,det,ll)
;
;    ;; bin size for histograms
;    ;bins = dblarr(10000)
;    ;for i = 0,10000-1 do bins[i] = scott(ll_struct[i].ll_dist)    
;    ;binsz = mean(bins)
;    binsz = 0.5
;    binmin = binc[0]-0.25
;    binmax = binc[-1]-0.25
;    hist_struct = {yhist:dblarr(10)}
;    hist_struct = replicate(hist_struct,10000)
;    ;; histogram each distribution
;    for i = 0,9999 do hist_struct[i].yhist = histogram(ll_struct[i].ll_dist,locations=tempx,bin=binsz,min=binmin,max=binmax)
;    xhist = tempx+0.25
;
;    ll_avg = mean(hist_struct.yhist,dim=2)
;    ll_med = median(hist_struct.yhist,dim=2)
;
;    loglx_non = ll_avg + loglir[inon]
;    lx_non = 10.^loglx_non
;    fx_non = lx_non / (4.*!const.pi*dl2[inon])
;endif


