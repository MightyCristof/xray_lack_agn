FUNCTION rl_surv, in_xis, $
                  in_cis, $
                  bin


;; KM estimator
;; S = ¹(ti < t) (1 - dj/yj)

;; Survival function
;; ti = distinct event (event) times
;; dj = # of event events at time ti
;; yj = # of survivals until time T (end), at least as big as ti

;; In astronomy terms
;; ti = RL of a detection
;; dj = # of detections at ti
;; yj = # of sources yet to be accounted for, both detections & non-detections

;; sort in ascending order of "time" AKA luminosity ratio
isort = sort(in_xis)
xis = in_xis[isort]
cis = in_cis[isort]
xis[where(cis eq -1)] *= -1
plothist, xis
minxi = floor(xis[0])
maxxi = ceil(xis[-1])
;; initial survival func
s = 1.
v = 0.
x = minxi

;; continuous variable
if (n_elements(bin) eq 0) then begin
    ievent = where(cis eq 0,nevents)

    xi = xis[ievent[0]]

    for i = 0,n_elements(xis)-1 do begin
        ;; iterate until an event (cis == 0)
        if (cis[i] eq -1) then continue
        ;; event time
        xi = xis[i]
        ;; survivals 
        xis_surv = xis[i:-1]
        cis_surv = cis[i:-1]    
        ;; find # events at xi
        idi = where(xis_surv eq xi and cis_surv eq 0,di)
        ;; find # of survivals past ti
        ni = 1.*n_elements(xis_surv)
        ;; cannot have more current than remaining sources (for variance calc.)
        if (di eq ni) then continue
        ;; survival function at xi
        si = (1. - di/ni)
        ;; product of
        x = [x,xi]
        v = [v,s[-1]^2*(di/ni^2)*1./si]
        s = [s,s[-1]*si]
    endfor
endif else $
;; discrete variable
if (n_elements(bin) gt 0) then begin
    xiint = [minxi:maxxi:bin]

    for i = 0,n_elements(xiint)-2 do begin
        ;; values in interval
        iin = where(xis ge xiint[i] and xis lt xiint[i+1],nin)
        if (nin eq 0) then continue
        ievent = where(cis[iin] eq 0,nevent)
        ;; find bin with detections
        if (nevent eq 0) then continue
        isurv = where(xis ge xiint[i+1],nsurv)
        ;; cannot have more current than remaining sources
        if (nin ge nsurv) then continue
        ;; ALL objects in bin, not just events
        di = 1.*nin
        ;; sources that survive, not including current invterval
        ni = 1.*nsurv
        ;print, di, ni, t[-1], s[-1]
        si = (1. - di/ni)
        x = [x,xiint[i]+bin/2.]
        v = [v,s[-1]^2*(di/ni^2)*1./si]
        s = [s,s[-1]*si]
    endfor
endif
sfunc = {x:x,s:s,v:v}

return, soa2aos(sfunc)


END



