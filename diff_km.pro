FUNCTION diff_km, surv, $
                  time, $
                  binsz, $
                  ntot, $
                  BINC = binc


bstart = -2d;floor(min(time))
bend = 3d;ceil(max(time))
bins = [bstart:bend:binsz]
binc = bins+binsz/2.
nbins = n_elements(bins)
diff = dblarr(nbins)

mass = dblarr(n_elements(surv))
mass[0] = 1.-surv[0]
for i = 1,n_elements(time)-1 do mass[i] = abs(surv[i]-surv[i-1])

for j = 0,nbins-2 do begin
    bin = bins[j:j+1]
    ii = where(time gt bin[0] and time le bin[1],ninbin)
    if (ninbin eq 0) then continue
    ;; if you get here, there is a source in the bin
    for i = 0,ninbin-1 do diff[j] += mass[ii[i]]
    ;; check for remaining censored data
    isurv = where(time gt max(time[ii]),nsurv)
    if ((nsurv eq 0) or (total(mass[isurv]) eq 0)) then break
endfor
diff *= ntot

;; if there are still censored data, distribute their data in the remaining bins
if (nsurv gt 0) then begin
    ;; remaining mass
    rmass = ntot-total(diff)
    ;; upper bound of remaning time
    remt = max(time[isurv])
    ibins = [j:(where(remt lt bins))[0]-1]
    ;; take 1, evenly distribute remaining censored sources
    nbins = n_elements(ibins)
    diff[ibins] = rmass/nbins
    ;; take 2, distribute as the fraction of censored sources remaining
    ;for i = 0,n_elements(ibins)-1 do begin
    ;    bin = bins[ibins[i]:ibins[i]+1]
    ;    ii = where(time gt bin[0] and time le bin[1],ninbin)
    ;    diff[ibins[i]] = rmass*(1.*ninbin/nsurv)
    ;endfor
endif

return, diff


END



