;; resample the RL 2 NH distribution to find the error on NH due to the LX-LIR relation
FUNCTION resamp_rlnh_err, ind, $
                          RET = ret, $
                          MODEL = model


;STOP
;; this needs to get updated because we're not using LSCAT in the RL calculations anymore
;; NO IT FUCKING DOESN'T, CHECK LINE#32 ASSHOLE
common _wac
common _xconv
common _agnlum
common _quality
common _combined

;; additional indexing aside from X-ray detected
if (n_elements(ind) eq 0) then ind = lonarr(nsrc)+1

;; RL2NH model must be assigned
if ~keyword_set(model) then begin
    print, '====================================='
    print, 'NO MODEL SET. RUNNING POWER LAW MODEL'
    print, '====================================='
    model = 'POWER'
endif

ifin = where(iiqual_det and ind,nf)
lxraw = dblarr(nsrc)

for i = 0,nfield-1 do begin
    re = execute('ivalid = where(lxraw eq 0. and IIQUAL_DET'+xfield[i]+' and ind)')
    re = execute('lxraw[ivalid] = lx'+xfield[i]+'[ivalid]')
endfor
loglxraw = alog10(lxraw[ifin])
loglxirf = loglxir[ifin]

llsamp = dblarr(1000,nf)
nhsamp = dblarr(1000,nf)
ll_std = dblarr(nf)
ll_mad = dblarr(nf)
nh_std = dblarr(nf)
nh_mad = dblarr(nf)

for i = 0,nf-1 do begin
    llsamp[*,i] = loglxraw[i]-(loglxirf[i]+randomn(seed,1000)*0.3)
    ll_std[i] = stddev(llsamp[*,i])
    ll_mad[i] = medabsdev(llsamp[*,i])
    if (strupcase(model) eq 'POWER') then nhsamp[*,i] = rl2nh(llsamp[*,i],model='POWER') else $
    if (strupcase(model) eq 'BORUS') then nhsamp[*,i] = rl2nh(llsamp[*,i],model='BORUS')
    nh_std[i] = stddev(nhsamp[*,i])
    nh_mad[i] = medabsdev(nhsamp[*,i])
endfor

case strupcase(ret) of
    'LL_STD': return, ll_std
    'LL_MAD': return, ll_mad
    'NH_STD': return, nh_std
    'NH_MAD': return, nh_mad
    ELSE: print, 'NO RETURN CASE INPUT'
endcase


END



