FUNCTION rl2nh, in_arr, $
                MODEL = model, $
                RL_OUT = rl_out, $
                SCAT = scat
                

if ~keyword_set(model) then begin
    print, '====================================='
    print, 'NO MODEL SET. RUNNING POWER LAW MODEL'
    print, '====================================='
    model = 'POWER'
endif
out_arr = dblarr(n_elements(in_arr))-9999.

xnh = [21.:25.:0.25];5:0.25]
;; choose XSPEC model
if (strupcase(model) eq 'BORUS') then begin
    xrl = [0.0000000,-0.0033660974,-0.0092815038,-0.019583447,-0.037246642, $
           -0.065215346,-0.10966440,-0.17617641,-0.26929564,-0.39619911, $
           -0.57459371,-0.83845970,-1.2304071,-1.7084625,-1.9592897, $
           -1.9936579,-1.9992526];,-2.0013390,-2.0017935]
endif else if (strupcase(model) eq 'POWER') then begin
    xrl = [0.0000000,-0.0034045604,-0.0093884460,-0.019812448,-0.037693350, $
           -0.067579243,-0.11551196,-0.18811179,-0.29194576,-0.43760300, $
           -0.64905592,-0.97255999,-1.4888534,-2.3339063,-3.7336613, $
           -6.0767360,-10.055204];,-16.903763,-28.840818]
endif

;; Calculate either NH->L/L -OR- L/L->NH
if keyword_set(rl_out) then begin
    ;; Calculate Lx/Lx(Lir)
    igd = where(in_arr ge xnh[0] and in_arr le xnh[-1],gdct)
    if (gdct gt 0.) then out_arr[igd] = interpol(xrl,xnh,in_arr[igd])
    ihi = where(in_arr gt xnh[-1],hict)
    if (hict gt 0.) then out_arr[ihi] = xrl[-1]
    ilo = where(in_arr lt xnh[0] and in_arr ne -9999.,loct)
    if (loct gt 0.) then out_arr[ilo] = xrl[0]
endif else begin
    ;; Calculate NH
    igd = where(in_arr ge xrl[-1] and in_arr le xrl[0],gdct)
    if (gdct gt 0.) then out_arr[igd] = interpol(xnh,xrl,in_arr[igd])
    ihi = where(in_arr gt xrl[0],hict)
    if (hict gt 0.) then out_arr[ihi] = xnh[0]
    ilo = where(in_arr lt xrl[-1] and in_arr ne -9999.,loct)
    if (loct gt 0.) then out_arr[ilo] = xnh[-1]
endelse

;; add scatter
if (n_elements(scat) gt 0) then out_arr += randomn(seed,n_elements(out_arr))*scat

return, out_arr


END





