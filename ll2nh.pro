FUNCTION ll2nh, in_arr, $
                LUM_OUT = lum_out, $
                MODEL = model, $
                VERBOSE = verbose


if ~keyword_set(model) then model = 'POWER'
out_arr = dblarr(n_elements(in_arr))-9999.

;; Following NH & LL computed with BORUS
;; NH used in BORUS modeling
nh = [18.:25.5:0.25]

if (strupcase(model) eq 'BORUS') then begin
    if keyword_set(verbose) then print, 'MODEL:  BORUS'
    ll = [0.0000000,-1.7494337e-06,-4.8438337e-06,-1.0352710e-05,-2.0148913e-05,-3.7565625e-05,-6.8537547e-05,-0.00012360928,-0.00022152840,-0.00039560872,-0.00070502834,-0.0012548026,-0.0022310111,-0.0039624272,-0.0070270606,-0.012432091,-0.021905311,-0.037425573,-0.063707769,-0.10707675,-0.17545840,-0.27836069,-0.42899658,-0.65166572,-0.98144732,-1.4321754,-1.8510858,-1.9891130,-2.0028169,-2.0046453,-2.0052188]
endif else if (strupcase(model) eq 'POWER') then begin
    if keyword_set(verbose) then print, 'MODEL:  POWLAW'
    ll = [0.0000000,-3.4647453e-06,-9.6075680e-06,-2.0538313e-05,-3.9974826e-05,-7.4533979e-05,-0.00013598264,-0.00024523369,-0.00043943864,-0.00078456120,-0.0013975610,-0.0024853603,-0.0044125784,-0.0078171388,-0.013801024,-0.024225027,-0.042105928,-0.071991821,-0.11992454,-0.19252437,-0.29635834,-0.44201558,-0.65346850,-0.97697257,-1.4932660,-2.3383189,-3.7380739,-6.0811486,-10.059617,-16.908176,-28.845231]
endif

;; Calculate either NH->L/L -OR- L/L->NH
if keyword_set(LUM_OUT) then begin
    ;; Calculate Lx/Lx(Lir)
    ;; need separate WHERE() calls; iiunb =/= ~iibnd due to -9999
    iibnd = (in_arr ge nh[0] and in_arr le nh[-1])
    iiunb = (in_arr lt nh[0] or in_arr gt nh[-1]) and in_arr ne -9999.
    ibnd = where(iibnd,ctbnd)
    iunb = where(iiunb,ctunb)
    if (ctbnd gt 0.) then out_arr[ibnd] = interpol(ll,nh,in_arr[ibnd])
    if (ctunb gt 0.) then begin
        out_arr[iunb[where(in_arr[iunb] lt nh[0],/NULL)]] = 0.
        out_arr[iunb[where(in_arr[iunb] gt nh[-1],/NULL)]] = -2.5
    endif
endif else begin
    ;; Calculate NH
    ;; need separate WHERE() calls; iiunb =/= ~iibnd due to -9999
    iibnd = (in_arr le ll[0] and in_arr ge ll[-1])
    iiunb = (in_arr gt ll[0] or in_arr lt ll[-1]) and in_arr ne -9999.
    ibnd = where(iibnd,ctbnd)
    iunb = where(iiunb,ctunb)
    if (ctbnd gt 0.) then out_arr[ibnd] = interpol(nh,ll,in_arr[ibnd])
    if (ctunb gt 0.) then begin
        out_arr[iunb[where(in_arr[iunb] gt ll[0],/NULL)]] = 18.
        out_arr[iunb[where(in_arr[iunb] lt ll[-1],/NULL)]] = 26.
    endif
endelse

return, out_arr


END





