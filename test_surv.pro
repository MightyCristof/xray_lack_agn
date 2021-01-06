



e = {xra:[-3,2],yra:[0.,1.05],stairstep:1,xtitle:'$!8R!7_L$'}
sfa = rl_surv(rl,xd)
sfw = rl_surv(rlw,xdw)
sfr = rl_surv(rlr,xdr)
p = plot(sfr.t,sfr.s,col='blue',_extra=e)
p = plot(sfw.t,sfw.s,col='red',/ov,_extra=e)
p = plot(sfa.t,sfa.s,col='black',/ov,_extra=e)
p.ytitle = 'Survival probability'

sfab = rl_surv(rl,xd,0.1)
sfwb = rl_surv(rlw,xdw,0.1)
sfrb = rl_surv(rlr,xdr,0.1)
p = plot(sfrb.t,sfrb.s,col='blue',_extra=e)
p = plot(sfwb.t,sfwb.s,col='red',/ov,_extra=e)
p = plot(sfab.t,sfab.s,col='black',/ov,_extra=e)
p.ytitle = 'Survival probability (binned)'

END

