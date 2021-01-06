FUNCTION equal_cenc, xd, rl, obsc, wagn, xdeq=xdeq, rleq=rleq, obsceq=obsceq, wagneq=wagneq


idet = where(xd eq 1,ndet)
inon = where(xd eq 0,nnon)
inon = inon[randomi(ndet,nnon-1,/nodup)]
xd_det = xd[idet]
xd_non = xd[inon]
rl_det = rl[idet]
rl_non = rl[inon]
obsc_det = obsc[idet]
obsc_non = obsc[inon]
wagn_det = wagn[idet]
wagn_non = wagn[inon]
xdeq = [xd_det,xd_non]
rleq = [rl_det,rl_non]
obsceq = [obsc_det,obsc_non]
wagneq = [wagn_det,wagn_non]

return,1



end
