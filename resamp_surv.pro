PRO resamp_surv


surv_analysis,fmt='R'
restore,'surv_anal.sav'

rem_col = [5,48,97]
wise_col = [178,24,43]
kmcol = [65,182,196]

;; KM ESTIMATOR
ekm = {xra:[-3.,2.],yra:[0.,1.05], $
       sym_size:1.,sym_filled:0, $
       stairstep:1, $
       xtitle:'$!8R!7_{!8L!7_X}$',ytitle:'KM Estimator', $
       font_name:'Times',font_size:14, $
       buffer:0}
if keyword_set(hide) then ekm.buffer = 1    
p = plot(-abin,akm,'-',_extra=e,/nodata)
pr = plot(-rbin,e_rkmhi,'-',col=rem_col,_extra=ekm,fill_background=1,fill_color=rem_col,fill_level=ekm.yra[0],fill_transparency=75,/ov)
pr = plot(-rbin,e_rkmlo,'-',col=rem_col,_extra=ekm,fill_background=1,fill_color='white',fill_level=ekm.yra[0],/ov)
pr = plot(-rbin,rkm,'-.',col=rem_col,_extra=ekm,thick=2,/ov,name='Residual')
;; white square to remove the extended lines
pclear = plot([ekm.xra[0],min(-rbin)-0.015],[1.,1.],col='white',fill_background=1,fill_color='white',fill_level=ekm.yra[0],/ov)
pw = plot(-wbin,e_wkmhi,'-',col=wise_col,_extra=ekm,fill_background=1,fill_color=wise_col,fill_level=ekm.yra[0],fill_transparency=75,/ov)
pw = plot(-wbin,e_wkmlo,'-',col=wise_col,_extra=ekm,fill_background=1,fill_color='white',fill_level=ekm.yra[0],/ov)
pw = plot(-wbin,wkm,'--',col=wise_col,_extra=ekm,thick=2,/ov,name='!8WISE!7 AGN')
;; white square to remove the extended lines
pclear = plot([ekm.xra[0],min(-wbin)-0.054],[0.6,0.6],col='white',fill_background=1,fill_color='white',fill_level=ekm.yra[0],/ov)
;; legend
l = legend(target=[pw,pr],position=[0.86,0.21],/relative,sample_width=0.15,horizontal_spacing=0.06,font_name='Times',font_size=14)


END