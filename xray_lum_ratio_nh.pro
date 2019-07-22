PRO xray_lum_ratio_nh


common _inf_fits
common _inf_nst
common _inf_xmm
common _inf_cha
common _det_nst
common _det_xmm
common _det_cha
common _det_wac
common _soft210
common _flx_lim
common _agn_lum
;common _comp







---------------------------------------------------------------------------------------
;; NH DISTRIBUTIONS -- RAW UNMODELED OBSCURATION
;;----------------------------------------------------------------------------------------
nhdet_210 = 'NHDET'+xfield+'_210'
nhlim_210 = 'NHLIM'+xfield+'_210'
eband = '2-10'
for i = 1,nfield-1 do begin
    re = execute(nhdet_210[i]+' = dblarr(nsrc)')
    re = execute(nhlim_210[i]+' = dblarr(nsrc)')
    re = execute(nhdet_210[i]+'['+iagn_det_210[i]+'] = ll2nh('+lldet_210[i]+'['+iagn_det_210[i]+'],eband)')
    re = execute(nhlim_210[i]+'['+iagn_lim_210[i]+'] = ll2nh('+lllim_210[i]+'['+iagn_lim_210[i]+'],eband)')
endfor

nhdet_1040 = 'NHDET'+xfield+'_1040'
nhlim_1040 = 'NHLIM'+xfield+'_1040'
eband = '10-40'
for i = 0,nfield-1 do begin
    re = execute(nhdet_1040[i]+' = dblarr(nsrc)')
    re = execute(nhlim_1040[i]+' = dblarr(nsrc)')
    re = execute(nhdet_1040[i]+'['+iagn_det[i]+'] = ll2nh('+lldet_1040[i]+'['+iagn_det[i]+'],eband)')
    re = execute(nhlim_1040[i]+'['+iagn_lim[i]+'] = ll2nh('+lllim_1040[i]+'['+iagn_lim[i]+'],eband)')
endfor

sav_vars = [sav_vars,'NHDET_1040','NHLIM_1040',nhdet_1040,nhlim_1040]
sav_inds = [sav_inds]
endif









;;----------------------------------------------------------------------------------------
;; PLOTS!!
;;----------------------------------------------------------------------------------------
if keyword_set(plt) then begin
    
    
    
    
    nh = [18.0000,18.2500,18.5000,18.7500,19.0000,19.2500,19.5000,19.7500,20.0000,20.2500,20.5000,20.7500,21.0000,21.2500,21.5000,21.7500,22.0000,22.2500,22.5000,22.7500,23.0000,23.2500,23.5000,23.7500,24.0000,24.2500,24.5000,24.7500,25.0000,25.2500,25.5000]
    ll = [0.0000000,-1.5057464e-06,-4.1664155e-06,-8.9047723e-06,-1.7329840e-05,-3.2310490e-05,-5.8948438e-05,-0.00010631492,-0.00019053757,-0.00034027833,-0.00060646268,-0.0010795051,-0.0019197445,-0.0034108858,-0.0060529924,-0.010721430,-0.018930153,-0.032439742,-0.055555079,-0.094315205,-0.15685371,-0.25342692,-0.39729097,-0.61041384,-0.92494804,-1.3586775,-1.7994156,-1.9806261,-2.0018679,-2.0042227,-2.0049609]    ;nh_lines = alog10([1e20,5e20,1e21,5e21,1e22,5e22,1e23,5e23,1e24,5e24,1e25,5e25])
    nh_lines = [5e23,1e24,5e24,1e25]
    ll_lines = transpose(rebin(interpol(ll,nh,alog10(nh_lines)),n_elements(nh_lines),2))

    ebv_lines = [-4,-3]
    
    e = {symbol:'o',sym_size:0.5,sym_filled:1,linestyle:'', $
         ;xra:[5e-3,1e2],xlog:1,yra:[-4,2],ylog:0, $
         xra:[-4.5,2.5],xlog:0,yra:[-3,1],ylog:0, $
         xtitle:'$log  \itE\rm(B-V)_{AGN}$',ytitle:'$log  \itL\rm_{2-10keV} / \itL\rm_{2-10keV}(\itL\rm_{IR})$', $
         buffer:0}
    pwlim = plot(lebv[iwlim],lll_210[iwlim],col='orange',_extra=e,name='AGN w.o X-ray (limits)',transparency=75)
    pwdet = plot(lebv[iwdet],lld_210[iwdet],col='dodger blue',/ov,_extra=e,name='AGN w. X-ray (detections)',transparency=50)
    p = plot(e.xra,[0,0],'--',/ov)
    nh_str = '$'+['5\times10^{23}','1\times10^{24}','5\times10^{24}','1\times10^{25}']+'$'
    ll_text = ll_lines[0,*]+[0.08,-0.22,0.08,-0.22]
    t = text(ebv_lines[0]-0.3,ll_text[0]+0.25,'$\itN\rm_H / cm^{-2}$',col='red',/data)
    for i = 0,n_elements(nh_lines)-1 do t = text(ebv_lines[0]+0.1,ll_text[i],nh_str[i],col='red',/data)
    for i = 0,n_elements(nh_lines)-1 do p = plot(ebv_lines,ll_lines[*,i],':r',thick=2,/ov)
    leg = legend(target=[pwdet,pwlim],position=[0.41,0.95],/auto_text_color,sample_width=0,font_size=12,/relative)
    
    
    
    
    
    
    
    plim = plot(lebv[ixlim],lll_210[ixlim],col='orange',_extra=e,name='AGN X-ray non-det.')
    pdet = plot(lebv[ixdet],lld_210[ixdet],col='dodger blue',/ov,_extra=e,name='AGN + X-ray')
    ;pwac = plot(lebv[iwlim],lllim[iwlim],'or',sym_size=0.25,sym_filled=1,/ov)
    ;plim.save,'ll_vs_ebv.png'

    
    
    
    
    
    e = {symbol:'o',sym_size:0.5,sym_filled:1,sym_transparency:50,linestyle:'', $
         ;xra:[5e-3,1e2],xlog:1,yra:[-4,2],ylog:0, $
         xra:[-3.5,2.5],xlog:0,yra:[-2,1],ylog:0, $
         xtitle:'$log E(B-V)_{AGN}$',ytitle:'$log( L_{2-10} keV} / L_{2-10 keV}(L_{IR}) )$', $
         buffer:0}
    plim = plot(lebv[ixlim],lllim_1040[ixlim],col='orange',_extra=e,name='AGN X-ray non-det.')
    pdet = plot(lebv[ixdet],lldet_1040[ixdet],col='dodger blue',/ov,_extra=e,name='AGN + X-ray')
    ;pwac = plot(lebv[iwlim],lllim[iwlim],'or',sym_size=0.25,sym_filled=1,/ov)
    ;plim.save,'ll_vs_ebv.png'
    plim = plot(lebv[iwlim],lllim_1040[iwlim],col='orange',_extra=e,name='AGN X-ray non-det.')
    pdet = plot(lebv[iwdet],lldet_1040[iwdet],col='dodger blue',/ov,_extra=e,name='AGN + X-ray')
    ;plim.save,'ll_vs_ebv_wac.png'
    
stop

    xnhlines = rebin([[-3.3],[-2.8]],4,2)
    ynhlines = rebin([-0.264805,-1.03190,-1.52835,-2.06607],4,2)
    ;ynhlines = rebin([-0.924948,-1.94364,-2.00187,-2.00555],4,2)
    nhlabel = ['$N_H=10^{24}cm^{-2}$','$    =5x10^{24}cm^{-2}$','$    =10^{25}cm^{-2}$','$    =5x10^{25}cm^{-2}$']
    for i = 0,3 do p = plot(xnhlines[i,*],ynhlines[i,*],':',thick=4,/ov)
    e = {font_size:8,font_style:'bold',data:1}
    for i = 0,3 do t = text(xnhlines[i,0]-0.05,ynhlines[i,0]-0.15,nhlabel[i],_extra=e)

    
endif



illlw = where(lll_210 ne -9999. and iiwac)
illxw = where(llx_210 ne -9999. and iiwac)
illl = where(lll_210 ne -9999.)
illx = where(llx_210 ne -9999.)

e = {symbol:'o',sym_size:0.5,sym_filled:1,sym_transparency:50,linestyle:'', $
     xra:[-3.5,2.5],xlog:0,yra:[-3,1],ylog:0, $
     ytitle:'$log  \itL\rm_{2-10 keV} / \itL\rm_{2-10 keV}(\itL\rm_{IR}) $', $
     font_size:14., $
     dimension:[1140,510],buffer:1}
pos = [[65,65,570,445],[570,65,1075,445]]

     
plim = plot(lebv[illlw],lll_210[illlw],col='orange',_extra=e,name='AGN X-ray non-det.',position=pos[*,0],/device)
pdet = plot(lebv[illxw],llx_210[illxw],col='dodger blue',/ov,_extra=e,name='AGN + X-ray')
pzero = plot(e.xra,[0.,0.],'--',/ov)
lab = text(0.05,0.12,'$\itWISE\rm  AGN$',/relative,target=plim,font_size=12)
plim = plot(lebv[illl],lll_210[illl],col='orange',_extra=e,name='AGN X-ray non-det.',position=pos[*,1],/device,/current)
pdet = plot(lebv[illx],llx_210[illx],col='dodger blue',/ov,_extra=e,name='AGN + X-ray')
pzero = plot(e.xra,[0.,0.],'--',/ov)
plim.axes[1].showtext=0
lab = text(0.05,0.12,'$All sample AGN$',/relative,target=plim,font_size=12)
xt = text(0.5,0.03,'$log  \itE\rm(\itB\rm-\itV\rm)_{AGN}$',alignment=0.5,/relative,font_size=14)
plim.save,'ll_vs_ebv.png'


END










