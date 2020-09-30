PRO surv_hist, TWOSAMP = twosamp, $
               PLT = plt



;; histogram bin size to match C20
nhbin = 0.4



;; pull bin size
readcol,'wagn.com',col,format='a'
binsz = double(col[9])

;; univariate, KM
;; WISE AGN SOURCES
readcol,'wagn.out',wbinlo,wbinhi,wkm,werr,format='x,d,x,d,d',/silent
;; add first line
wbinhi = [wbinlo[0],wbinhi]
wbinlo = [-3.,wbinlo]
wkm = [1.,wkm]
werr = [0.,werr]
readcol,'wagn.out',wbinc,wd,format='d,d',/silent
;; remove bad line
wbinc = wbinc[1:-1]
wd = wd[1:-1]
;; REMAINING SOURCES
readcol,'ragn.out',rbinlo,rbinhi,rkm,rerr,format='x,d,x,d,d',/silent
;; add first line
rbinhi = [rbinlo[0],rbinhi]
rbinlo = [-3.,rbinlo]
rkm = [1.,rkm]
rerr = [0.,rerr]
readcol,'ragn.out',rbinc,rd,format='d,d',/silent
;; remove bad line
rbinc = rbinc[1:-1]
rd = rd[1:-1]


;; recreate luminosity ratio arrays
;; use central luminosity ratio value and distribute randomly around it using histogram frequency
wrl = []
rrl = []
offset = binsz/2.
for i = 0,n_elements(wbinc)-1 do begin
    if (wd[i] gt 0.) then wrl = [wrl,wbinc[i]+(randomu(seed,round(wd[i]))*binsz-binsz/2.)]
    if (rd[i] gt 0.) then rrl = [rrl,rbinc[i]+(randomu(seed,round(rd[i]))*binsz-binsz/2.)]
endfor

wnh_borus = rl2nh(wrl,model='BORUS')
wnh_power = rl2nh(wrl,model='POWER')
rnh_borus = rl2nh(rrl,model='BORUS')
rnh_power = rl2nh(rrl,model='POWER')


ywnh_borus = histogram(wnh_borus,locations=xwnh_borus,bin=nhbin)
ywnh_power = histogram(wnh_power,locations=xwnh_power,bin=nhbin)
yrnh_borus = histogram(rnh_borus,locations=xrnh_borus,bin=nhbin)
yrnh_power = histogram(rnh_power,locations=xrnh_power,bin=nhbin)

xwnh_borus = [xwnh_borus[0]-nhbin,xwnh_borus,xwnh_borus[-1]+nhbin]
xwnh_power = [xwnh_power[0]-nhbin,xwnh_power,xwnh_power[-1]+nhbin]
xrnh_borus = [xrnh_borus[0]-nhbin,xrnh_borus,xrnh_borus[-1]+nhbin]
xrnh_power = [xrnh_power[0]-nhbin,xrnh_power,xrnh_power[-1]+nhbin]
ywnh_borus = [0.,ywnh_borus,0.]
ywnh_power = [0.,ywnh_power,0.]
yrnh_borus = [0.,yrnh_borus,0.]
yrnh_power = [0.,yrnh_power,0.]



if keyword_set(twosamp) then begin
    
    ;; pull bin size
    readcol,'wagn2.com',col,format='a'
    binsz = double(col[11])
    
    ;; WISE AGN SOURCES
    readcol,'wagn2.out',wbinlo1,wbinhi1,wkm1,werr1,format='x,d,x,d,d,d';,/silent,skipline=wskip11,numline=(wnum11-wskip11)
    ;; split samples
    ii = where(diff(wbinlo1) lt 0.)
    wbinlo2 = wbinlo1[ii+1:-1]
    wbinlo1 = wbinlo1[0:ii]
    wbinhi2 = wbinhi1[ii+1:-1]
    wbinhi1 = wbinhi1[0:ii]
    wkm2 = wkm1[ii+1:-1]
    wkm1 = wkm1[0:ii]
    werr2 = werr1[ii+1:-1]
    werr1 = werr1[0:ii]
    ;; correct first lines
    wbinhi1 = [wbinlo1[0],wbinhi1]
    wbinlo1 = [-3.,wbinlo1]
    wkm1 = [1.,wkm1]
    werr1 = [0.,werr1]
    wbinhi2 = [wbinlo2[0],wbinhi2]
    wbinlo2 = [-3.,wbinlo2]
    wkm2 = [1.,wkm2]
    werr2 = [0.,werr2]
    
    readcol,'wagn2.out',wbinc1,wd1,format='d,d',/silent
    ;; remove bad lines
    wbinc1 = wbinc1[1:-1]
    wd1 = wd1[1:-1]
    ii = (where(diff(wbinc1) lt 0))[0]
    wbinc2 = wbinc1[ii+2:-1]
    wd2 = wd1[ii+2:-1]
    wbinc1 = wbinc1[0:ii]
    wd1 = wd1[0:ii]
    
    ;; REMAINING SOURCES
    readcol,'ragn2.out',rbinlo1,rbinhi1,rkm1,rerr1,format='x,d,x,d,d,d';,/silent,skipline=wskip11,numline=(wnum11-wskip11)
    ;; split samples
    ii = where(diff(rbinlo1) lt 0.)
    rbinlo2 = rbinlo1[ii+1:-1]
    rbinlo1 = rbinlo1[0:ii]
    rbinhi2 = rbinhi1[ii+1:-1]
    rbinhi1 = rbinhi1[0:ii]
    rkm2 = rkm1[ii+1:-1]
    rkm1 = rkm1[0:ii]
    rerr2 = rerr1[ii+1:-1]
    rerr1 = rerr1[0:ii]
    ;; correct first lines
    rbinhi1 = [rbinlo1[0],rbinhi1]
    rbinlo1 = [-3.,rbinlo1]
    rkm1 = [1.,rkm1]
    rerr1 = [0.,rerr1]
    rbinhi2 = [rbinlo2[0],rbinhi2]
    rbinlo2 = [-3.,rbinlo2]
    rkm2 = [1.,rkm2]
    rerr2 = [0.,rerr2]
    
    readcol,'ragn2.out',rbinc1,rd1,format='d,d',/silent
    ;; remove bad lines
    rbinc1 = rbinc1[1:-1]
    rd1 = rd1[1:-1]
    ii = (where(diff(rbinc1) lt 0))[0]
    rbinc2 = rbinc1[ii+2:-1]
    rd2 = rd1[ii+2:-1]
    rbinc1 = rbinc1[0:ii]
    rd1 = rd1[0:ii]
    
    ;; recreate luminosity ratio arrays
    ;; use central luminosity ratio value and distribute randomly around it using histogram frequency
    wrl1 = []
    wrl2 = []
    rrl1 = []
    rrl2 = []
    for i = 0,n_elements(wbinc1)-1 do begin
        if (wd1[i] gt 0.) then wrl1 = [wrl1,wbinc1[i]+(randomu(seed,round(wd1[i]))*binsz-binsz/2.)]
        if (wd2[i] gt 0.) then wrl2 = [wrl2,wbinc2[i]+(randomu(seed,round(wd2[i]))*binsz-binsz/2.)]
        if (rd1[i] gt 0.) then rrl1 = [rrl1,rbinc1[i]+(randomu(seed,round(rd1[i]))*binsz-binsz/2.)]
        if (rd2[i] gt 0.) then rrl2 = [rrl2,rbinc2[i]+(randomu(seed,round(rd2[i]))*binsz-binsz/2.)]
    endfor

    ;; WISE AGN SOURCES
    wnh1_borus = rl2nh(wrl1,model='BORUS')
    wnh1_power = rl2nh(wrl1,model='POWER')
    wnh2_borus = rl2nh(wrl2,model='BORUS')
    wnh2_power = rl2nh(wrl2,model='POWER')
    
    ywnh1_borus = histogram(wnh1_borus,locations=xwnh1_borus,bin=nhbin)
    ywnh1_power = histogram(wnh1_power,locations=xwnh1_power,bin=nhbin)
    ywnh2_borus = histogram(wnh2_borus,locations=xwnh2_borus,bin=nhbin)
    ywnh2_power = histogram(wnh2_power,locations=xwnh2_power,bin=nhbin)

    xwnh1_borus = [xwnh1_borus[0]-nhbin,xwnh1_borus,xwnh1_borus[-1]+nhbin]
    xwnh1_power = [xwnh1_power[0]-nhbin,xwnh1_power,xwnh1_power[-1]+nhbin]
    xwnh2_borus = [xwnh2_borus[0]-nhbin,xwnh2_borus,xwnh2_borus[-1]+nhbin]
    xwnh2_power = [xwnh2_power[0]-nhbin,xwnh2_power,xwnh2_power[-1]+nhbin]
    ywnh1_borus = [0.,ywnh1_borus,0.]
    ywnh1_power = [0.,ywnh1_power,0.]
    ywnh2_borus = [0.,ywnh2_borus,0.]
    ywnh2_power = [0.,ywnh2_power,0.]

    ;; REMAINING SOURCES
    rnh1_borus = rl2nh(rrl1,model='BORUS')
    rnh1_power = rl2nh(rrl1,model='POWER')
    rnh2_borus = rl2nh(rrl2,model='BORUS')
    rnh2_power = rl2nh(rrl2,model='POWER')
    
    yrnh1_borus = histogram(rnh1_borus,locations=xrnh1_borus,bin=nhbin)
    yrnh1_power = histogram(rnh1_power,locations=xrnh1_power,bin=nhbin)
    yrnh2_borus = histogram(rnh2_borus,locations=xrnh2_borus,bin=nhbin)
    yrnh2_power = histogram(rnh2_power,locations=xrnh2_power,bin=nhbin)

    xrnh1_borus = [xrnh1_borus[0]-nhbin,xrnh1_borus,xrnh1_borus[-1]+nhbin]
    xrnh1_power = [xrnh1_power[0]-nhbin,xrnh1_power,xrnh1_power[-1]+nhbin]
    xrnh2_borus = [xrnh2_borus[0]-nhbin,xrnh2_borus,xrnh2_borus[-1]+nhbin]
    xrnh2_power = [xrnh2_power[0]-nhbin,xrnh2_power,xrnh2_power[-1]+nhbin]
    yrnh1_borus = [0.,yrnh1_borus,0.]
    yrnh1_power = [0.,yrnh1_power,0.]
    yrnh2_borus = [0.,yrnh2_borus,0.]
    yrnh2_power = [0.,yrnh2_power,0.]
    
endif


save,file='surv_nh.sav'


if keyword_set(plt) then begin

    e = {xra:[21.,25.5],yra:[0.,100.], $
         stairstep:1,fill_background:1, $
         dimension:[1130,880], $
         font_name:'Times'}
    e.yra = [0.,1.2]
    
    ;;;;;;;;;;
    ;; BORUS
    ;;;;;;;;;;
    nh_bound = alog10([1e21,1.5e24])
    rl_bound = rl2nh(nh_bound,model='BORUS',/lum_out)
    ;; WISE AGN
    pw = plot(xwnh_borus,nm(ywnh_borus),_extra=e,position=[75,445,585,825],/device,/nodata)
    pw.axes[0].showtext=0
    pw.ytickvalues = [0.0:1.2:0.2]
    ;; CT shading
    p = plot([alog10(1.5e24),e.xra[1]],e.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=e.yra[0],fill_color='light grey',fill_transparency=80,/ov)
    ;; data
    hwobs = plot(xwnh_borus,nm(ywnh_borus),_extra=e,col='black',thick=2,linestyle='__',fill_color='grey',fill_transparency=75,/ov,name=' SED unobscured')
    ;; CT lines
    p = plot([1,1]*alog10(1.5e24),[e.yra[1]*0.875,e.yra[1]],'-',thick=2,/ov)
    p = plot([1,1]*alog10(1.5e24),[e.yra[0],e.yra[1]*0.125],'-',thick=2,/ov)
    xmod = text(e.xra[0]+diff(e.xra)/3.,e.yra[1]*0.89,'BORUS model',target=pw,/data,font_size=14,font_name='Times',font_style='Bold')
    ;; REMAINING SOURCES
    pr = plot(xrnh_borus,(yrnh_borus),_extra=e,current=1,position=[585,445,1095,825],/device,/nodata)
    pr.axes[0].showtext=0
    pr.axes[1].showtext=0
    p = plot([alog10(1.5e24),e.xra[1]],e.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=e.yra[0],fill_color='light grey',fill_transparency=80,/ov)
    hr = plot(xrnh_borus,nm(yrnh_borus),_extra=e,col='black',thick=2,linestyle='__',fill_color='grey',fill_transparency=75,/ov,name=' X-ray non-det.')
    ;; CT lines
    p = plot([1,1]*alog10(1.5e24),[e.yra[1]*0.875,e.yra[1]],'-',thick=2,/ov)
    p = plot([1,1]*alog10(1.5e24),[e.yra[0],e.yra[1]*0.125],'-',thick=2,/ov)

    ;;;;;;;;;;;;;;
    ;; POWER LAW
    ;;;;;;;;;;;;;;
    nh_bound = alog10([1e21,1.5e24])
    ll_bound = rl2nh(nh_bound,model='POWER',/lum_out)
    ;; WISE AGN
    pw = plot(xwnh_power,nm(ywnh_power),_extra=e,current=1,position=[75,65,585,445],/device,/nodata)
    pw.xtickvalues = [21.:24.:1.]
    pw.ytickvalues = [0.0:1.:0.2]
    ;; CT shading
    p = plot([alog10(1.5e24),e.xra[1]],e.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=e.yra[0],fill_color='light grey',fill_transparency=80,/ov)
    ;; data
    hw = plot(xwnh_power,nm(ywnh_power),_extra=e,col='black',thick=2,linestyle='__',fill_color='grey',fill_transparency=75,/ov,name=' X-ray non-det.')
    ;; CT lines
    p = plot([1,1]*alog10(1.5e24),[e.yra[1]*0.875,e.yra[1]],'-',thick=2,/ov)
    p = plot([1,1]*alog10(1.5e24),[e.yra[0],e.yra[1]*0.125],'-',thick=2,/ov)
    xmod = text(e.xra[0]+diff(e.xra)/3.,e.yra[1]*0.89,'Power law model',target=pw,/data,font_size=14,font_name='Times',font_style='Bold')
    ;; REMAINING SOURCES
    ;; make plot
    pr = plot(xrnh_power,nm(yrnh_power),_extra=e,current=1,position=[585,65,1095,445],/device,/nodata)
    pr.xtickvalues = [22.:25.:1.]
    pr.axes[1].showtext=0
    ;; CT shading
    p = plot([alog10(1.5e24),e.xra[1]],e.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=e.yra[0],fill_color='light grey',fill_transparency=80,/ov)
    ;; data
    hr = plot(xrnh_power,nm(yrnh_power),_extra=e,col='black',thick=2,linestyle='__',fill_color='grey',fill_transparency=75,/ov,name=' X-ray non-det.')
    ;; CT lines
    p = plot([1,1]*alog10(1.5e24),[e.yra[1]*0.875,e.yra[1]],'-',thick=2,/ov)
    p = plot([1,1]*alog10(1.5e24),[e.yra[0],e.yra[1]*0.125],'-',thick=2,/ov)
    ct = text(24.22,e.yra[1]*0.095,'Compton',col='black',target=pr,font_name='Times',font_style='Bold',/data)          
    ct = text(24.22,e.yra[1]*0.053,'thick',col='black',target=pr,font_name='Times',font_style='Bold',/data)            
    ;; legend
    ;l = legend(target=[hw,hw],/normal,/auto_text_color,sample_width=0.1,horizontal_spacing=0.06,font_name='Times')
    ;l.position = [0.745,0.49]
    ;; axes
    xt = text(0.52,0.03,'$log  !8N!7_{H} [cm^{-2}]$',alignment=0.5,vertical_alignment=0.5,font_size=14,font_name='Times')            
    yt = text(0.015,0.5,'Frequency [normalized]',orientation=90.,alignment=0.5,vertical_alignment=0.5,font_size=14,font_name='Times')
    p.save,'nhdist_km.png'
    
    
    if keyword_set(twosamp) then begin
        ;;;;;;;;;;
        ;; BORUS
        ;;;;;;;;;;
        nh_bound = alog10([1e21,1.5e24])
        rl_bound = rl2nh(nh_bound,model='BORUS',/lum_out)
        ;; WISE AGN
        pw = plot(xwnh2_borus,ywnh2_borus,_extra=e,position=[75,445,585,825],/device,/nodata)
        pw.axes[0].showtext=0
        pw.ytickvalues = [0.0:1.2:0.2]
        ;; CT shading
        p = plot([alog10(1.5e24),e.xra[1]],e.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=e.yra[0],fill_color='light grey',fill_transparency=80,/ov)
        ;; data
        hwobs = plot(xwnh2_borus,nm(ywnh2_borus),_extra=e,col='red',thick=2,linestyle='__',fill_color='red',fill_transparency=75,/ov)
        hwuno = plot(xwnh1_borus,nm(ywnh1_borus),_extra=e,col='blue',thick=2,fill_color='blue',fill_transparency=75,/ov)
        ;; CT lines
        p = plot([1,1]*alog10(1.5e24),[e.yra[1]*0.875,e.yra[1]],'-',thick=2,/ov)
        p = plot([1,1]*alog10(1.5e24),[e.yra[0],e.yra[1]*0.125],'-',thick=2,/ov)
        ;; catalog label
        t = text(e.xra[0]+0.14,e.yra[1]*0.902,'!16WISE!15 AGN',col='red',target=pw,/data,font_size=12,font_name='Times',font_style='Bold')
        t = text(e.xra[0]+0.14,e.yra[1]*0.86,'sources',col='red',target=pw,/data,font_size=12,font_name='Times',font_style='Bold')
        xmod = text(e.xra[0]+diff(e.xra)/3.,e.yra[1]*0.89,'BORUS model',target=pw,/data,font_size=14,font_name='Times',font_style='Bold')
        ;; REMAINING SOURCES
        pr = plot(xrnh2_borus,yrnh2_borus,_extra=e,current=1,position=[585,445,1095,825],/device,/nodata)
        pr.axes[0].showtext=0
        pr.axes[1].showtext=0
        p = plot([alog10(1.5e24),e.xra[1]],e.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=e.yra[0],fill_color='light grey',fill_transparency=80,/ov)
        hrobs = plot(xrnh2_borus,nm(yrnh2_borus),_extra=e,col='red',thick=2,linestyle='__',fill_color='red',fill_transparency=75,/ov)
        hruno = plot(xrnh1_borus,nm(yrnh1_borus),_extra=e,col='blue',thick=2,fill_color='blue',fill_transparency=75,/ov)
        ;; CT lines
        p = plot([1,1]*alog10(1.5e24),[e.yra[1]*0.875,e.yra[1]],'-',thick=2,/ov)
        p = plot([1,1]*alog10(1.5e24),[e.yra[0],e.yra[1]*0.125],'-',thick=2,/ov)
        ;; catalog label
        t = text(e.xra[0]+0.14,e.yra[1]*0.902,'Remaining',col='red',target=pr,/data,font_size=12,font_name='Times',font_style='Bold')
        t = text(e.xra[0]+0.14,e.yra[1]*0.86,'sources',col='red',target=pr,/data,font_size=12,font_name='Times',font_style='Bold')

        ;;;;;;;;;;;;;;
        ;; POWER LAW
        ;;;;;;;;;;;;;;
        nh_bound = alog10([1e21,1.5e24])
        ll_bound = rl2nh(nh_bound,model='POWER',/lum_out)
        ;; WISE AGN
        pw = plot(xwnh2_power,ywnh2_power,_extra=e,current=1,position=[75,65,585,445],/device,/nodata)
        pw.xtickvalues = [21.:24.:1.]
        pw.ytickvalues = [0.0:1.:0.2]
        ;; CT shading
        p = plot([alog10(1.5e24),e.xra[1]],e.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=e.yra[0],fill_color='light grey',fill_transparency=80,/ov)
        ;; data
        hwobs = plot(xwnh2_power,nm(ywnh2_power),_extra=e,col='red',thick=2,linestyle='__',fill_color='red',fill_transparency=75,/ov)
        hwnob = plot(xwnh1_power,nm(ywnh1_power),_extra=e,col='blue',thick=2,fill_color='blue',fill_transparency=75,/ov)
        ;; CT lines
        p = plot([1,1]*alog10(1.5e24),[e.yra[1]*0.875,e.yra[1]],'-',thick=2,/ov)
        p = plot([1,1]*alog10(1.5e24),[e.yra[0],e.yra[1]*0.125],'-',thick=2,/ov)
        xmod = text(e.xra[0]+diff(e.xra)/3.,e.yra[1]*0.89,'Power law model',target=pw,/data,font_size=14,font_name='Times',font_style='Bold')
        ;; REMAINING SOURCES
        ;; make plot
        pr = plot(xrnh2_power,yrnh2_power,_extra=e,current=1,position=[585,65,1095,445],/device,/nodata)
        pr.xtickvalues = [22.:e.xra[1]:1.]
        pr.axes[1].showtext=0
        ;; CT shading
        p = plot([alog10(1.5e24),e.xra[1]],e.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=e.yra[0],fill_color='light grey',fill_transparency=80,/ov)
        ;; data
        hrobs = plot(xrnh2_power,nm(yrnh2_power),_extra=e,col='red',thick=2,linestyle='__',fill_color='red',fill_transparency=75,/ov,name=' SED obscured')
        hruno = plot(xrnh1_power,nm(yrnh1_power),_extra=e,col='blue',thick=2,fill_color='blue',fill_transparency=75,/ov,name=' SED unobscured')
        ;; CT lines
        p = plot([1,1]*alog10(1.5e24),[e.yra[1]*0.875,e.yra[1]],'-',thick=2,/ov)
        p = plot([1,1]*alog10(1.5e24),[e.yra[0],e.yra[1]*0.125],'-',thick=2,/ov)
        ct = text(24.22,e.yra[1]*0.095,'Compton',col='black',target=pr,font_name='Times',font_style='Bold',/data)          
        ct = text(24.22,e.yra[1]*0.053,'thick',col='black',target=pr,font_name='Times',font_style='Bold',/data)            
        ;; legend
        l = legend(target=[hrobs,hruno],/normal,/auto_text_color,sample_width=0.1,horizontal_spacing=0.06,font_name='Times')
        l.position = [0.745,0.49]
        ;; axes
        xt = text(0.52,0.03,'$log  !8N!7_{H} [cm^{-2}]$',alignment=0.5,vertical_alignment=0.5,font_size=14,font_name='Times')            
        yt = text(0.015,0.5,'Frequency [normalized]',orientation=90.,alignment=0.5,vertical_alignment=0.5,font_size=14,font_name='Times')
        p.save,'nhdist_2samp.png'
    endif
endif




END




;  ii = where(iiinf_cha,ncha)
;  data = {sdss_id:0LL,ra:0d,dec:0d,z:0d,zerr:0d,l6um:0d,e_l6um:0d,fx210:0d,e_fx210:0d,iix:0b,iidet:0b,iifinal:0b}
;  data = replicate(data,ncha)
;  tags = ['objid','ra','dec','z','zerr','lir','e_lir','fx_cha','e_fx_cha','iix_cha','iidet_cha','iifinal']
;  for i = 0,n_tags(data)-1 do re = execute('data.(i) = '+tags[i]+'[ii]')
;  mwrfits,data,'src_inf_cha.fits',/create





