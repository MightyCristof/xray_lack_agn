PRO plot_xray_lack_agn, PHOT_SPEC = phot_spec, $
                        ZHIST = zhist, $
                        FLUX_LIMIT = flux_limit, $
                        SEDS = seds, $
                        LX_LIR = lx_lir, $
                        LUM_RATIO = lum_ratio, $
                        NH_DIST = nh_dist, $
                        HIDE = hide, $
                        NEW = new, $
                        SAV = sav, $
                        LOW_RES = low_res       


;; load data
common _fits    
common _resamp
common _comp    
common _inf_cha 
common _inf_xmm 
common _inf_nst 
common _det_cha 
common _det_xmm 
common _det_nst 
common _det_wac 
common _xconv   
common _fxlim 
common _agnlum 
;common _clean_cha
;common _clean_xmm
;common _clean_nst
common _quality  
common _combined
common _nhdist


if keyword_set(low_res) then res = 100 else res = 600

file_mkdir,'figures'
;;----------------------------------------------------------------------------------------
;; z-PHOT vs z-SPEC
;;----------------------------------------------------------------------------------------
if keyword_set(phot_spec) then begin
    print, '    PREPARING REDSHIFT COMPARISON'
    
    file = file_search()
    ifile = where(strmatch(file,'phot_spec_srcs.sav'),nfile)
    if keyword_set(new) then nfile = 0
    if (nfile eq 0) then begin
        ;; AGES
        dir_ages = '/Users/ccarroll/Research/surveys/NDWFS/AGES-survey-J_ApJS_200_8_sources.fits'
        ages = mrdfits(dir_ages,1)
        ;; has spec z
        izs = where(finite(ages.z1),zslen)
        if (zslen gt 0) then ages = ages[izs]
        
        ;; SDSS DR12
        dir_sdss = '/Users/ccarroll/Research/surveys/SDSS/DR14/sdss-dr14-cat-part'+['37','38']+'.fits.gz'
        sdss = [mrdfits(dir_sdss[0],1),mrdfits(dir_sdss[1],1)]
        ;; has photo z
        izp = where(finite(sdss.zp) and sdss.zp ne -9999.,zplen)
        if (zplen gt 0) then sdss = sdss[izp]
    
        ;; match
        spherematch,sdss.ra_sdss,sdss.dec_sdss,ages._raj2000,ages._dej2000,3./3600.,isdss,iages,sep
        ages = ages[iages]
        sdss = sdss[isdss]
        sep *= 3600.
        
        isep = where(sep lt 1.,ct)
        if (ct gt 0.) then begin
            ages = ages[isep]
            sdss = sdss[isep]
            sep  = sep[isep]
        endif
        
        save,ages,sdss,sep,/compress,file='phot_spec_srcs.sav'
    endif else restore, file[ifile]
    ;; quality cuts for catalog sources
    iiqages = (ages.gal eq 1 or ages.qso eq 1 or ages.agn eq 1) and ages.s_n1 gt 3. and ages.z1 gt 0.
    iiqsdss = sdss.zp gt 0. and sdss.zp le 0.8 and sdss.photoerrorclass ge -1 and sdss.photoerrorclass le 3
    iq = where(iiqages and iiqsdss)    
    qages = ages[iq]
    qsdss = sdss[iq]

    ;; separate sample from non-sample
    spherematch,ra,dec,qsdss.ra_sdss,qsdss.dec_sdss,1./3600.,isamp,isdss,sep_samp
    iiin = bytarr(n_elements(qsdss))
    iiin[isdss] = 1
    iin = where(iiin eq 1)
    iex = where(iiin eq 0)
    ages_ex = qages[iex]
    sdss_ex = qsdss[iex]
    ages_in = qages[iin]
    sdss_in = qsdss[iin]

    ;; AGES z-spec to replace SDSS z-phot in sample
    ;zages = dblarr(nsrc)
    ;zages[isamp] = ages_in.z1

    ;; variables
    zs_ex = ages_ex.z1
    zp_ex = sdss_ex.zp
    delz_ex = (zp_ex-zs_ex)/(1+zs_ex)
    zexsig = stddev(delz_ex)
    nex = n_elements(zs_ex)
    zs_in = ages_in.z1
    zp_in = sdss_in.zp
    delz_in = (zp_in-zs_in)/(1+zs_in)
    zinsig = stddev(delz_in)
    nin = n_elements(zs_in)
    
    ;; PLOT ZP vs ZS
    e = {xra:[0.,1.],yra:[0.,1.],aspect_ratio:1, $
         sym_size:0.5,sym_filled:1,transparency:75, $
         xtitle:'$!8z!7_{spec} (VLT)$',ytitle:'$!8z!7_{phot} (SDSS)$', $
         font_name:'Times',font_size:14, $
         dimension:[700,800],buffer:0}
    if keyword_set(hide) then e.buffer = 1
    col_ex = [0,114,178]
    col_in = [240,228,66]
    ;; [204,121,167],[213,94,0],[0,158,115],[0,114,178],[240,228,66]
    p = plot(zs_ex,zp_ex,'o',color=col_ex,_extra=e)
    p.position = [0.12735492,0.33393553,0.769520,0.83583007]
    p = plot([0.,1.],[0.,1.],'--',thick=2,/ov)
    p = plot(zs_in,zp_in,'o',color=col_in,sym_size=0.5,sym_filled=1,/ov)
    p = plot(zs_in,zp_in,'o',color=black,sym_size=0.5,sym_thick=1.,sym_filled=0,transparency=75,/ov)
    iexoff = where(zs_ex gt e.xra[1],nexoff)
    iinoff = where(zs_in gt e.xra[1],ninoff)
    p.axes[0].showtext = 0
    num_str = '$Num. matches = '+commas(strtrim(nex,2))+' ('+commas(strtrim(nin,2))+')'+'$'
    t = text(target=p,e.xra[1]*0.95,e.yra[1]*0.05,num_str,/data,alignment=1,font_size=14,font_style='Bold',font_name='Times',fill_background=1,fill_color='white')
    ;; PLOT RESIDUALS
    er = {xra:e.xra,yra:[-0.8,0.8],aspect_ratio:0, $
          sym_size:0.5,sym_filled:1,transparency:75, $
          xtitle:'$!8z!7_{spec} (MMT)$',ytitle:'$\Delta!8z!7 / (1+!8z!7_{spec})$', $
          font_name:'Times',font_size:14, $
          dimension:[700,800],buffer:0}
    pr = plot(zs_ex,delz_ex,'o',color=col_ex,_extra=er,/current)
    pr.position = [p.position[0],0.13195312,p.position[2],p.position[1]]
    pr = plot(p.xra,[0,0],'--',thick=2,/ov)
    pr = plot(zs_in,delz_in,'o',color=col_in,sym_size=0.5,sym_filled=1,/ov)
    pr = plot(zs_in,delz_in,'o',color=black,sym_size=0.5,sym_thick=1.,sym_filled=0,transparency=75,/ov)
    pr.ytickvalues = [-0.4,0.,0.4]
    pr.yminor = 3.
    zsig_str = '$\sigma_{\Delta z} = '+string(rnd(zexsig,5),format='(d5.3)')+'('+string(rnd(zinsig,5),format='(d5.3)')+')'+'$'
    t = text(target=pr,er.xra[1]*0.95,er.yra[1]*0.62,zsig_str,/data,alignment=1,font_size=14,font_style='Bold',font_name='Times',fill_background=1,fill_color='white')
    ;; PLOT RESIDUALS HISTOGRAM
    yh_ex = histogram(delz_ex,bin=scott(delz_ex),locations=xh_ex)
    yh_in = histogram(delz_in,bin=scott(delz_in),locations=xh_in)
    eh = {xra:[0.,ceil(max(yh_ex)/100.)*100.],yra:[-0.3,0.3], $
          stairstep:1,fill_background:1,fill_transparency:10, $
          xtitle:'Frequency', $
          font_name:'Times',font_size:14}
    ph = plot(yh_ex,xh_ex,/current,_extra=eh,fill_color=col_ex)
    ph = plot(yh_in,xh_in,/ov,_extra=eh,fill_color=col_in)
    ph = plot(ph.xra,[0.,0.],'--',thick=2,/ov)
    ph.position = [p.position[2],pr.position[1],0.91303385,pr.position[3]]
    ph.axes[1].showtext = 0
    ph.xtickvalues = [ph.xra[1]/2.,ph.xra[1]]
    ph.xminor = 3.
    ;; PLOT ZS HISTOGRAM
    zbin = 0.025;scott(zs)>scott(zp)
    yhzs_ex = histogram(zs_ex,bin=zbin,locations=xhzs_ex,min=0.,max=1.)
    yhzs_in = histogram(zs_in,bin=zbin,locations=xhzs_in,min=0.,max=1.)
    xhzs_ex = xhzs_ex+zbin/2.
    xhzs_in = xhzs_in+zbin/2.
    eh_zs = {xra:e.xra,yra:eh.xra, $
            stairstep:1,fill_background:1,fill_transparency:10, $
            ytitle:'Frequency', $
            font_name:'Times',font_size:14}
    ph_zs = plot(xhzs_ex,yhzs_ex,/current,_extra=eh_zs,fill_color=col_ex)
    ph_zs = plot(xhzs_in,yhzs_in,/ov,_extra=eh_zs,fill_color=col_in)
    ph_zs.position = [p.position[0],p.position[3],p.position[2],0.95303385]
    ph_zs.axes[0].showtext = 0
    ph_zs.ytickvalues = ph.xtickvalues
    ph_zs.yminor = 3.
    ;; PLOT ZP HISTOGRAM
    yhzp_ex = histogram(zp_ex,bin=zbin,locations=xhzp_ex,min=0.,max=1.)
    yhzp_in = histogram(zp_in,bin=zbin,locations=xhzp_in,min=0.,max=1.)
    ;xhzp_ex = xhzp_ex+zbin/2.
    ;xhzp_in = xhzp_in+zbin/2.
    eh_zp = {xra:eh.xra,yra:e.yra, $
            stairstep:1,fill_background:1,fill_transparency:10, $
            font_name:'Times',font_size:14}
    ph_zp = plot(yhzp_ex,xhzp_ex,/current,_extra=eh_zp,fill_color=col_ex)
    ph_zp = plot(yhzp_in,xhzp_in,/ov,_extra=eh_zp,fill_color=col_in)
    ph_zp.position = [p.position[2],p.position[1],ph.position[2],p.position[3]]
    ph_zp.axes[1].showtext = 0
    ph_zp.axes[0].showtext = 0
    ph_zp.xtickvalues = ph.xtickvalues
    ph_zp.xminor = 3.
        
    print, 'NUM. SDSS:        ' + strtrim(nex,2)
    print, 'NUM. SAMP:        ' + strtrim(nin,2)
    print, 'NUM. SDSS MISSES: ' + strtrim(nexoff,2)
    print, 'NUM. SAMP MISSES: ' + strtrim(ninoff,2)
    print, 'SIGMA ZSDSS:      ' + strtrim(zexsig,2)
    print, 'SIGMA ZSAMP:      ' + strtrim(zinsig,2)
        
    save,nex,nin,nexoff,ninoff,zexsig,zinsig,file='phot_spec_errs.sav'
    ;save,zages,file='zages.sav'

    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/phot_v_spec.eps',/BITMAP else $
                                                     p.save,'figures/phot_v_spec.png',resolution=res
    endif
endif



;;----------------------------------------------------------------------------------------
;; z HIST
;;----------------------------------------------------------------------------------------
if keyword_set(zhist) then begin

    izs = where(strmatch(ztype,'ZS*'))
    izp = where(strmatch(ztype,'ZP'))
    izx = where(strmatch(ztype,'PEAKZ'))
    
    zbin = 0.02d
    yz = [0.,histogram(z,bin=zbin,location=xz,min=0.,max=max(z)),0.]
    yzs = [0.,histogram(z[izs],bin=zbin,location=xzs,min=0.,max=max(z)),0.]
    yzp = [0.,histogram(z[izp],bin=zbin,location=xzp,min=0.,max=max(z)),0.]
    yzx = [0.,histogram(z[izx],bin=zbin,location=xzx,min=0.,max=max(z)),0.]
    xz = [-0.02d,xz,max(z)+zbin]+zbin/2.
    xzs = [-0.02d,xzs,max(z)+zbin]+zbin/2.
    xzp = [-0.02d,xzp,max(z)+zbin]+zbin/2.
    xzx = [-0.02d,xzx,max(z)+zbin]+zbin/2.

    col = [[27,158,119],[217,95,2],[117,112,179],[231,41,138],[102,166,30]]
    col = col[*,[0,2,3]]
    ;col = [[227,26,28],[166,206,227]]
    
    e = {xtitle:'!8z!7',ytitle:'Frequency', $
         font_size:14,font_name:'Times', $
         xra:[-zbin,max(z)+zbin], $
         yra:[0.,ceil((max(yz))/100.)*100.], $
         stairstep:1,thick:2, $
         fill_background:1,fill_transparency:75}
    p = plot(xz,yz,_extra=e,linestyle='-',/nodata,yra=[0,2300])
    pzp = plot(xzp,yzp,_extra=e,col=col[*,0],fill_color=col[*,0],linestyle='__',/ov,name='!8z!7$_{phot}$ (SDSS DR12)')    
    pzs = plot(xzs,yzs,_extra=e,col=col[*,1],fill_color=col[*,1],linestyle='-',/ov,name='!8z!7$_{spec}$') 
    pzx = plot(xzx,yzx,_extra=e,col=col[*,2],fill_color=col[*,2],linestyle='-.',/ov,name='!8z!7$_{phot}$ (!8XDQSOz!7)')
    l = legend(target=[pzs,pzp,pzx],position=[0.86,0.95],/relative,/auto_text_color,sample_width=0.16,horizontal_spacing=0.06,font_name='Times',font_size=14)
    
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/zhist.eps',/BITMAP else $
                                                     p.save,'figures/zhist.png',resolution=res
    endif
endif



;;----------------------------------------------------------------------------------------
;; Flux Limit
;;----------------------------------------------------------------------------------------
if keyword_set(flux_limit) then begin
    multi_sn = 1
    print, '    PREPARING PAPER FLUX LIMIT'
    dim = [1320,510]
    src_pos = [[80,75,465,455],[465,75,850,455],[850,75,1235,455]]
    axis_src = 0
    which_axis = 1
    leg_src = nfield-1
    leg_pos = [0.755,0.78]
    xpos = [0.5,0.035]
    ypos = [0.025,0.5]
    
    e = {yra:[-15.5,-9.], $
         xlog:1,ylog:0, $
         xtitle:'$exposure time [s]$',ytitle:'$log  !8F!7_{2-10 keV} [erg s^{-1}cm^{-2}]$', $
         font_name:'Times',font_size:14, $ $
         buffer:0}
    if keyword_set(hide) then e.buffer = 1
;    xra = [[8e2,max(texp_cha)],[8e2,max(texp_xmm)],[2e3,max(texp_nst)>3e5]]
    xra = [[8e2,3e5],[8e2,3e5],[2e3,3e5]]
    
    cat = ['Chandra','XMM-Newton','NuSTAR']
    
    for i = 0,nfield-1 do begin
        if (i eq 0) then current = 0 else current = 1
        re = execute('p = plot(CAT_EXP'+xfield[i]+',alog10(CAT_FLX'+xfield[i]+'),_extra=e,xra=xra[*,i],dimension=dim,position=src_pos[*,i],/device,current=current,/nodata)')
        if (i ne 0) then p.axes[1].showtext=0
        if (i ne 1) then p.xtitle = ''
        pcat_dark = plot([1e-5,1e-5],[1e5,1e6],'s',col='dark grey',sym_filled=1,sym_size=0.5,transparency=75,/ov,name='X-ray Catalog(s)')
        
        ;; average error bar
        xpt = [1.5e3,1.6e3,3.63e3]
        ypt = [-14.7]
        perr = plot([xpt[i]],ypt,'s',sym_size=4,sym_thick=2,/ov)
        re = execute('yer = [median(CAT_ERR'+xfield[i]+'/(alog(10.)*CAT_FLX'+xfield[i]+'))]')
        perr = errorplot([xpt[i]],ypt,yer,errorbar_capsize=0.1,linestyle='',/ov)
        perr = plot([xpt[i]],ypt,'S',col='white',sym_size=1.,sym_filled=1,/ov)
        perr = plot([xpt[i]],ypt,'S',col='dodger blue',sym_size=1.,sym_filled=0,/ov)
        ;if (i eq 2) then terr = text(xpt[i]*1.4,ypt,target=perr,/data,'$Median \sigma_{log !8F!7}$',vertical_alignment=0.5,font_size=12,font_name='Times')
        ;if (i ne axis_src) then p.axes[which_axis].showtext = 0
        t = text(0.06,0.88,cat[i],target=p,/RELATIVE,font_size=16,font_style='Bold italic',alignment=0.,font_name='Times')

        re = execute('ncat = n_elements(CAT_EXP'+xfield[i]+')')
        ;if (ncat gt 1e5) then re = execute('irand = round(randomu(seed,n_elements(CAT_EXP'+xfield[i]+')/10)*n_elements(CAT_EXP'+xfield[i]+'))') else $
                              re = execute('irand = lindgen(n_elements(CAT_EXP'+xfield[i]+'))')
        re = execute('pcat_lite = plot(CAT_EXP'+xfield[i]+'[irand],alog10(CAT_FLX'+xfield[i]+'[irand]),"s",col="light grey",sym_filled=1,sym_size=0.5,transparency=85,/ov)')
        ; plot open symbols first (makes for a cleaner plot)
        re = execute('pnon = plot(TEXP'+xfield[i]+'[where(iiqual_non'+xfield[i]+')],logfxir[where(iiqual_non'+xfield[i]+')],"o",sym_filled=1,sym_size=0.5,color="orange",transparency=50,/ov,name="X-ray non-det. [$!8F!7_{X,expected}$]")')
        ;; off plot range
        ;re = execute('ioff = where(iiqual_non'+xfield[i]+' and TEXP'+xfield[i]+' gt xra[1,i],noff)')
        ;if (noff gt 0) then pwoff = arrow(transpose([[make_array(noff,value=xra[1,i]*0.75)],[make_array(noff,value=xra[1,i])]]),transpose([[logfxir[ioff]],[logfxir[ioff]]]),color='orange',/ov,/data,thick=2,head_size=0.5,target=pnon)
        ;; stretch X-ray flux limit for plot
        re = execute('lxr = alog10([min(EXP'+xfield[i]+'[where(EXP'+xfield[i]+' gt 0.)])<min(TEXP'+xfield[i]+'[where(TEXP'+xfield[i]+' gt 0.)]),max(EXP'+xfield[i]+'[where(EXP'+xfield[i]+' gt 0.)])>max(TEXP'+xfield[i]+'[where(TEXP'+xfield[i]+' gt 0.)])])')
        xx = [lxr[0]:lxr[1]:diff(minmax(lxr))/1000.]
        xx_str = '['+strtrim(indgen(degr[i]+1),2)+']*xx^'+strtrim(indgen(degr[i]+1),2)+'.'
        re = execute('fxlim = FXLIM_CS'+xfield[i])
        re = execute('yy = '+strjoin(+'fxlim'+xx_str,' + '))
        pxlim = plot(10.^xx,yy,'--r',thick=2,/ov,name="X-ray flux limit")
        if keyword_set(multi_sn) then begin
            re = execute('ilim = where(CAT_EXP5'+xfield[i]+' ge min(10.^xx) and CAT_EXP5'+xfield[i]+' le max(10.^xx))')
            re = execute('p = plot(CAT_EXP5'+xfield[i]+'[ilim],alog10(CAT_LIM5'+xfield[i]+'[ilim]),"--",thick=2,/ov)')
            if (i eq nfield-1) then begin
                snt = text(2.5e3,-12.3,'$!8S/N!7 \geq  5.0$',target=p,/DATA,font_size=12,font_name='Times')
                snt = text(2.5e3,-13.1,'$!8S/N!7 \geq  3.0$',target=p,/DATA,font_size=12,font_name='Times',col='red')
            endif
        endif        
        
        ;; ensure axes range
        p.xra = xra[*,i]
        
        re = execute('pdet = plot(EXP'+xfield[i]+'[where(iiqual_det'+xfield[i]+')],alog10(FX'+xfield[i]+'[where(iiqual_det'+xfield[i]+')]),"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov,name="X-ray detected [$!8F!7_{X,observed}$]")')
        re = execute('pdo = plot(EXP'+xfield[i]+'[where(iiqual_det'+xfield[i]+')],alog10(FX'+xfield[i]+'[where(iiqual_det'+xfield[i]+')]),"S",sym_filled=0,sym_size=1.5,transparency=85,/ov)')
        if (i eq nfield-1) then begin
            re = execute('pright = plot(TEXP'+xfield[i]+'[where(iiqual_non'+xfield[i]+')],logfxir[where(iiqual_non'+xfield[i]+')],/nodata,/current,position=src_pos[*,i],axis_style=0)')
            aright = axis('y',target=pnon,location=[pnon.xra[1],0,0],textpos=1,tickdir=1,title='$log  !8F!7_{2-10 keV}(expected) [erg s^{-1}cm^{-2}]$',tickfont_name='Times',tickfont_size=14)
        endif
        
        if (i eq leg_src) then l = legend(target=[pdet,pnon,pcat_dark],position=leg_pos,/normal,/auto_text_color,sample_width=0,horizontal_spacing=0.06,vertical_alignment=0.5,horizontal_alignment=0)
    endfor
    
    ;xt = text(xpos[0],xpos[1],'$exposure time [s]$',alignment=0.5,font_size=14,font_name='Times')
    ;yt = text(ypos[0],ypos[1],'$log  !8F!7_{2-10 keV} [erg s^{-1}cm^{-2}]$',orientation=90,alignment=0.5,font_size=14,font_name='Times')    
    
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/fx_limits.eps',/BITMAP else $
                                                     p.save,'figures/fx_limits.png',resolution=res
    endif
endif



;;----------------------------------------------------------------------------------------
;; SEDs
;;----------------------------------------------------------------------------------------
if keyword_set(seds) then begin
    print, '    PREPARING PAPER SED'
    igal = ['1237652900229415056','1237655370892116222','1237651272441725119']
    ids = ['1237679254133735725','1237662663208206976','1237679437739524714','1237679321788317973',igal[0:1]]
    ;; others     '1237664291009921910'
    nplt = n_elements(ids)
    inds = []
    for i = 0,nplt-1 do inds = [inds,where(objid eq ids[i])]
    print, 'DETECT:  ', iiqual_det[inds]
    print, 'NON-DET: ', iiqual_non[inds]
    print, 'WISE AGN:', iidet_wac[inds]
    ;; determine which template components
    components = tag_names(comp)
    ;; all possible templates (SED modeling procedure can handle max=5 templates)
    temps = ['AGN','ELL','SFG','IRR','DST']   
    ;; colors for plotting
    col = [[204,121,167],[213,94,0],[0,158,115],[0,114,178],[240,228,66]]
    linestyle = ['-','--','__','-.','__']
    ;col = ['magenta','magenta','dark orange','sea green','sea green', 'blue','brown']
    ;; match input components (use MATCH2.PRO to keep named order of TEMPS; MATCH.PRO alphabetizes; important for plotting purposes)
    match2,components,temps,icomp,itemp
    ;; ensure we contain at least one valid template and sort
    if (total(itemp ne -1) le 0) then stop		           
    temps = temps[where(itemp ne -1)]
    col = col[*,where(itemp ne -1)]
    linestyle = linestyle[where(itemp ne -1)]
    ntemps = n_elements(temps)
    ;; extract indices of sources to plot
    plot_flux = flux[*,inds]
    plot_err = e_flux[*,inds]
    plot_bin = bin[*,inds]
    plot_ra = ra[inds]
    plot_dec = dec[inds]
    plot_id = objid[inds]
    plot_z = z[inds]
    plot_zerr = zerr[inds]
    plot_fits = param[*,inds]
    plot_ebvsig = sig_ebv[1,inds]>0.
    ;; extract model parameters
    plot_ebv = plot_fits[0,*]
    coeff = plot_fits[2:2+ntemps-1,*]
    plot_chi = plot_fits[-2:-1,*]
    obswav = wave
    objwav = rebin(obswav,n_elements(obswav),nplt)
    objnu = (!const.c*1e6)/objwav
    tempwav = comp.wav#reform(1+plot_z)
    tempnu = (!const.c*1e6)/tempwav
    ;; covert data from plot_flux density [microjansky] to plot_flux [erg/s/cm2]
    plot_err *= 1e-29 * objnu         
    plot_flux *= 1e-29 * objnu
    ;; reconstruct models
    ;; convert models from plot_flux density [microjansky] to plot_flux [erg/s/cm2]
    agn = 1e-29 * tempnu * (coeff[0,*]##comp.(where(strmatch(tag_names(comp),'AGN*')))) * 10.^(-0.4 * comp.kap # plot_ebv)                         ;; AGN model
    for i = 1,ntemps-1 do re = execute(temps[i]+' = 1e-29 * tempnu * (coeff[i,*]##comp.'+temps[i]+')')  ;; galaxy models
    re = execute('model = '+strjoin(temps,"+"))                                                         ;; coadded models
    ;; convert to log scale
    plot_err = abs((plot_err)/(plot_flux*alog(10)))
    plot_flux = alog10(plot_flux)
    for i = 0,ntemps-1 do re = execute(temps[i]+' = alog10('+temps[i]+')')
    model = alog10(model)
    ;print, agn[value_locate(tempwav[*,0],6.),0]
    ;; string variables for call to TEXT()
    ;plot_pos = strtrim(plot_ra,2)+','
    ;plot_pos[where(plot_dec gt 0.)] += '+'
    ;plot_pos[*] += strtrim(plot_dec,2)
    plot_ra = strtrim(plot_ra,2)
    str_dec = strarr(nplt)+'+'
    str_dec[where(plot_dec gt 0.)] += strtrim(plot_dec[where(plot_dec gt 0.)],2)
    str_dec[where(plot_dec lt 0.)] = strtrim(plot_dec[where(plot_dec lt 0.)],2)
    plot_dec = str_dec
    plot_z = strtrim(string(plot_z,format='(d5.3)'),2)
    plot_zerr = strtrim(string(plot_zerr,format='(d5.3)'),2)
    plot_ebv = strtrim(string(plot_ebv,format='(d5.2)'),2)
    plot_ebvsig = strtrim(string(plot_ebvsig,format='(d4.2)'),2)
    coeff = reform(strtrim(string(coeff,format='(e10.3)'),2),ntemps,nplt)
    plot_chi = strtrim(string(plot_chi[0,*],format='(d0.2)'),2)+' / '+strtrim(string(plot_chi[1,*],format='(i)'),2)
    ;; plot SEDs
    e = {xr:[0.05,30.],yra:[-15.5,-8.5], $     ;[floor(min(plot_flux[where(finite(plot_flux))]))-1.5,ceil(max(plot_flux[where(finite(plot_flux))]))+2.0], $
         xlog:1, $
         font_name:'Times',font_size:14,$
         ;xtitle:'$\lambda (observed) [ \mum ]$',ytitle:'$log  \nu!8F!7_\nu  [erg s^{-1}cm^{-2}]$', $
         ;nodata:1,dimension:[1140,890], $       pos = [[80,455,585,835],[585,455,1090,835],[80,75,585,455],[585,75,1090,455]]
         nodata:1,dimension:[1140,1270], $      
         buffer:0}
    if keyword_set(hide) then e.buffer = 1
    pos = [[80,835,585,1215],[585,835,1090,1215],[80,455,585,835],[585,455,1090,835],[80,75,585,455],[585,75,1090,455]]
    ;label = transpose([['         ID : '+strtrim(plot_id,2)],['$            !8z!7 : $'+plot_z],['$!8E!7(!8B-V!7)$ : '+plot_ebv],['$\chi^2 / DoF$ : '+chi]])
    ;label = transpose([['ID :          '+strtrim(plot_id,2)],['$!8z!7 :             $'+plot_z],['$!8E!7(!8B-V!7)$ : '+plot_ebv],['$\chi^2 / DoF$ : '+chi]])
    ;label = transpose([['ObjID: '+strtrim(plot_id,2)],['$!8z!7 = $'+plot_z+'$\pm$'+plot_zerr],['$!8E!7(!8B-V!7)$ = '+plot_ebv+'$\pm$'+plot_ebvsig],['$\chi^2 / DoF$ = '+plot_chi]])
    ;label = transpose([['RA/Dec: '+plot_pos],['$!8z!7: $'+plot_z+'$\pm$'+plot_zerr],['$!8E!7(!8B-V!7)$: '+plot_ebv+'$\pm$'+plot_ebvsig],['$\chi^2 / DoF$: '+plot_chi]])
    label = transpose([['RA: '+plot_ra],['Dec: '+plot_dec],['$!8z!7: $'+plot_z+'$\pm$'+plot_zerr],['$!8E!7(!8B-V!7)$: '+plot_ebv+'$\pm$'+plot_ebvsig],['$\chi^2 / DoF$: '+plot_chi]])
    for i = 0,nplt-1 do begin
        if (i eq 0) then current = 0 else current = 1
        ig = where(plot_bin[*,i],/null)
        p = plot(objwav[ig,i],plot_flux[ig,i],_extra=e,position=pos[*,i],/DEVICE,current=current)
        if (i lt 4) then p.axes[0].showtext = 0
        if (i mod 2) then p.axes[1].showtext = 0
        ;for t = 0,ntemps-1 do re = execute('p = plot(tempwav[*,i],'+temps[t]+'[*,i],col=col[*,t],thick=2,linestyle=linestyle[t],/ov)')   ;; plot models
        p_agn = plot(tempwav[*,i],agn[*,i],col=col[*,0],thick=2,linestyle=linestyle[0],/ov,name=' AGN')
        p_ell = plot(tempwav[*,i],ell[*,i],col=col[*,1],thick=2,linestyle=linestyle[1],/ov,name=' ELL')
        p_sfg = plot(tempwav[*,i],sfg[*,i],col=col[*,2],thick=2,linestyle=linestyle[2],/ov,name=' SFG')
        p_irr = plot(tempwav[*,i],irr[*,i],col=col[*,3],thick=2,linestyle=linestyle[3],/ov,name=' IRR')
        p = plot(tempwav[*,i],model[*,i],col='dark slate grey',thick=2,/ov)                                                           ;; plot coadded models
        p = errorplot(objwav[ig,i],plot_flux[ig,i],plot_err[ig,i],'o',SYM_FILLED=1,LINESTYLE='',sym_size=1.5,errorbar_thick=2,/OV)          ;; plot photometry
        for t = 0,n_elements(label[*,0])-1 do begin
            lab = text(0.08,e.yra[1]-0.2-(1.5+t)/2.4,label[t,i],font_size=16,/DATA,target=p,font_name='Times')
            ;fit = text(1.95,yp-t*0.35,temps[t]+': '+coeff[t,i],col=col[*,t],font_size=12,/DATA,TARGET=p,font_name='Times')
        endfor
        if (i eq 1) then l = legend(target=[p_agn,p_ell,p_sfg,p_irr],position=[0.91,0.94],/auto_text_color,sample_width=0.07,horizontal_spacing=0.06,font_name='Times')
    endfor
    xt = text(0.52,0.025,'$\lambda (observed) [ \mum ]$',alignment=0.5,font_size=16,font_name='Times')
    yt = text(0.028,0.52,'$log  \nu!8F!7_\nu  [erg s^{-1}cm^{-2}]$',orientation=90,alignment=0.5,font_size=14,font_name='Times')
    
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/sed_models.eps',/BITMAP else $
                                                     p.save,'figures/sed_models.png',resolution=res
    endif
endif



;;----------------------------------------------------------------------------------------
;; LX vs LIR
;;----------------------------------------------------------------------------------------
if keyword_set(lx_lir) then begin
    print, '    PREPARING PAPER LX vs LIR'
    
    file = file_search('chen17_bw?.png')
    if (file[0] eq '') then begin
        spawn,'cp /Users/ccarroll/Research/projects/xray_lack_agn/workspace/data_prep/chen17_bw*.png .'
    endif
    
    ;; Chen+17 LX-LIR
    xrel_chen = [40.:50.:0.01]
    yrel_chen = dblarr(n_elements(xrel_chen))
    iichen = xrel_chen lt 44.79
    yrel_chen[where(iichen)] = 0.84*(xrel_chen[where(iichen)]-45.)+44.60
    yrel_chen[where(~iichen)] = 0.40*(xrel_chen[where(~iichen)]-45.)+44.51   
    ;; Fiore+09 LX-LIR
    xrel_fiore = [40.:50.:0.01]
    yrel_fiore = dblarr(n_elements(xrel_fiore))
    iifiore = xrel_fiore lt 43.04
    yrel_fiore[where(iifiore)] = xrel_fiore[where(iifiore)] - 0.3
    yrel_fiore[where(~iifiore)] = 43.574 + 0.72*(xrel_fiore[where(~iifiore)] - 44.2)
    
    ;lx = dblarr(nsrc)
    ;lx[where(lx eq 0. and iiqual_det_cha)] = lx_cha[where(lx eq 0. and iiqual_det_cha)]
    ;lx[where(lx eq 0. and iiqual_det_xmm)] = lx_xmm[where(lx eq 0. and iiqual_det_xmm)]
    ;lx[where(lx eq 0. and iiqual_det_nst)] = lx_nst[where(lx eq 0. and iiqual_det_nst)]

    iqso1 = where(iiqual_det and ebv le 0.15)
    iqso2 = where(iiqual_det and ebv gt 0.15)
    
    ;; test (1+z) correction
    ;loglirz = alog10((10.^loglir)/(1+z))
    ;loglirz[where(loglir eq -9999)] = -9999
    ;loglir = loglirz

    binsz = 0.15
    hist1 = hist_2d(loglir[iqso1],loglx[iqso1],bin1=binsz,min1=41.,max1=47.5,bin2=binsz,min2=41.,max2=46.5)
    hist2 = hist_2d(loglir[iqso2],loglx[iqso2],bin1=binsz,min1=41.,max1=47.5,bin2=binsz,min2=41.,max2=46.5)
    ind1 = array_indices(hist1,where(hist1 eq max(hist1)))
    ind2 = array_indices(hist2,where(hist2 eq max(hist2)))
    xlum = [41.:47.5:binsz]
    ylum = [41.:46.5:binsz]
    
    e = {xra:[41.,47.5],yra:[41.,46.5], $
         font_name:'Times',font_size:14, $
         xtitle:'$log  !8L!7_{6 \mu m} [erg s^{-1}cm^{-2}]$', $
         buffer:0}
    if keyword_set(hide) then e.buffer = 1
    rel_lin = ['-','-']
    ;rel_col = [[105,105,105],[0,0,0]]
    rel_col = [[0,158,115],[213,94,0]]
    ;; [0,158,115]  [204,121,167]
    im1 = image('chen17_bw1.png',transparency=50,dimension=[640,877],position=[80,440,587,810],buffer=buff,/device)
    pnodata = plot(xrel_chen,yrel_chen,_extra=e,/current,position=im1.position)
    cqso1 = contour(hist1,xlum,ylum,c_thick=4,rgb_table=colortable(49,/reverse),/fill,/ov,c_label_show=0,transparency=0,name='$!8E(B-V)!7 \leq  0.15$')
    prel_chen = plot(xrel_chen,yrel_chen,linestyle=rel_lin[0],col=rel_col[*,0],thick=4,_extra=e,/ov,position=im1.position,name='Chen+17')
    prel_fiore = plot(xrel_fiore,yrel_fiore,linestyle=rel_lin[1],col=rel_col[*,1],thick=4,/ov,name='Fiore+09')
    chent = text(41.15,42.3,/data,target=prel_chen,'Chen+17',col=rel_col[*,0],font_size=14,font_style='Bold',font_name='Times',fill_background=1,fill_color='white')
    fioret = text(41.95,41.3,/data,target=prel_fiore,'Fiore+09',col=rel_col[*,1],font_size=14,font_style='Bold',font_name='Times',fill_background=1,fill_color='white')
    prel_chen.axes[0].showtext = 0
    ;p = plot([median(xlum[ind1[0,*]])],[median(ylum[ind1[1,*]])],'o',col='black',sym_size=1.5,/ov)
    ;p = plot([median(xlum[ind1[0,*]])],[median(ylum[ind1[1,*]])],'S',col='black',sym_size=1.0,sym_filled=1,/ov)
    ;p = plot(xlum[ind1[0,*]],ylum[ind1[1,*]],'X',sym_size=1.5,sym_thick=1,/sym_filled,col=[0,158,115],/ov)
    im2 = image('chen17_bw2.png',transparency=50,/current,position=[80,75,587,445],/device)
    pnodata = plot(xrel_chen,yrel_chen,_extra=e,/current,position=im2.position)    
    cqso2 = contour(hist2,xlum,ylum,c_thick=4,rgb_table=colortable(62,/reverse),/fill,/ov,c_label_show=0,transparency=0,name='$!8E(B-V)!7 > 0.15$')
    prel_chen = plot(xrel_chen,yrel_chen,'-',col=rel_col[*,0],thick=4,_extra=e,/ov,position=im2.position,name='Chen+17')    
    prel_fiore = plot(xrel_fiore,yrel_fiore,'-',col=rel_col[*,1],thick=4,/ov,name='Fiore+09')
    ;p = plot([median(xlum[ind2[0,*]])],[median(ylum[ind2[1,*]])],'o',col='black',sym_size=1.5,/ov)
    ;p = plot([median(xlum[ind2[0,*]]),0,0,0],[median(ylum[ind2[1,*]]),0,0,0],'S',col='black',sym_size=1.0,sym_filled=1,_extra=e,/ov)
    ;p = plot(xlum[ind2[0,*]],ylum[ind2[1,*]],'X',sym_size=1.5,sym_thick=3,/sym_filled,col=[0,158,115],/ov)
    leg = legend(target=[cqso1,cqso2],position=[0.15,0.41],/normal,vertical_alignment=0.,horizontal_alignment=0.,font_name='Times')   
    ;leg.position = [0.14,0.39]

    ;xt = text(0.52,0.03,'$log  !8L!7_{6 \mu m} [erg s^{-1}cm^{-2}]$',alignment=0.5,font_size=14,font_name='Times')
    yt = text(0.05,0.51,'$log  !8L!7_{2-10 keV} [erg s^{-1}cm^{-2}]$',orientation=90,alignment=0.5,font_size=14,font_name='Times')
        
    ;peak_circ = text(0.72,0.13,'$\U(25EF)$',col='dark slate grey',font_size=12,font_name='Times')
    ;peak_star = text(0.721,0.128,'$\U(2605)$',col='dark slate grey',font_size=14,font_name='Times')
    ;peak_text = text(0.76,0.128,'Peak Freq.',col='dark slate grey',font_size=12,font_name='Times')
    
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then im1.save,'figures/lx_v_lir.eps',/BITMAP else $
                                                     im1.save,'figures/lx_v_lir.png',resolution=res
    endif
endif




;;----------------------------------------------------------------------------------------
;; Luminosity Ratio
;;----------------------------------------------------------------------------------------
if keyword_set(lum_ratio) then begin
    print, '    PREPARING LUMINOSITY RATIOS'
    
    ;iiqual_non = iiqual_non and (ztype ne 'ZP' or (photoerrorclass eq 1)); and photoerrorclass le 3))
    ;iiqual_det = iiqual_det and (ztype ne 'ZP' or (photoerrorclass eq 1)); and photoerrorclass le 3))
    
    model = 'POWER'
    ;model = 'BORUS'
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; PLOT L/L VS E(B-V)
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ell = {xra:[-3.8,2.2],yra:[-3.,1.], $
           font_name:'Times',font_size:14, $
           dimension:[1180,880], $
           buffer:0}
    if keyword_set(hide) then ell.buffer = 1

    nh = [21.:25.5:0.25]
    ll = rl2nh(nh,model=model,/lum_out)
    nh_values = alog10([1e23,5e23,1e24,1.5e24,3e24,1e25])
    nhtext = ['$1\times10^{23}$','$5\times10^{23}$','$1\times10^{24}$','$1.5\times10^{24}$','$3\times10^{24}$','$1\times10^{25}$']

    nh_lines = [ell.xra[0],ell.xra[0]+0.9]
    ll_lines = interpol(ll,nh,nh_values)
    ypos = ll_lines+0.08
    
    ;; [unobscured,Compton-thick]     
    nh_bound = alog10([1e21,1.5e24])
    ll_bound = rl2nh(nh_bound,model=model,/lum_out)

    ;; WISE AGN Catalog sources
    pw = plot(logebv[where(iiqual_non)],llnon[where(iiqual_non)],_extra=ell,position=[85,445,595,825],/device,/nodata,ytickvalues=[-3.0,-2.0,-1.0,0.0,1.0],ytickformat='(d4.1)',ytitle='$log  !8L!7_X / !8L!7_X(!8L!7_{IR})  (2$-$10 keV)$')
    ;; unobscured shading
    pw = plot(ell.xra,[1,1]*ell.yra[1],linestyle='',fill_background=1,fill_level=0.,fill_color='light grey',fill_transparency=0,/ov)
    ;; CT shading
    pw = plot(ell.xra,ll_bound[1]*[1,1],linestyle='',fill_background=1,fill_level=ell.yra[0],fill_color='light grey',fill_transparency=00,/ov)
    ;; plot data
    pwllnon = arrow(transpose([[logebv[where(iiqual_non and iidet_wac)]],[logebv[where(iiqual_non and iidet_wac)]]]),transpose([[llnon[where(iiqual_non and iidet_wac)]],[llnon[where(iiqual_non and iidet_wac)]-0.12]]),color="orange",fill_transparency=75,target=[pw],/data,head_size=0.3,clip=1)
    ;pwllnon = plot(logebv[where(iiqual_non and iidet_wac)],llnon[where(iiqual_non and iidet_wac)],"td",sym_filled=1,sym_size=1.,color="orange",transparency=50,/ov,_extra=e,name=" X-ray non-det.",/nodata)    
    pwlldet = plot(logebv[where(iiqual_det and iidet_wac)],lldet[where(iiqual_det and iidet_wac)],"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov,name=" X-ray detected")
    p = plot(logebv[where(iiqual_det and iidet_wac)],lldet[where(iiqual_det and iidet_wac)],"S",sym_filled=0,sym_size=1.5,transparency=85,/ov)
    ;; non-detections off plot range
    ;ioff = where(iiqual_non and iidet_wac and llnon ge ell.yra[1],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[1]-0.12)],[make_array(noff,value=ell.yra[1])-0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;ioff = where(iiqual_non and iidet_wac and llnon le ell.yra[0],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[0]+0.12)],[make_array(noff,value=ell.yra[0])+0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;;; detections off plot range
    ;ioff = where(iiqual_det and iidet_wac and lldet gt ell.yra[1],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[1]-0.12)],[make_array(noff,value=ell.yra[1])-0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;ioff = where(iiqual_det and iidet_wac and lldet lt ell.yra[0],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[0]+0.12)],[make_array(noff,value=ell.yra[0])+0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; median error bars
    perr = plot([-2.2],[ell.yra[0]-(ell.yra[0]-ll_bound[1])/2.],'s',sym_size=4,sym_thick=2,/ov)
    perr = errorplot([-2.2],[ell.yra[0]-(ell.yra[0]-ll_bound[1])/2.],[median(e_lldet[where(iiqual_det and iidet_wac)])],errorbar_capsize=0.1,linestyle="",/ov)
    perr = plot([-2.2],[ell.yra[0]-(ell.yra[0]-ll_bound[1])/2.],"S",sym_filled=1,sym_size=1.,col='white',/ov)
    perr = plot([-2.2],[ell.yra[0]-(ell.yra[0]-ll_bound[1])/2.] ,"S",sym_filled=0,sym_size=1.,col='dodger blue',/ov)
    ;; CT lines
    p = plot(ell.xra[1]+[-0.7,0.],ll_bound[0]*[1,1],'-',thick=2,/ov)
    p = plot(ell.xra[1]+[-0.7,0.],ll_bound[1]*[1,1],'-',thick=2,/ov)
    ;; catalog sources text
    t = text(ell.xra[0]+0.22,0.50,'!16WISE!15 AGN',col='red',/data,target=[pwllnon],font_size=12,font_name='Times',font_style='Bold')
    t = text(ell.xra[0]+0.22,0.25,'sources',col='red',/data,target=[pwllnon],font_size=12,font_name='Times',font_style='Bold')
    t = text(ell.xra[1]-0.2,ll_bound[0]+0.14,'Unobscured',col='black',target=pwllnon,alignment=1.,font_size=12,font_name='Times',font_style='Bold',/data)
    t = text(ell.xra[1]-0.2,ll_bound[1]-0.25,'Compton',col='black',target=pwllnon,alignment=1.,font_size=12,font_name='Times',font_style='Bold',/data)
    t = text(ell.xra[1]-0.2,ll_bound[1]-0.45,'thick',col='black',target=pwllnon,alignment=1.,font_size=12,font_name='Times',font_style='Bold',/data)
    
    ;; NH lines + text
    t = text(ell.xra[0]+0.2,ell.yra[0]+0.4,'$!8N!7_H  [cm^{-2}]$',col='black',/data,target=[pwllnon],font_size=11,font_name='Times')
    for i = 0,n_elements(ll_lines)-1 do begin
        p = plot(nh_lines,ll_lines[i]*[1.,1.],'--',col='black',thick=2,/ov)
        t = text(ell.xra[0]+0.2,ypos[i],nhtext[i],col='black',/data,target=[pwllnon],font_size=11,font_name='Times')
    endfor
    ;; remove axes
    pwlldet.axes[0].showtext=0
    ;; model choice
    mt = text(ell.xra[1]-diff(ell.xra)/2.,ell.yra[0]+0.5,'Power law model',target=pwlldet,alignment=0.5,vertical_alignment=0.5,font_size=14,font_name='Times',font_style='Bold',/data)

    ;; Remaining sources
    pr = plot(logebv[where(iiqual_non)],llnon[where(iiqual_non)],_extra=ell,position=[595,445,1105,825],current=1,/device,/nodata)
    pr.axes[1].showtext = 0 & pr.axes[0].showtext = 0
    ;; unobscured shading
    pr = plot(ell.xra,[1,1]*ell.yra[1],linestyle='',fill_background=1,fill_level=0.,fill_color='light grey',fill_transparency=0,/ov)
    ;; CT shading
    pr = plot(ell.xra,ll_bound[1]*[1,1],linestyle='',fill_background=1,fill_level=ell.yra[0],fill_color='light grey',fill_transparency=0,/ov)
    ;; plot data
    prllnon = arrow(transpose([[logebv[where(iiqual_non and ~iidet_wac)]],[logebv[where(iiqual_non and ~iidet_wac)]]]),transpose([[llnon[where(iiqual_non and ~iidet_wac)]],[llnon[where(iiqual_non and ~iidet_wac)]-0.12]]),color="orange",fill_transparency=75,target=[pr],/data,head_size=0.3,clip=1)
    prlldet = plot(logebv[where(iiqual_det and ~iidet_wac)],lldet[where(iiqual_det and ~iidet_wac)],"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov,name=" X-ray detected")
    p = plot(logebv[where(iiqual_det and ~iidet_wac)],lldet[where(iiqual_det and ~iidet_wac)],"S",sym_filled=0,sym_size=1.5,transparency=85,/ov)
    prllnon = plot(logebv[where(iiqual_non and ~iidet_wac)],llnon[where(iiqual_non and ~iidet_wac)],"td",sym_filled=1,sym_size=1.,color="orange",/current,/nodata,name=" X-ray non-det.",position=[595,445,1105,825],_extra=ell,axis_style=0)
    a_rllnon = axis('y',target=pr,location=[ell.xra[1],0,0],textpos=1,tickdir=1,title='$log  !8L!7_X(!8F!7_X lim) / !8L!7_X(!8L!7_{IR})  (2$-$10 keV)$',tickformat='(d4.1)',tickfont_name='Times',tickfont_size=14)
    ;; non-detections off plot range
    ;ioff = where(iiqual_non and ~iidet_wac and llnon ge ell.yra[1],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[1]-0.12)],[make_array(noff,value=ell.yra[1])-0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;ioff = where(iiqual_non and ~iidet_wac and llnon le ell.yra[0],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[0]+0.12)],[make_array(noff,value=ell.yra[0])+0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;;; detections off plot range
    ;ioff = where(iiqual_det and ~iidet_wac and lldet gt ell.yra[1],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[1]-0.12)],[make_array(noff,value=ell.yra[1])-0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;ioff = where(iiqual_det and ~iidet_wac and lldet lt ell.yra[0],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[0]+0.12)],[make_array(noff,value=ell.yra[0])+0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; CT lines
    p = plot(ell.xra[1]+[-0.7,0.],ll_bound[0]*[1,1],'-',thick=2,/ov)
    p = plot(ell.xra[1]+[-0.7,0.],ll_bound[1]*[1,1],'-',thick=2,/ov)
    ;; Catalog sources text
    t = text(ell.xra[0]+0.22,0.50,'Remaining',col='red',/data,target=[prlldet],font_size=12,font_name='Times',font_style='Bold')
    t = text(ell.xra[0]+0.22,0.25,'sources',col='red',/data,target=[prlldet],font_size=12,font_name='Times',font_style='Bold')
    ;; NH lines + text
    t = text(ell.xra[0]+0.2,ell.yra[0]+0.4,'$!8N!7_H  [cm^{-2}]$',col='black',/data,target=[prllnon],font_size=11,font_name='Times')
    for i = 0,n_elements(ll_lines)-1 do begin
        p = plot(nh_lines,ll_lines[i]*[1.,1.],'--',col='black',thick=2,/ov)
        t = text(ell.xra[0]+0.2,ypos[i],nhtext[i],col='black',/data,target=[prllnon],alignment=0.,font_size=11,font_name='Times')
    endfor
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; PLOT NH VS E(B-V)
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    enh = {xra:[-3.8,2.2],yra:[20.5,25.5], $
           font_name:'Times',font_size:14, $
           dimension:[1130,880], $
           buffer:0}
    if keyword_set(hide) then enh.buffer = 1
    
    nhxdet = rl2nh(lldet,model=model)
    nhxnon = rl2nh(llnon,model=model)
    
    ;; [unobscured,Compton-thick]     
    nh_bound = alog10([1e21,1.5e24])

    ;; WISE AGN Catalog sources
    pw = plot(logebv[where(iiqual_non and iidet_wac)],llnon[where(iiqual_non and iidet_wac)],_extra=enh,position=[85,65,595,445],current=1,/device,/nodata,ytitle='$log  !8N!7_H [cm^{-2}]$',ytickvalues=[21.,22.,23.,24.,25.],ytickformat='(i2)')
    ;; CT shading
    pw = plot(enh.xra,nh_bound[1]*[1,1],linestyle='',fill_background=1,fill_level=enh.yra[1],fill_color='light grey',fill_transparency=0,/ov)
    ;; plot data
    pwnhnon = arrow(transpose([[logebv[where(iiqual_non and iidet_wac)]],[logebv[where(iiqual_non and iidet_wac)]]]),transpose([[nhxnon[where(iiqual_non and iidet_wac)]],[nhxnon[where(iiqual_non and iidet_wac)]+0.12]]),color="orange",fill_transparency=75,target=[pw],/data,head_size=0.3,clip=1)
    ;pwnhnon = plot(logebv[where(iiqual_non and iidet_wac)],nhxnon[where(iiqual_non and iidet_wac)],"o",sym_filled=1,sym_size=0.5,color="orange",transparency=50,/ov,name=" X-ray non-det.")
    pwnhdet = plot(logebv[where(iiqual_det and iidet_wac)],nhxdet[where(iiqual_det and iidet_wac)],"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov,name=" X-ray detected")
    p = plot(logebv[where(iiqual_det and iidet_wac)],nhxdet[where(iiqual_det and iidet_wac)],"S",sym_filled=0,sym_size=1.5,transparency=85,/ov)
    ;; non-detections off plot range
    ;ioff = where(iiqual_non and iidet_wac and nhxnon ge enh.yra[1],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[1]-0.12)],[make_array(noff,value=enh.yra[1])-0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;ioff = where(iiqual_non and iidet_wac and nhxnon le enh.yra[0],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[0]+0.12)],[make_array(noff,value=enh.yra[0])+0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;;; detections off plot range-0.2
    ;ioff = where(iiqual_det and iidet_wac and nhxdet gt enh.yra[1],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[1]-0.12)],[make_array(noff,value=enh.yra[1])-0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;ioff = where(iiqual_det and iidet_wac and nhxdet lt enh.yra[0],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[0]+0.12)],[make_array(noff,value=enh.yra[0])+0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; median error bars
    perr = plot([-2.2],[enh.yra[1]-(enh.yra[1]-nh_bound[1])/2.],'s',sym_size=4,sym_thick=2,/ov)
    perr = errorplot([-2.2],[enh.yra[1]-(enh.yra[1]-nh_bound[1])/2.],[mean(resamp_rlnh(ret='nh_mad'))],errorbar_capsize=0.1,linestyle="",/ov)
    perr = plot([-2.2],[enh.yra[1]-(enh.yra[1]-nh_bound[1])/2.],"S",sym_filled=1,sym_size=1.,col='white',/ov)
    perr = plot([-2.2],[enh.yra[1]-(enh.yra[1]-nh_bound[1])/2.],"S",sym_filled=0,sym_size=1.,col='dodger blue',/ov)
    ;; CT lines
    p = plot(enh.xra[0]+[0.,0.7],nh_bound[1]*[1,1],'-',thick=2,/ov)
    ;; catalog sources text
    t = text(enh.xra[0]+0.22,nh_bound[1]+0.40,'Compton',col='black',target=pwnhnon,font_size=12,font_name='Times',font_style='Bold',/data)
    t = text(enh.xra[0]+0.22,nh_bound[1]+0.15,'thick',col='black',target=pwnhnon,font_size=12,font_name='Times',font_style='Bold',/data)

    ;; Remaining sources
    pr = plot(logebv[where(iiqual_non and ~iidet_wac)],nhxnon[where(iiqual_non and ~iidet_wac)],_extra=enh,position=[595,65,1105,445],current=1,/device,/nodata,ytickvalues=[21.,22.,23.,24.,25.],ytickformat='(i2)')
    pr.axes[1].showtext=0
    ;; CT shading
    pr = plot(enh.xra,nh_bound[1]*[1,1],linestyle='',fill_background=1,fill_level=enh.yra[1],fill_color='light grey',fill_transparency=0,/ov)
    ;; plot data
    prnhnon = arrow(transpose([[logebv[where(iiqual_non and ~iidet_wac)]],[logebv[where(iiqual_non and ~iidet_wac)]]]),transpose([[nhxnon[where(iiqual_non and ~iidet_wac)]],[nhxnon[where(iiqual_non and ~iidet_wac)]+0.12]]),color="orange",fill_transparency=75,target=[pr],/data,head_size=0.3,clip=1)
    prnhdet = plot(logebv[where(iiqual_det and ~iidet_wac)],nhxdet[where(iiqual_det and ~iidet_wac)],"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov,name=" X-ray detected")
    p = plot(logebv[where(iiqual_det and ~iidet_wac)],nhxdet[where(iiqual_det and ~iidet_wac)],"S",sym_filled=1,sym_size=1.5,transparency=85,/ov)
    prnhnon = plot(logebv[where(iiqual_non and ~iidet_wac)],nhxnon[where(iiqual_non and ~iidet_wac)],"td",sym_filled=1,sym_size=1.,color="orange",/current,/nodata,position=[595,65,1105,445],_extra=ell,axis_style=0)
    a_rnhnon = axis('y',target=pr,location=[enh.xra[1],0,0],textpos=1,tickdir=1,title='$log  !8N!7_H (lim)[cm^{-2}]$',tickfont_name='Times',tickfont_size=14)    
    ;; non-detections off plot range
    ;ioff = where(iiqual_non and ~iidet_wac and nhxnon ge enh.yra[1],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[1]-0.12)],[make_array(noff,value=enh.yra[1])-0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;ioff = where(iiqual_non and ~iidet_wac and nhxnon le enh.yra[0],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[0]+0.12)],[make_array(noff,value=enh.yra[0])+0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;;; detections off plot range-0.2
    ;ioff = where(iiqual_det and ~iidet_wac and nhxdet gt enh.yra[1],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[1]-0.12)],[make_array(noff,value=enh.yra[1])-0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;ioff = where(iiqual_det and ~iidet_wac and nhxdet lt enh.yra[0],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[0]+0.12)],[make_array(noff,value=enh.yra[0])+0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; CT lines
    p = plot(enh.xra[0]+[0.,0.7],nh_bound[1]*[1,1],'-',thick=2,/ov)    
    ;; legend
    l = legend(target=[prlldet,prllnon],/normal,/auto_text_color,sample_width=0.,horizontal_spacing=0.06,font_name='Times')
    l.position = [0.67,0.49]

    xt = text(0.52,0.02,'$log  !8E!7(!8B-V!7)$',alignment=0.5,font_size=14,font_name='Times')
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/rlum_b.eps',/BITMAP else $
                                                     p.save,'figures/rlum_b.png',resolution=res
    endif

    ;model = 'POWER'
    model = 'BORUS'
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; ;;;;;;;;;;;;;;;;;;;;;;
    ;; PLOT L/L VS E(B-V)
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    nh = [18.0000,18.2500,18.5000,18.7500,19.0000,19.2500,19.5000,19.7500,20.0000,20.2500,20.5000,20.7500,21.0000,21.2500,21.5000,21.7500,22.0000,22.2500,22.5000,22.7500,23.0000,23.2500,23.5000,23.7500,24.0000,24.2500,24.5000,24.7500,25.0000,25.2500,25.5000]
    ll = rl2nh(nh,model=model,/lum_out)
    nh_values = alog10([1e23,5e23,1.5e24,1e25])
    nhtext = ['$1\times10^{23}$','$5\times10^{23}$','$1.5\times10^{24}$','$1\times10^{25}$']

    nh_lines = [ell.xra[0],ell.xra[0]+0.9]
    ll_lines = interpol(ll,nh,nh_values)
    ypos = ll_lines+0.08
    
    ;; [unobscured,Compton-thick]     
    nh_bound = alog10([1e21,1.5e24])
    ll_bound = rl2nh(nh_bound,model=model,/lum_out)

    ;; WISE AGN Catalog sources
    pw = plot(logebv[where(iiqual_non)],llnon[where(iiqual_non)],_extra=ell,position=[85,445,595,825],/device,/nodata,ytickvalues=[-3.0,-2.0,-1.0,0.0,1.0],ytickformat='(d4.1)',ytitle='$log  !8L!7_X / !8L!7_X(!8L!7_{IR})  (2$-$10 keV)$')
    ;; unobscured shading
    pw = plot(ell.xra,[1,1],linestyle='',fill_background=1,fill_level=0.,fill_color='light grey',fill_transparency=0,/ov)
    ;; CT shading
    pw = plot(ell.xra,ll_bound[1]*[1,1],linestyle='',fill_background=1,fill_level=ell.yra[0],fill_color='light grey',fill_transparency=0,/ov)
    ;; plot data
    pwllnon = arrow(transpose([[logebv[where(iiqual_non and iidet_wac)]],[logebv[where(iiqual_non and iidet_wac)]]]),transpose([[llnon[where(iiqual_non and iidet_wac)]],[llnon[where(iiqual_non and iidet_wac)]-0.12]]),color="orange",fill_transparency=75,target=[pw],/data,head_size=0.3,clip=1)
    ;pwllnon = plot(logebv[where(iiqual_non and iidet_wac)],llnon[where(iiqual_non and iidet_wac)],"o",sym_filled=1,sym_size=0.5,color="orange",transparency=50,/ov,name=" X-ray non-det.")
    pwlldet = plot(logebv[where(iiqual_det and iidet_wac)],lldet[where(iiqual_det and iidet_wac)],"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov)
    p = plot(logebv[where(iiqual_det and iidet_wac)],lldet[where(iiqual_det and iidet_wac)],"S",sym_filled=0,sym_size=1.5,transparency=85,/ov)
    ;; non-detections off plot range
    ;ioff = where(iiqual_non and iidet_wac and llnon ge ell.yra[1],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[1]-0.12)],[make_array(noff,value=ell.yra[1])-0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;ioff = where(iiqual_non and iidet_wac and llnon le ell.yra[0],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[0]+0.12)],[make_array(noff,value=ell.yra[0])+0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;;; detections off plot range
    ;ioff = where(iiqual_det and iidet_wac and lldet gt ell.yra[1],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[1]-0.12)],[make_array(noff,value=ell.yra[1])-0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;ioff = where(iiqual_det and iidet_wac and lldet lt ell.yra[0],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[0]+0.12)],[make_array(noff,value=ell.yra[0])+0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; median error bars
    perr = plot([-2.2],[ell.yra[0]-(ell.yra[0]-ll_bound[1])/2.],'s',sym_size=4,sym_thick=2,/ov)
    perr = errorplot([-2.2],[ell.yra[0]-(ell.yra[0]-ll_bound[1])/2.],[median(e_lldet[where(iiqual_det and iidet_wac)])],errorbar_capsize=0.1,linestyle="",/ov)
    perr = plot([-2.2],[ell.yra[0]-(ell.yra[0]-ll_bound[1])/2.],"S",sym_filled=1,sym_size=1.,col='white',/ov)
    perr = plot([-2.2],[ell.yra[0]-(ell.yra[0]-ll_bound[1])/2.],"S",sym_filled=0,sym_size=1.,col='dodger blue',/ov)
    ;; CT lines
    p = plot(ell.xra[1]+[-0.7,0.],ll_bound[0]*[1,1],'-',thick=2,/ov)
    p = plot(ell.xra[1]+[-0.7,0.],ll_bound[1]*[1,1],'-',thick=2,/ov)
    ;; catalog sources text
    t = text(ell.xra[0]+0.22,0.50,'!16WISE!15 AGN',col='red',/data,target=[pwllnon],font_size=12,font_name='Times',font_style='Bold')
    t = text(ell.xra[0]+0.22,0.25,'sources',col='red',/data,target=[pwllnon],font_size=12,font_name='Times',font_style='Bold')
    t = text(ell.xra[1]-0.2,ll_bound[0]+0.14,'Unobscured',col='black',target=pwllnon,alignment=1.,font_size=12,font_name='Times',font_style='Bold',/data)
    t = text(ell.xra[1]-0.2,ll_bound[1]-0.25,'Compton',col='black',target=pwllnon,alignment=1.,font_size=12,font_name='Times',font_style='Bold',/data)
    t = text(ell.xra[1]-0.2,ll_bound[1]-0.45,'thick',col='black',target=pwllnon,alignment=1.,font_size=12,font_name='Times',font_style='Bold',/data)
    ;; NH lines + text
    t = text(ell.xra[0]+0.2,ell.yra[0]+0.4,'$!8N!7_H  [cm^{-2}]$',col='black',/data,target=[pwllnon],font_size=11,font_name='Times')
    for i = 0,n_elements(ll_lines)-1 do begin
        p = plot(nh_lines,ll_lines[i]*[1.,1.],'--',col='black',thick=2,/ov)
        t = text(ell.xra[0]+0.2,ypos[i],nhtext[i],col='black',/data,target=[pwllnon],font_size=11,font_name='Times')
    endfor
    ;; remove axes
    pwlldet.axes[0].showtext=0
    ;; model choice
    mt = text(ell.xra[1]-diff(ell.xra)/2.,ell.yra[0]+0.5,'BORUS model',target=pwlldet,alignment=0.5,vertical_alignment=0.5,font_size=14,font_name='Times',font_style='Bold',/data)

    ;; Remaining sources
    pr = plot(logebv[where(iiqual_non)],llnon[where(iiqual_non)],_extra=ell,position=[595,445,1105,825],current=1,/device,/nodata)
    ;; unobscured shading
    pr = plot(ell.xra,[1,1],linestyle='',fill_background=1,fill_level=0.,fill_color='light grey',fill_transparency=0,/ov)
    ;; CT shading
    pr = plot(ell.xra,ll_bound[1]*[1,1],linestyle='',fill_background=1,fill_level=ell.yra[0],fill_color='light grey',fill_transparency=0,/ov)
    ;; data
    prllnon = arrow(transpose([[logebv[where(iiqual_non and ~iidet_wac)]],[logebv[where(iiqual_non and ~iidet_wac)]]]),transpose([[llnon[where(iiqual_non and ~iidet_wac)]],[llnon[where(iiqual_non and ~iidet_wac)]-0.12]]),color="orange",fill_transparency=75,target=[pr],/data,head_size=0.3,clip=1)
    prlldet = plot(logebv[where(iiqual_det and ~iidet_wac)],lldet[where(iiqual_det and ~iidet_wac)],"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov,name=" X-ray detected")
    p = plot(logebv[where(iiqual_det and ~iidet_wac)],lldet[where(iiqual_det and ~iidet_wac)],"S",sym_filled=0,sym_size=1.5,transparency=85,/ov)
    prllnon = plot(logebv[where(iiqual_non and ~iidet_wac)],llnon[where(iiqual_non and ~iidet_wac)],"td",sym_filled=1,sym_size=1.,color="orange",/current,/nodata,name=" X-ray non-det.",position=[595,445,1105,825],_extra=ell,axis_style=0)
    a_rllnon = axis('y',target=pr,location=[ell.xra[1],0,0],textpos=1,tickdir=1,title='$log  !8L!7_X(!8F!7_X lim) / !8L!7_X(!8L!7_{IR})  (2$-$10 keV)$',tickformat='(d4.1)',tickfont_name='Times',tickfont_size=14)
    ;; non-detections off plot range
    ;ioff = where(iiqual_non and ~iidet_wac and llnon ge ell.yra[1],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[1]-0.12)],[make_array(noff,value=ell.yra[1])-0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;ioff = where(iiqual_non and ~iidet_wac and llnon le ell.yra[0],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[0]+0.12)],[make_array(noff,value=ell.yra[0])+0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;;; detections off plot range
    ;ioff = where(iiqual_det and ~iidet_wac and lldet gt ell.yra[1],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[1]-0.12)],[make_array(noff,value=ell.yra[1])-0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;ioff = where(iiqual_det and ~iidet_wac and lldet lt ell.yra[0],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[0]+0.12)],[make_array(noff,value=ell.yra[0])+0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; CT lines
    p = plot(ell.xra[1]+[-0.7,0.],ll_bound[0]*[1,1],'-',thick=2,/ov)
    p = plot(ell.xra[1]+[-0.7,0.],ll_bound[1]*[1,1],'-',thick=2,/ov)
    ;; Catalog sources text
    t = text(ell.xra[0]+0.22,0.50,'Remaining',col='red',/data,target=[prlldet],font_size=12,font_name='Times',font_style='Bold')
    t = text(ell.xra[0]+0.22,0.25,'sources',col='red',/data,target=[prlldet],font_size=12,font_name='Times',font_style='Bold')
    ;; NH lines + text
    t = text(ell.xra[0]+0.2,ell.yra[0]+0.4,'$!8N!7_H  [cm^{-2}]$',col='black',/data,target=[prllnon],font_size=11,font_name='Times')
    for i = 0,n_elements(ll_lines)-1 do begin
        p = plot(nh_lines,ll_lines[i]*[1.,1.],'--',col='black',thick=2,/ov)
        t = text(ell.xra[0]+0.2,ypos[i],nhtext[i],col='black',/data,target=[prllnon],alignment=0.,font_size=11,font_name='Times')
    endfor

    ;; remove axes
    prlldet.axes[1].showtext=0
    prlldet.axes[0].showtext=0
    
    ;xt = text(0.5,0.04,'$log  !8E!7(!8B-V!7)$',alignment=0.5,font_size=14,font_name='Times')
    ;yt = text(0.02,0.75,'$log  !8L!7_X / !8L!7_X(!8L!7_{IR})  (2$-$10 keV)$',orientation=90,alignment=0.5,font_size=14,font_name='Times')        
    ;yt.position = [0.007524362,0.59069969,0.026571130,0.85248213]

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; PLOT NH VS E(B-V)
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    nhxdet = rl2nh(lldet,model=model)
    nhxnon = rl2nh(llnon,model=model)
    
    ;; [unobscured,Compton-thick]     
    nh_bound = alog10([1e21,1.5e24])

    ;; WISE AGN Catalog sources
    pw = plot(logebv[where(iiqual_non and iidet_wac)],llnon[where(iiqual_non and iidet_wac)],_extra=enh,position=[85,65,595,445],current=1,/device,/nodata,ytitle='$log  !8N!7_H [cm^{-2}]$')
    ;; CT shading
    pw = plot(enh.xra,nh_bound[1]*[1,1],linestyle='',fill_background=1,fill_level=enh.yra[1],fill_color='light grey',fill_transparency=0,/ov)
    ;; plot data
    pwnhnon = arrow(transpose([[logebv[where(iiqual_non and iidet_wac)]],[logebv[where(iiqual_non and iidet_wac)]]]),transpose([[nhxnon[where(iiqual_non and iidet_wac)]],[nhxnon[where(iiqual_non and iidet_wac)]+0.12]]),color="orange",fill_transparency=75,target=[pw],/data,head_size=0.3,clip=1)
    ;pwnhnon = plot(logebv[where(iiqual_non and iidet_wac)],nhxnon[where(iiqual_non and iidet_wac)],"o",sym_filled=1,sym_size=0.5,color="orange",transparency=50,/ov,name=" X-ray non-det.")
    pwnhdet = plot(logebv[where(iiqual_det and iidet_wac)],nhxdet[where(iiqual_det and iidet_wac)],"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov)
    p = plot(logebv[where(iiqual_det and iidet_wac)],nhxdet[where(iiqual_det and iidet_wac)],"S",sym_filled=0,sym_size=1.5,transparency=85,/ov)
    ;; non-detections off plot range
    ;ioff = where(iiqual_non and iidet_wac and nhxnon ge enh.yra[1],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[1]-0.12)],[make_array(noff,value=enh.yra[1])-0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;ioff = where(iiqual_non and iidet_wac and nhxnon le enh.yra[0],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[0]+0.12)],[make_array(noff,value=enh.yra[0])+0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;;; detections off plot range
    ;ioff = where(iiqual_det and iidet_wac and nhxdet gt enh.yra[1],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[1]-0.12)],[make_array(noff,value=enh.yra[1])-0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;ioff = where(iiqual_det and iidet_wac and nhxdet lt enh.yra[0],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[0]+0.12)],[make_array(noff,value=enh.yra[0])+0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; median error bars
    perr = plot([-2.2],[enh.yra[1]-(enh.yra[1]-nh_bound[1])/2.],'s',sym_size=4,sym_thick=2,/ov)
    perr = errorplot([-2.2],[enh.yra[1]-(enh.yra[1]-nh_bound[1])/2.],[mean(resamp_rlnh(ret='nh_mad'))],errorbar_capsize=0.1,linestyle="",/ov)
    perr = plot([-2.2],[enh.yra[1]-(enh.yra[1]-nh_bound[1])/2.],"S",sym_filled=1,sym_size=1.,col='white',/ov)
    perr = plot([-2.2],[enh.yra[1]-(enh.yra[1]-nh_bound[1])/2.],"S",sym_filled=0,sym_size=1.,col='dodger blue',/ov)
    ;; CT lines
    p = plot(enh.xra[0]+[0.,0.7],nh_bound[1]*[1,1],'-',thick=2,/ov)
    ;; catalog sources text
    t = text(enh.xra[0]+0.22,nh_bound[1]+0.40,'Compton',col='black',target=pwnhnon,font_size=12,font_name='Times',font_style='Bold',/data)
    t = text(enh.xra[0]+0.22,nh_bound[1]+0.15,'thick',col='black',target=pwnhnon,font_size=12,font_name='Times',font_style='Bold',/data)

    ;; Remaining sources
    pr = plot(logebv[where(iiqual_non and ~iidet_wac)],nhxnon[where(iiqual_non and ~iidet_wac)],_extra=enh,position=[595,65,1105,445],current=1,/device,/nodata)
    ;; CT shading
    pr = plot(enh.xra,nh_bound[1]*[1,1],linestyle='',fill_background=1,fill_level=enh.yra[1],fill_color='light grey',fill_transparency=0,/ov)
    ;; plot data
    prnhnon = arrow(transpose([[logebv[where(iiqual_non and ~iidet_wac)]],[logebv[where(iiqual_non and ~iidet_wac)]]]),transpose([[nhxnon[where(iiqual_non and ~iidet_wac)]],[nhxnon[where(iiqual_non and ~iidet_wac)]+0.12]]),color="orange",fill_transparency=75,target=[pr],/data,head_size=0.3,clip=2)
    prnhdet = plot(logebv[where(iiqual_det and ~iidet_wac)],nhxdet[where(iiqual_det and ~iidet_wac)],"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov)
    p = plot(logebv[where(iiqual_det and ~iidet_wac)],nhxdet[where(iiqual_det and ~iidet_wac)],"S",sym_filled=1,sym_size=1.5,transparency=85,/ov)
    prnhnon = plot(logebv[where(iiqual_non and ~iidet_wac)],nhxnon[where(iiqual_non and ~iidet_wac)],"td",sym_filled=1,sym_size=1.,color="orange",/current,/nodata,position=[595,65,1105,445],_extra=ell,axis_style=0)
    a_rnhnon = axis('y',target=pr,location=[enh.xra[1],0,0],textpos=1,tickdir=1,title='$log  !8N!7_H (lim)[cm^{-2}]$',tickfont_name='Times',tickfont_size=14)
    ;; non-detections off plot range
    ;ioff = where(iiqual_non and ~iidet_wac and nhxnon ge enh.yra[1],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[1]-0.12)],[make_array(noff,value=enh.yra[1])-0.02]]),color='orange',/ov,thick=2,/data,fill_transparency=95,head_size=0.8,target=p)
    ;ioff = where(iiqual_non and ~iidet_wac and nhxnon le enh.yra[0],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[0]+0.12)],[make_array(noff,value=enh.yra[0])+0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;;; detections off plot range
    ;ioff = where(iiqual_det and ~iidet_wac and nhxdet gt enh.yra[1],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[1]-0.12)],[make_array(noff,value=enh.yra[1])-0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;ioff = where(iiqual_det and ~iidet_wac and nhxdet lt enh.yra[0],noff)
    ;if (noff gt 0) then pwoff = arrow(transpose([[logebv[ioff]],[logebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[0]+0.12)],[make_array(noff,value=enh.yra[0])+0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; CT lines
    p = plot(enh.xra[0]+[0.,0.7],nh_bound[1]*[1,1],'-',thick=2,/ov)
    ;; remove axes
    prnhdet.axes[1].showtext=0    
    ;; legend
    l = legend(target=[prlldet,prllnon],/normal,/auto_text_color,sample_width=0.,horizontal_spacing=0.06,font_name='Times')
    l.position = [0.67,0.49]
    
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/rlum_a.eps',/BITMAP else $
                                                     p.save,'figures/rlum_a.png',resolution=res
    endif
endif



;;----------------------------------------------------------------------------------------
;; NH Distribution
;;----------------------------------------------------------------------------------------
if keyword_set(nh_dist) then begin
    print, '    PREPARING NH DISTRIBUTION'
    
    wnormalize = 0
    rnormalize = 1   
    
    ;; BORUS MODEL
    ;; shaded plot boundaries
    ;; [unobscured,Compton-thick]
    model = 'BORUS'     
    nh_bound = alog10([1e21,1.5e24])
    ll_bound = rl2nh(nh_bound,model=model,/lum_out)
    
    ;; WISE AGN source detections 
    ew = {xra:[21.,25.],yra:[0.,ceil((max(yhist_wnon_borus)>max(yhist_wnon_power)>max(yhist_wdet_borus)>max(yhist_wdet_power))/10.)*10.], $
          stairstep:1,fill_background:1, $
          dimension:[1190,880], $
          font_name:'Times',font_size:14, $
          buffer:0}
    if keyword_set(hide) then ew.buffer = 1    
    ;if (wnormalize eq 1) then ew.yra = [0.,1.2]
    
    ;; WISE AGN Catalog source numbers
    ctwd_borus = commas(fix(total(yhist_wdet_borus)))
    ctwn_borus = commas(fix(total(yhist_wnon_borus)))
    
    ;; arrow x/y span
    yarr = [1.:5.]
    yarr = transpose([[yarr/10.],[yarr/10.]])
    ;if (normalize eq 0) then yarr*=ew.yra[1]
    ;xarr = transpose([[intarr(5)+alog10(1.5e24)+0.1],[intarr(5)+alog10(1.5e24)+0.3]])
    xarr = transpose([[intarr(n_elements(yarr[0,*]))-0.1],[intarr(n_elements(yarr[0,*]))+0.1]])
    print, yarr

    ;; make plot
    pw = plot(xhist_wnon_borus,yhist_wnon_borus,_extra=ew,position=[85,445,595,825],/device,/nodata)
    ;pw.xtickvalues = [ceil(e.xra[0]):floor(e.xra[1])]
    pw.axes[0].showtext=0
    ;if (normalize eq 1) then pw.ytickvalues = [0.0:1.2:0.2]
    ;; CT shading
    p = plot([alog10(1.5e24),ew.xra[1]],ew.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=ew.yra[0],fill_color='light grey',fill_transparency=0,/ov)
    ;; data
    if (wnormalize eq 1) then begin
        ywnb = nm(yhist_wnon_borus)
        ywdb = nm(yhist_wdet_borus)
    endif else begin
        ywnb = yhist_wnon_borus
        ywdb = yhist_wdet_borus
    endelse
    hwnon = plot(xhist_wnon_borus,ywnb,_extra=ew,col='orange',thick=2,linestyle='__',fill_color='orange',fill_transparency=75,/ov)
    hwdet = plot(xhist_wdet_borus,ywdb,_extra=ew,col='dodger blue',thick=2,fill_color='dodger blue',fill_transparency=75,/ov)
    re = gaussfit(xhist_wnon_borus,ywnb,a,nterms=3)
    hwlim = arrow(xarr+a[1],yarr*ew.yra[1],color='orange',/ov,/data,thick=4,fill_transparency=75,head_size=0.8,target=hwnon)
    ;; CT lines
    p = plot([1,1]*alog10(1.5e24),[ew.yra[1]*0.875,ew.yra[1]],'-',thick=2,/ov)
    p = plot([1,1]*alog10(1.5e24),[ew.yra[0],ew.yra[1]*0.125],'-',thick=2,/ov)

    ;; print sources in plot range
    t = text(0.25,0.70,'X-ray det:',col='dodger blue',target=pw,/relative,alignment=0.5,font_name='Times',font_style='Bold')
    t = text(0.25,0.65,strtrim(ctwd_borus,2)+' sources',col='dodger blue',target=pw,/relative,alignment=0.5,font_name='Times',font_style='Bold')
    t = text(0.25,0.55,'X-ray non-det:',col='orange',target=pw,/relative,alignment=0.5,font_name='Times',font_style='Bold')
    t = text(0.25,0.50,strtrim(ctwn_borus,2)+' sources',col='orange',target=pw,/relative,alignment=0.5,font_name='Times',font_style='Bold')
    ;; Catalog sources
    t = text(0.05,0.90,'!16WISE!15 AGN',col='red',target=pw,/relative,font_size=12,font_name='Times',font_style='Bold')
    t = text(0.05,0.85,'sources',col='red',target=pw,/relative,font_size=12,font_name='Times',font_style='Bold')
    xmod = text(ew.xra[0]+diff(ew.xra)/3.,ew.yra[1]*0.89,'BORUS model',target=pw,/data,font_size=14,font_name='Times',font_style='Bold')

    ;; Remaining sources
    er = {xra:[21.,25.],yra:[0.,ceil((max(yhist_rnon_borus)>max(yhist_rnon_power)>max(yhist_rdet_borus)>max(yhist_rdet_power))/200.)*200.], $
          stairstep:1,fill_background:1, $
          dimension:[1130,880], $
          font_name:'Times',font_size:14, $
          buffer:0}
    if keyword_set(hide) then er.buffer = 1    
    if (rnormalize eq 1) then er.yra = [0.,1.2]
    ;; Remaining source numbers
    ctrd_borus = commas(fix(total(yhist_rdet_borus)))
    ctrn_borus = commas(fix(total(yhist_rnon_borus)))
    
    ;; data
    if (rnormalize eq 1) then begin
        yrnb = nm(yhist_rnon_borus)
        yrdb = nm(yhist_rdet_borus)
    endif else begin
        yrnb = yhist_rnon_borus
        yrdb = yhist_rdet_borus
    endelse
    ;; make plot
    pr = plot(xhist_rnon_borus,yrnb,_extra=er,current=1,position=[595,445,1105,825],/device,/nodata)
    pr.axes[0].showtext=0
    pr.axes[1].showtext=0
    ;pr.xtickvalues = [22.:24.]
    ;; CT shading
    p = plot([alog10(1.5e24),er.xra[1]],er.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=er.yra[0],fill_color='light grey',fill_transparency=0,/ov)
    hrnon = plot(xhist_rnon_borus,yrnb,_extra=er,col='orange',thick=2,linestyle='__',fill_color='orange',fill_transparency=75,/ov)
    hrdet = plot(xhist_rdet_borus,yrdb,_extra=er,col='dodger blue',thick=2,fill_color='dodger blue',fill_transparency=75,/ov)
    a_rnon = axis('y',target=hrnon,location=[er.xra[1],0,0],textpos=1,tickdir=1,tickfont_name='Times',tickfont_size=14)    
    re = gaussfit(xhist_rnon_borus,yrnb,a,nterms=3)
    hrlim = arrow(xarr+a[1],yarr*er.yra[1],color='orange',/ov,/data,thick=4,fill_transparency=75,head_size=0.8,target=hrnon)
    ;; CT lines
    p = plot([1,1]*alog10(1.5e24),[er.yra[1]*0.875,er.yra[1]],'-',thick=2,/ov)
    p = plot([1,1]*alog10(1.5e24),[er.yra[0],er.yra[1]*0.125],'-',thick=2,/ov)
    ;ct = text(24.25,1.10,'Compton',col='black',target=pr,font_name='Times',font_style='Bold',/data)            
    ;ct = text(24.25,1.05,'thick',col='black',target=pr,font_name='Times',font_style='Bold',/data)            

    ;; sources in plot range
    t = text(0.25,0.70,'X-ray det:',col='dodger blue',target=pr,/relative,alignment=0.5,font_name='Times',font_style='Bold')
    t = text(0.25,0.65,strtrim(ctrd_borus,2)+' sources',col='dodger blue',target=pr,/relative,alignment=0.5,font_name='Times',font_style='Bold')
    t = text(0.25,0.55,'X-ray non-det:',col='orange',target=pr,/relative,alignment=0.5,font_name='Times',font_style='Bold')
    t = text(0.25,0.50,strtrim(ctrn_borus,2)+' sources',col='orange',target=pr,/relative,alignment=0.5,font_name='Times',font_style='Bold')
    ;; Catalog sources
    t = text(0.05,0.90,'Remaining',col='red',target=pr,/relative,font_size=12,font_name='Times',font_style='Bold')
    t = text(0.05,0.85,'sources',col='red',target=pr,/relative,font_size=12,font_name='Times',font_style='Bold')
    
    ;; POWER-LAW MODEL
    ;; shaded plot boundaries
    ;; [unobscured,Compton-thick]
    model = 'POWER'     
    nh_bound = alog10([1e21,1.5e24])
    ll_bound = rl2nh(nh_bound,model=model,/lum_out)

    ;; WISE AGN Catalog sources
    pw = plot(xhist_wnon_power,yhist_wnon_power,_extra=ew,current=1,position=[85,65,595,445],/device,/nodata)
    pw.xtickvalues = [21.:24.:1.]
    pw.ytickvalues = pw.ytickvalues[0:-2]
    ;if (normalize eq 1) then pw.ytickvalues = [0.0:1.:0.2]
    
    ;; CT shading
    p = plot([alog10(1.5e24),ew.xra[1]],ew.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=ew.yra[0],fill_color='light grey',fill_transparency=0,/ov)
    ;; data
    if (wnormalize eq 1) then begin
        ywnp = nm(yhist_wnon_power)
        ywdp = nm(yhist_wdet_power)
    endif else begin
        ywnp = yhist_wnon_power
        ywdp = yhist_wdet_power
    endelse
    hwnon = plot(xhist_wnon_power,ywnp,_extra=ew,col='orange',thick=2,linestyle='__',fill_color='orange',fill_transparency=75,/ov)
    hwdet = plot(xhist_wdet_power,ywdp,_extra=ew,col='dodger blue',thick=2,fill_color='dodger blue',fill_transparency=75,/ov)
    re = gaussfit(xhist_wnon_power,yhist_wnon_power,a,nterms=3)
    hwlim = arrow(xarr+a[1],yarr*ew.yra[1],color='orange',/ov,/data,thick=4,fill_transparency=75,head_size=0.8,target=hwnon)
    ;; CT lines
    p = plot([1,1]*alog10(1.5e24),[ew.yra[1]*0.875,ew.yra[1]],'-',thick=2,/ov)
    p = plot([1,1]*alog10(1.5e24),[ew.yra[0],ew.yra[1]*0.125],'-',thick=2,/ov)
    xmod = text(ew.xra[0]+diff(ew.xra)/3.,ew.yra[1]*0.89,'Power law model',target=pw,/data,font_size=14,font_name='Times',font_style='Bold')
    
    ;; Remaining sources
    ctrd_power = commas(fix(total(yhist_rdet_power)))
    ctrn_power = commas(fix(total(yhist_rnon_power)))
    
    ;; make plot
    pr = plot(xhist_rnon_power,yhist_rnon_power,_extra=er,current=1,position=[595,65,1105,445],/device,/nodata)
    pr.xtickvalues = [22.:25.:1.]
    pr.axes[1].showtext=0
    ;; CT shading
    p = plot([alog10(1.5e24),er.xra[1]],er.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=er.yra[0],fill_color='light grey',fill_transparency=0,/ov)
    if (rnormalize eq 1) then begin
        yrnp = nm(yhist_rnon_power)
        yrdp = nm(yhist_rdet_power)
    endif else begin
        yrnp = yhist_rnon_power
        yrdp = yhist_rdet_power
    endelse
    hrnon = plot(xhist_rnon_power,yrnp,_extra=er,col='orange',thick=2,linestyle='__',fill_color='orange',fill_transparency=75,/ov,name=' X-ray non-det.')
    hrdet = plot(xhist_rdet_power,yrdp,_extra=er,col='dodger blue',thick=2,fill_color='dodger blue',fill_transparency=75,/ov,name=' X-ray detected')
    a_rnon = axis('y',target=hrnon,location=[er.xra[1],0,0],textpos=1,tickdir=1,tickfont_name='Times',tickfont_size=14)    
    a_rnon.tickvalues = a_rnon.tickvalues[0:-2]
    re = gaussfit(xhist_rnon_power,yrnp,a,nterms=3)
    hrlim = arrow(xarr+a[1],yarr*er.yra[1],color='orange',/ov,/data,thick=4,fill_transparency=75,head_size=0.8,target=hrnon)
    ;; CT lines
    p = plot([1,1]*alog10(1.5e24),[er.yra[1]*0.875,er.yra[1]],'-',thick=2,/ov)
    p = plot([1,1]*alog10(1.5e24),[er.yra[0],er.yra[1]*0.125],'-',thick=2,/ov)
    ;; CT text
    ct = text(0.805,0.09,'Compton',col='black',target=pr,font_size=12,font_name='Times',font_style='Bold',/relative)          
    ct = text(0.805,0.03,'thick',col='black',target=pr,font_size=12,font_name='Times',font_style='Bold',/relative)            

    xt = text(0.5,0.032,'$log  !8N!7_{H} [cm^{-2}]$',alignment=0.5,vertical_alignment=0.5,font_size=14, font_name='Times')            
    if (wnormalize eq 1) then yt = text(0.02,0.5,'Frequency (!8WISE!7 AGN sources) [normalized]',orientation=90.,alignment=0.5,vertical_alignment=0.5,font_size=14,font_name='Times') else $
                              yt = text(0.02,0.5,'Frequency (!8WISE!7 AGN sources)',orientation=90.,alignment=0.5,vertical_alignment=0.5,font_size=14,font_name='Times')
    if (rnormalize eq 1) then yt2 = text(0.975,0.5,'Frequency (remaining sources) [normalized]',orientation=90.,alignment=0.5,vertical_alignment=0.5,font_size=14,font_name='Times') else $
                              yt2 = text(0.975,0.5,'Frequency (remaining sources)',orientation=90.,alignment=0.5,vertical_alignment=0.5,font_size=14,font_name='Times')

    l = legend(target=[hrdet,hrnon],position=[0.70,0.49],/normal,/auto_text_color,sample_width=0.1,horizontal_spacing=0.06,font_name='Times')

    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/nh_dist.eps',/BITMAP else $
                                                     p.save,'figures/nh_dist.png',resolution=res
                             
    endif
endif


END




