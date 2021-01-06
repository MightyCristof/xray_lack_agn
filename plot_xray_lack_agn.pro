PRO plot_xray_lack_agn, PHOT_SPEC = phot_spec, $
                        Z_DIST = z_dist, $
                        FX_LIMIT = fx_limit, $
                        SEDS = seds, $
                        LX_LIR = lx_lir, $
                        LUM_RATIO = lum_ratio, $
                        NH_DIST = nh_dist, $
                        SURV_ANAL = surv_anal, $
                        XSTACK = xstack, $
                        GSTACK = gstack, $
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
common _wac 
common _xconv   
common _fxlim 
common _agnlum 
common _quality  
common _combined
common _nhdist
common _surv
common _xstack


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
         font_name:'Times',font_size:16, $
         dimension:[700,800],buffer:0}
    if keyword_set(hide) then e.buffer = 1
    
    col = [[204,121,167],[213,94,0],[0,158,115],[0,114,178],[240,228,66]]
    col_ex = col[*,3]
    col_in = col[*,4]

    p = plot(zs_ex,zp_ex,'o',color=col_ex,_extra=e)
    p.position = [0.14735492,0.33393553,0.789520,0.83583007]
    p = plot([0.,1.],[0.,1.],'--',thick=2,/ov)
    p = plot(zs_in,zp_in,'o',color=col_in,sym_size=0.5,sym_filled=1,/ov)
    p = plot(zs_in,zp_in,'o',color=black,sym_size=0.5,sym_thick=1.,sym_filled=0,transparency=90,/ov)
    iexoff = where(zs_ex gt e.xra[1],nexoff)
    iinoff = where(zs_in gt e.xra[1],ninoff)
    p.axes[0].showtext = 0
    num_str = '$Num. matches = '+commas(strtrim(nex,2))+' ('+commas(strtrim(nin,2))+')'+'$'
    t = text(target=p,e.xra[1]*0.95,e.yra[1]*0.05,num_str,/data,alignment=1,font_size=14,font_style='Bold',font_name='Times',fill_background=1,fill_color='white')
    ;; PLOT SecondaryS
    er = {xra:e.xra,yra:[-0.8,0.8],aspect_ratio:0, $
          sym_size:0.5,sym_filled:1,transparency:75, $
          xtitle:'$!8z!7_{spec} (MMT)$',ytitle:'$\Delta!8z!7 / (1+!8z!7_{spec})$', $
          font_name:'Times',font_size:16, $
          dimension:[700,800],buffer:0}
    pr = plot(zs_ex,delz_ex,'o',color=col_ex,_extra=er,/current)
    pr.position = [p.position[0],0.13195312,p.position[2],p.position[1]]
    pr = plot(p.xra,[0,0],'--',thick=2,/ov)
    pr = plot(zs_in,delz_in,'o',color=col_in,sym_size=0.5,sym_filled=1,/ov)
    pr = plot(zs_in,delz_in,'o',color=black,sym_size=0.5,sym_thick=1.,sym_filled=0,transparency=90,/ov)
    pr.ytickvalues = [-0.4,0.,0.4]
    pr.yminor = 3.
    zsig_str = '$\sigma_{\Delta z} = '+string(rnd(zexsig,5),format='(d5.3)')+' ('+string(rnd(zinsig,5),format='(d5.3)')+')'+'$'
    t = text(target=pr,er.xra[1]*0.95,er.yra[1]*0.62,zsig_str,/data,alignment=1,font_size=14,font_style='Bold',font_name='Times',fill_background=1,fill_color='white')
    ;; PLOT SecondaryS HISTOGRAM
    yh_ex = histogram(delz_ex,bin=scott(delz_ex),locations=xh_ex)
    yh_in = histogram(delz_in,bin=scott(delz_in),locations=xh_in)
    eh = {xra:[0.,ceil(max(yh_ex)/100.+1)*100.],yra:[-0.3,0.3], $
          stairstep:1,fill_background:1,fill_transparency:30, $
          xtitle:'Frequency', $
          font_name:'Times',font_size:16}
    ph = plot(yh_ex,xh_ex,/current,_extra=eh,fill_color=col_ex)
    ph = plot(yh_in,xh_in,/ov,_extra=eh,fill_color=col_in)
    ph = plot(ph.xra,[0.,0.],'--',thick=2,/ov)
    ph.position = [p.position[2],pr.position[1],0.93303385,pr.position[3]]
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
            stairstep:1,fill_background:1,fill_transparency:30, $
            ytitle:'Frequency', $
            font_name:'Times',font_size:16}
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
            stairstep:1,fill_background:1,fill_transparency:30, $
            font_name:'Times',font_size:16}
    ph_zp = plot(yhzp_ex,xhzp_ex,/current,_extra=eh_zp,fill_color=col_ex)
    ph_zp = plot(yhzp_in,xhzp_in,/ov,_extra=eh_zp,fill_color=col_in)
    ph_zp.position = [p.position[2],p.position[1],ph.position[2],p.position[3]]
    ph_zp.axes[1].showtext = 0
    ph_zp.axes[0].showtext = 0
    ph_zp.xtickvalues = ph.xtickvalues
    ph_zp.xminor = 3.
        
    ;print, 'NUM. SDSS:        ' + strtrim(nex,2)
    ;print, 'NUM. SAMP:        ' + strtrim(nin,2)
    ;print, 'NUM. SDSS MISSES: ' + strtrim(nexoff,2)
    ;print, 'NUM. SAMP MISSES: ' + strtrim(ninoff,2)
    ;print, 'SIGMA ZSDSS:      ' + strtrim(zexsig,2)
    ;print, 'SIGMA ZSAMP:      ' + strtrim(zinsig,2)
        
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
if keyword_set(z_dist) then begin

    izs = where(strmatch(ztype,'ZS*'))
    izp = where(strmatch(ztype,'ZP*'))
    ;izx = where(strmatch(ztype,'ZPXD'))
    
    zbin = 0.02d
    yz = [0.,histogram(z,bin=zbin,location=xz,min=0.,max=max(z)),0.]
    yzs = [0.,histogram(z[izs],bin=zbin,location=xzs,min=0.,max=max(z)),0.]
    yzp = [0.,histogram(z[izp],bin=zbin,location=xzp,min=0.,max=max(z)),0.]
    xz = [-0.02d,xz,max(z)+zbin]+zbin/2.
    xzs = [-0.02d,xzs,max(z)+zbin]+zbin/2.
    xzp = [-0.02d,xzp,max(z)+zbin]+zbin/2.

    col = [[204,121,167],[213,94,0],[0,158,115],[0,114,178],[240,228,66]]
    col_zp = col[*,3]
    col_zs = col[*,4]

    ;; percent of redshift
    pzp = string(rnd(total(strmatch(ztype,'ZP*'))/nsrc*100.,0),format='(i2)')+'%'
    pzs = string(rnd(total(strmatch(ztype,'ZS*'))/nsrc*100.,0),format='(i2)')+'%'
    nzp = commas(string(total(strmatch(ztype,'ZP*')),format='(i0)'))
    nzs = commas(string(total(strmatch(ztype,'ZS*')),format='(i0)'))
    
    e = {xtitle:'$!8z!7$',ytitle:'Frequency', $
         font_size:16,font_name:'Times', $
         xra:[-zbin,max(z)+zbin], $
         yra:[0.,ceil((max(yz))/100.)*100.], $
         stairstep:1,thick:2, $
         fill_background:1,fill_transparency:30}
    dim = [800,700]
    mrg = [100,65]
    p = plot(xz,yz,_extra=e,linestyle='-',dim=dim,/device,position=[mrg[0],mrg[1],(dim-mrg)[0]+30,(dim-mrg)[1]],/nodata)
    pzp = plot(xzp,yzp,_extra=e,fill_color=col_zp,linestyle='-',/ov,name='$!8z!7_{phot}'+' ('+nzp+'-'+pzp+'$)');'!8z!7$_{phot}$ (SDSS+!8XDQSOz!7)')    
    pzs = plot(xzs,yzs,_extra=e,fill_color=col_zs,linestyle='__',/ov,name='$!8z!7_{spec}'+' ('+nzs+'-'+pzs+'$)');'!8z!7$_{spec}$ (SDSS+AGES+supp.)')
    l = legend(target=[pzp,pzs],position=[0.84,0.94],/relative,sample_width=0.16,horizontal_spacing=0.06,horizontal_alignment=1.,vertical_alignment=1.,font_name='Times',font_size=14)
    
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/z_dist.eps',/BITMAP else $
                                                     p.save,'figures/z_dist.png',resolution=res
    endif
endif



;;----------------------------------------------------------------------------------------
;; Flux Limit
;;----------------------------------------------------------------------------------------
if keyword_set(fx_limit) then begin
    multi_sn = 1
    print, '    PREPARING PAPER FLUX LIMIT'
    dim = [1280,510]
    mrg = [90,65]
    src_pos = [[mrg[0],mrg[1],mrg[0]+380,(dim-mrg)[1]],[mrg[0]+380,mrg[1],mrg[0]+380*2.,(dim-mrg)[1]],[850,mrg[1],mrg[0]+380*3.,(dim-mrg)[1]]]
    
    e = {yra:[-15.5,-10.], $
         xlog:1,ylog:0, $
         xtitle:'$exposure time [s]$',ytitle:'$log !8F!7_{2-10 keV} [erg s^{-1}cm^{-2}]$', $
         font_name:'Times',font_size:16, $ $
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
        perr = plot([xpt[i]],ypt,'S',col=[12,123,220],sym_size=1.,sym_filled=0,/ov)
        ;; text
        t = text(0.06,0.88,cat[i],target=p,/RELATIVE,font_size=14,font_style='Bold italic',alignment=0.,font_name='Times')

        re = execute('ncat = n_elements(CAT_EXP'+xfield[i]+')')
        re = execute('irand = lindgen(n_elements(CAT_EXP'+xfield[i]+'))')
        re = execute('pcat_lite = plot(CAT_EXP'+xfield[i]+'[irand],alog10(CAT_FLX'+xfield[i]+'[irand]),"s",col="light grey",sym_filled=1,sym_size=0.5,transparency=85,/ov)')
        ; plot open symbols first (makes for a cleaner plot)
        re = execute('pnon = plot(TEXP'+xfield[i]+'[where(iiqual_non'+xfield[i]+')],logfxir[where(iiqual_non'+xfield[i]+')],"o",sym_filled=1,sym_size=0.5,color="orange",transparency=50,/ov,name="X-ray non-det. (expected)")')
        re = execute('pdet = plot(EXP'+xfield[i]+'[where(iiqual_det'+xfield[i]+')],alog10(FX'+xfield[i]+'[where(iiqual_det'+xfield[i]+')]),"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov,name="X-ray detected (observed)")')
        re = execute('pdo = plot(EXP'+xfield[i]+'[where(iiqual_det'+xfield[i]+')],alog10(FX'+xfield[i]+'[where(iiqual_det'+xfield[i]+')]),"S",sym_filled=0,sym_size=1.,transparency=90,/ov)')
        ;; stretch X-ray flux limit for plot
        re = execute('lxr = alog10([min(EXP'+xfield[i]+'[where(EXP'+xfield[i]+' gt 0.)])<min(TEXP'+xfield[i]+'[where(TEXP'+xfield[i]+' gt 0.)]),max(EXP'+xfield[i]+'[where(EXP'+xfield[i]+' gt 0.)])>max(TEXP'+xfield[i]+'[where(TEXP'+xfield[i]+' gt 0.)])])')
        xx = [lxr[0]:lxr[1]:diff(minmax(lxr))/1000.]
        xx_str = '['+strtrim(indgen(degr[i]+1),2)+']*xx^'+strtrim(indgen(degr[i]+1),2)+'.'
        re = execute('fxlim = FXLIM_CS'+xfield[i])
        re = execute('yy = '+strjoin(+'fxlim'+xx_str,' + '))
        pxlim = plot(10.^xx,yy,'--r',thick=2,/ov,name="X-ray flux limit")
        re = execute('ilim = where(CAT_EXP5'+xfield[i]+' ge min(10.^xx) and CAT_EXP5'+xfield[i]+' le max(10.^xx))')
        re = execute('p = plot(CAT_EXP5'+xfield[i]+'[ilim],alog10(CAT_LIM5'+xfield[i]+'[ilim]),"--",thick=2,/ov)')
        if (i eq nfield-1) then begin
            snt = text(2.5e3,-12.3,'$SNR \geq  5.0$',target=p,/DATA,font_size=14,font_name='Times',font_style='Bold')
            snt = text(2.5e3,-13.3,'$SNR \geq  3.0$',target=p,/DATA,font_size=14,font_name='Times',font_style='Bold',col='red')
        endif
        
        ;; ensure axes range
        p.xra = xra[*,i]
        
        ;if (i eq nfield-1) then begin
        ;    re = execute('pright = plot(TEXP'+xfield[i]+'[where(iiqual_non'+xfield[i]+')],logfxir[where(iiqual_non'+xfield[i]+')],/nodata,/current,position=src_pos[*,i],axis_style=0)')
        ;    aright = axis('y',target=pnon,location=[pnon.xra[1],0,0],textpos=1,tickdir=1,title='$log !8F!7_{2-10 keV} (expected) [erg s^{-1}cm^{-2}]$',tickfont_name='Times',tickfont_size=16)
        ;endif
    endfor
        
    l = legend(target=[pdet,pnon,pcat_dark],position=[0.98,0.85],/relative,sample_width=0,horizontal_spacing=0.06,vertical_alignment=0.5,horizontal_alignment=1.0)

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
    print, 'WISE AGN:', iiwac[inds]
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
    ;; strings for plot
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
    e = {xr:[0.05,30.],yra:[-15.5,-8.5], $
         xlog:1,thick:2,font_name:'Times',font_size:16, $
         buffer:0}
    dim = [880,968]
    mrg = [80,65]
    psp = (dim[1]-2.*mrg[1])/3.
    pos = [[mrg[0]+20,mrg[1]+psp*2,dim[0]/2+20,(dim-mrg)[1]], [dim[0]/2+20,mrg[1]+psp*2,(dim-mrg)[0]+20,(dim-mrg)[1]], $
           [mrg[0]+20,mrg[1]+psp,dim[0]/2+20,mrg[1]+psp*2],   [dim[0]/2+20,mrg[1]+psp,(dim-mrg)[0]+20,mrg[1]+psp*2], $
           [mrg[0]+20,mrg[1],dim[0]/2+20,mrg[1]+psp],         [dim[0]/2+20,mrg[1],(dim-mrg)[0]+20,mrg[1]+psp]]       
    label = transpose([['RA: '+plot_ra],['Dec: '+plot_dec],['$!8z!7: $'+plot_z+'$\pm$'+plot_zerr],['$!8E!7(!8B-V!7)$: '+plot_ebv+'$\pm$'+plot_ebvsig],['$\chi^2 / DoF$: '+plot_chi]])

    if keyword_set(hide) then e.buffer = 1           
    for i = 0,nplt-1 do begin
        if (i eq 0) then current = 0 else current = 1
        ig = where(plot_bin[*,i],/null)
        p = plot(objwav[ig,i],plot_flux[ig,i],_extra=e,dim=dim,position=pos[*,i],/DEVICE,current=current,/nodata)
        if (i lt 4) then p.axes[0].showtext = 0
        if (i mod 2) then p.axes[1].showtext = 0
        p_agn = plot(tempwav[*,i],agn[*,i],col=col[*,0],_extra=e,linestyle=linestyle[0],/ov,name=' AGN')
        p_ell = plot(tempwav[*,i],ell[*,i],col=col[*,1],_extra=e,linestyle=linestyle[1],/ov,name=' ELL')
        p_sfg = plot(tempwav[*,i],sfg[*,i],col=col[*,2],_extra=e,linestyle=linestyle[2],/ov,name=' SFG')
        p_irr = plot(tempwav[*,i],irr[*,i],col=col[*,3],_extra=e,linestyle=linestyle[3],/ov,name=' IRR')
        p = plot(tempwav[*,i],model[*,i],col='dark slate grey',thick=2,/ov)                                                           ;; plot coadded models
        p = errorplot(objwav[ig,i],plot_flux[ig,i],plot_err[ig,i],'o',SYM_FILLED=1,LINESTYLE='',sym_size=1.5,errorbar_thick=2,/OV)          ;; plot photometry
        for t = 0,n_elements(label[*,0])-1 do begin
            lab = text(0.08,e.yra[1]-0.2-(1.5+t)/2.4,label[t,i],font_size=14,/DATA,target=p,font_name='Times')
        endfor
        if (i eq 1) then l = legend(target=[p_agn,p_ell,p_sfg,p_irr],position=[0.98,0.94],/relative,sample_width=0.09,horizontal_spacing=0.02,font_name='Times',font_size=14)
        
        if (i eq 2) then p.ytitle = '$log \nu!8F!7_\nu  [erg s^{-1}cm^{-2}]$'
        ;if (i eq 5) then p.xtitle = '$\lambda (observed) [ \mum ]$'
    endfor
    xt = text(0.5,0.016,'$\lambda (observed) [ \mum ]$',/relative,alignment=0.5,font_size=16,font_name='Times')
    
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
    ;; Stern+15 LX-LIR
    xrel_stern = [40.:50.:0.01]
    yrel_stern = dblarr(n_elements(xrel_stern))
    yrel_stern[*] = 40.981 + 1.024*(xrel_stern-41.) - 0.047*(xrel_stern-41.)^2.
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

    ;; quality sources
    ;iqso1 = where(iiqual_det and ebv le 0.15)
    ;iqso2 = where(iiqual_det and ebv gt 0.15)
    ;; final sample
    iqso1 = where(iilir and iidet and ebv le 0.15)
    iqso2 = where(iilir and iidet and ebv gt 0.15)
    
    ;; test (1+z) correction
    ;loglirz = alog10((10.^loglir)/(1+z))
    ;loglirz[where(loglir eq -9999)] = -9999
    ;loglir = loglirz

    lxbinsz = 0.15
    hist1 = hist_2d(loglir[iqso1],loglx[iqso1],bin1=lxbinsz,min1=41.,max1=47.5,bin2=lxbinsz,min2=41.,max2=46.5)
    hist2 = hist_2d(loglir[iqso2],loglx[iqso2],bin1=lxbinsz,min1=41.,max1=47.5,bin2=lxbinsz,min2=41.,max2=46.5)
    ind1 = array_indices(hist1,where(hist1 eq max(hist1)))
    ind2 = array_indices(hist2,where(hist2 eq max(hist2)))
    xlum = [41.:47.5:lxbinsz]
    ylum = [41.:46.5:lxbinsz]
    
    e = {xra:[41.,47.5],yra:[41.,46.5], $
         font_name:'Times',font_size:16, $
         xtitle:'$log !8L!7_{6 \mu m} [erg s^{-1}]$', $
         buffer:0}
    if keyword_set(hide) then e.buffer = 1
    rel_lin = ['-.','-','__']
    ;rel_col = [[105,105,105],[0,0,0]]
    ;rel_col = [[0,158,115],[204,121,167],[213,94,0]]
    rel_col = [[30,39,44],[30,39,44],[30,39,44]]
    ;; [0,158,115]  [204,121,167]
    im1 = image('chen17_bw1.png',transparency=50,dimension=[640,877],position=[80,440,587,810],buffer=buff,/device)
    pnodata = plot(xrel_chen,yrel_chen,/nodata,_extra=e,/current,position=im1.position)
    cqso1 = contour(hist1,xlum,ylum,c_thick=4,rgb_table=colortable(49,/reverse),/fill,/ov,c_label_show=0,transparency=0,name='$!8E(B-V)!7 \leq  0.15$')
    prel_chen = plot(xrel_chen,yrel_chen,linestyle=rel_lin[0],col=rel_col[*,0],thick=4,_extra=e,/ov,position=im1.position)
    prel_stern = plot(xrel_stern,yrel_stern,linestyle=rel_lin[1],col=rel_col[*,1],thick=4,_extra=e,/ov)
    prel_fiore = plot(xrel_fiore,yrel_fiore,linestyle=rel_lin[2],col=rel_col[*,2],thick=4,/ov)
    ;chent = text(41.15,42.3,/data,target=prel_chen,'Chen+17',col=rel_col[*,0],font_size=14,font_style='Bold',font_name='Times',fill_background=1,fill_color='white')
    ;sternt = text(41.15,42.3,/data,target=prel_chen,'Stern+15',col=rel_col[*,1],font_size=14,font_style='Bold',font_name='Times',fill_background=1,fill_color='white')
    ;fioret = text(41.95,41.3,/data,target=prel_fiore,'Fiore+09',col=rel_col[*,2],font_size=14,font_style='Bold',font_name='Times',fill_background=1,fill_color='white')
    prel_chen.axes[0].showtext = 0
    ;p = plot([median(xlum[ind1[0,*]])],[median(ylum[ind1[1,*]])],'o',col='black',sym_size=1.5,/ov)
    ;p = plot([median(xlum[ind1[0,*]])],[median(ylum[ind1[1,*]])],'S',col='black',sym_size=1.0,sym_filled=1,/ov)
    ;p = plot(xlum[ind1[0,*]],ylum[ind1[1,*]],'X',sym_size=1.5,sym_thick=1,/sym_filled,col=[0,158,115],/ov)
    im2 = image('chen17_bw2.png',transparency=50,/current,position=[80,75,587,445],/device)
    pnodata = plot(xrel_chen,yrel_chen,/nodata,_extra=e,/current,position=im2.position)    
    cqso2 = contour(hist2,xlum,ylum,c_thick=4,rgb_table=colortable(62,/reverse),/fill,/ov,c_label_show=0,transparency=0,name='$!8E(B-V)!7 > 0.15$')
    prel_chen = plot(xrel_chen,yrel_chen,linestyle=rel_lin[0],col=rel_col[*,0],thick=4,_extra=e,/ov,position=im2.position,name='Chen+2017')    
    prel_stern = plot(xrel_stern,yrel_stern,linestyle=rel_lin[1],col=rel_col[*,1],thick=4,_extra=e,/ov,name="Stern 2015")
    prel_fiore = plot(xrel_fiore,yrel_fiore,linestyle=rel_lin[2],col=rel_col[*,2],thick=4,/ov,name='Fiore+2009')
    ;p = plot([median(xlum[ind2[0,*]])],[median(ylum[ind2[1,*]])],'o',col='black',sym_size=1.5,/ov)
    ;p = plot([median(xlum[ind2[0,*]]),0,0,0],[median(ylum[ind2[1,*]]),0,0,0],'S',col='black',sym_size=1.0,sym_filled=1,_extra=e,/ov)
    ;p = plot(xlum[ind2[0,*]],ylum[ind2[1,*]],'X',sym_size=1.5,sym_thick=3,/sym_filled,col=[0,158,115],/ov)
    leg_qso1 = legend(target=[cqso1],position=[0.85,0.53],/normal,vertical_alignment=0.,horizontal_alignment=1.,font_size=14,font_name='Times')   
    leg_qso2 = legend(target=[cqso2],position=[0.85,0.115],/normal,vertical_alignment=0.,horizontal_alignment=1.,font_size=14,font_name='Times')   
    leg_rel = legend(target=[prel_chen,prel_stern,prel_fiore],position=[0.16,0.39],/normal,sample_width=0.12,horizontal_spacing=0.06,vertical_alignment=0.,horizontal_alignment=0.,font_size=14,font_name='Times')   
    ;leg.position = [0.14,0.39]

    ;xt = text(0.52,0.03,'$log !8L!7_{6 \mu m} [erg s^{-1}cm^{-2}]$',alignment=0.5,font_size=14,font_name='Times')
    yt = text(0.040,0.51,'$log !8L!7_{2-10 keV} [erg s^{-1}]$',orientation=90,alignment=0.5,font_size=16,font_name='Times')
        
    ;peak_circ = text(0.72,0.13,'$\U(25EF)$',col='dark slate grey',font_size=14,font_name='Times')
    ;peak_star = text(0.721,0.128,'$\U(2605)$',col='dark slate grey',font_size=14,font_name='Times')
    ;peak_text = text(0.76,0.128,'Peak Freq.',col='dark slate grey',font_size=14,font_name='Times')
    
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

    nh = [21.:25.:0.25]
    nh_values = alog10([1e23,5e23,1e24,1.5e24,1e25])
    nhtext = ['23.00','23.70','24.00','24.18','25.00']
    
    models = ['POWER','BORUS']
    modlab = ['Power Law Model','BORUS Model']
    savfile = ['rlum_a','rlum_b']
    print, n_elements(models)
    for i = 0,n_elements(models)-1 do begin
        model = models[i]
        ll = rl2nh(nh,model=model,/rl_out)
        re = execute('nhdet = NHDET_'+model)
        re = execute('e_nhdet = E_NHDET_'+model)
        re = execute('nhnon = NHNON_'+model)

        ;; RL
        ell = {xra:[-3.3,2.2],yra:[-2.9,1.], $
               xtickformat:'(i0)',ytickformat:'(i0)', $
               font_name:'Times',font_size:16, $
               buffer:0}
        dim = [1000,800]
        mrg = [80,65]
        if keyword_set(hide) then ell.buffer = 1

        nh_lines = [ell.xra[0],ell.xra[0]+0.7]
        ll_lines = rebin(interpol(ll,nh,nh_values),n_elements(nh_values),2)

        ;; WISE AGN
        pw = plot(logebv[where(iiqual_non)],llnon[where(iiqual_non)],_extra=ell,dim=dim,/device,position=[mrg[0],dim[1]/2,dim[0]/2,(dim-mrg)[1]],/nodata,ytitle='$log !8L!7_X / !8L!7_X(!8L!7_{MIR})  (2$-$10 keV)$')
        ;; shading
        nh_bound = alog10([1e21,1.5e24])
        ll_bound = rl2nh(nh_bound,model=model,/rl_out)
        pw = plot(ell.xra,[1,1]*ell.yra[1],linestyle='',fill_background=1,fill_level=0.,fill_color='light grey',fill_transparency=0,/ov)
        pw = plot(ell.xra,ll_bound[1]*[1,1],linestyle='',fill_background=1,fill_level=ell.yra[0],fill_color='light grey',fill_transparency=0,/ov)
        ;; plot data
        pwllnon = arrow(transpose([[logebv[where(iiqual_non and iiwac)]],[logebv[where(iiqual_non and iiwac)]]]),transpose([[llnon[where(iiqual_non and iiwac)]],[llnon[where(iiqual_non and iiwac)]-0.1]]),color="orange",fill_transparency=75,target=[pw],/data,head_size=0.3,clip=1)
        pwlldet = plot(logebv[where(iiqual_det and iiwac)],lldet[where(iiqual_det and iiwac)],"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov,name=" X-ray detected")
        p = plot(logebv[where(iiqual_det and iiwac)],lldet[where(iiqual_det and iiwac)],"S",sym_filled=0,sym_size=1.5,transparency=85,/ov)
        ;; median error bars
        perr = plot([-2.2],[ell.yra[0]-(ell.yra[0]-ll_bound[1])/2.],'s',sym_size=4,sym_thick=2,/ov)
        perr = errorplot([-2.2],[ell.yra[0]-(ell.yra[0]-ll_bound[1])/2.],[median(e_lldet[where(iiqual_det and iiwac)])],errorbar_capsize=0.1,linestyle="",/ov)
        perr = plot([-2.2],[ell.yra[0]-(ell.yra[0]-ll_bound[1])/2.],"S",sym_filled=1,sym_size=1.,col='white',/ov)
        perr = plot([-2.2],[ell.yra[0]-(ell.yra[0]-ll_bound[1])/2.] ,"S",sym_filled=0,sym_size=1.,col=[12,123,220],/ov)
        ;; text
        t = text(ell.xra[1]-diff(ell.xra)/2.,ell.yra[0]+0.5,modlab[i],target=pwlldet,alignment=0.5,vertical_alignment=0.5,font_size=14,font_name='Times',font_style='Bold',/data)
        t = text(0.05,0.88,'!16WISE!15 AGN',col='black',/relative,target=[pwllnon],font_size=14,font_name='Times',font_style='Bold')
        t = text(ell.xra[1]-0.2,ll_bound[0]+0.12,'unobscured',col='black',target=pwllnon,alignment=1.,font_size=14,font_name='Times',font_style='Bold',/data)
        t = text(ell.xra[1]-0.2,ll_bound[1]-0.12,'CT',col='black',target=pwllnon,alignment=1.,vertical_alignment=1.,font_size=14,font_name='Times',font_style='Bold',/data)
        t = text(nh_lines[1],ll_bound[0]+0.1,'$log !8N!7_H$',col='black',/data,alignment=1.,target=[pwllnon],font_size=11,font_name='Times')
        for n = 0,n_elements(nh_values)-1 do begin
            t = text(nh_lines[1],ll_lines[n,*]+0.05,nhtext[n],col='black',/data,alignment=1.,target=[pwllnon],font_size=11,font_name='Times')
            p = plot(nh_lines,ll_lines[n,*],linestyle=':',col='black',thick=2,/ov)
        endfor
        ;; remove axes
        pwlldet.axes[0].showtext=0
        ;; model choice

        ;; SECONDARY SOURCES
        pr = plot(logebv[where(iiqual_non)],llnon[where(iiqual_non)],_extra=ell,position=[dim[0]/2,dim[1]/2,(dim-mrg)[0],(dim-mrg)[1]],current=1,/device,/nodata)
        pr.axes[1].showtext = 0 & pr.axes[0].showtext = 0
        ;; shading
        pr = plot(ell.xra,[1,1]*ell.yra[1],linestyle='',fill_background=1,fill_level=0.,fill_color='light grey',fill_transparency=0,/ov)
        pr = plot(ell.xra,ll_bound[1]*[1,1],linestyle='',fill_background=1,fill_level=ell.yra[0],fill_color='light grey',fill_transparency=0,/ov)
        ;; plot data
        prllnon = arrow(transpose([[logebv[where(iiqual_non and ~iiwac)]],[logebv[where(iiqual_non and ~iiwac)]]]),transpose([[llnon[where(iiqual_non and ~iiwac)]],[llnon[where(iiqual_non and ~iiwac)]-0.1]]),color="orange",fill_transparency=75,target=[pr],/data,head_size=0.3,clip=1)
        prlldet = plot(logebv[where(iiqual_det and ~iiwac)],lldet[where(iiqual_det and ~iiwac)],"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov,name=" X-ray detected")
        p = plot(logebv[where(iiqual_det and ~iiwac)],lldet[where(iiqual_det and ~iiwac)],"S",sym_filled=0,sym_size=1.5,transparency=85,/ov)
        prllnon = plot(logebv[where(iiqual_non and ~iiwac)],llnon[where(iiqual_non and ~iiwac)],"td",sym_filled=1,sym_size=1.,color="orange",/current,/nodata,name=" X-ray non-det.",position=[dim[0]/2,dim[1]/2,(dim-mrg)[0],(dim-mrg)[1]],_extra=ell,axis_style=0)
        a_rllnon = axis('y',target=pr,location=[ell.xra[1],0,0],textpos=1,tickdir=1,title='$log !8L!7_{X,lim} / !8L!7_X(!8L!7_{MIR})  (2$-$10 keV)$',tickformat='(i0)',tickfont_name='Times',tickfont_size=16)
        ;; text
        t = text(0.05,0.88,0.50,'Secondary',col='black',/relative,target=[prlldet],font_size=14,font_name='Times',font_style='Bold')

        ;; NH
        enh = {xra:ell.xra,yra:[20.8,25.2], $
               xtickformat:'(i0)',ytickformat:'(i0)', $
               font_name:'Times',font_size:16, $
               buffer:0}
        
        ;; WISE AGN Catalog sources
        pw = plot(logebv[where(iiqual_non and iiwac)],llnon[where(iiqual_non and iiwac)],_extra=enh,current=1,/device,position=[mrg[0],mrg[1],dim[0]/2,dim[1]/2],/nodata,ytitle='$log !8N!7_H [cm^{-2}]$')
        ;; shading
        pw = plot(enh.xra,nh_bound[1]*[1,1],linestyle='',fill_background=1,fill_level=enh.yra[1],fill_color='light grey',fill_transparency=0,/ov)
        ;; plot data
        pwnhnon = arrow(transpose([[logebv[where(iiqual_non and iiwac)]],[logebv[where(iiqual_non and iiwac)]]]),transpose([[nhnon[where(iiqual_non and iiwac)]],[nhnon[where(iiqual_non and iiwac)]+0.1]]),color="orange",fill_transparency=75,target=[pw],/data,head_size=0.3,clip=1)
        pwnhdet = plot(logebv[where(iiqual_det and iiwac)],nhdet[where(iiqual_det and iiwac)],"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov,name=" X-ray detected")
        p = plot(logebv[where(iiqual_det and iiwac)],nhdet[where(iiqual_det and iiwac)],"S",sym_filled=0,sym_size=1.5,transparency=85,/ov)
        ;; error bars
        perr = plot([-2.2],[enh.yra[1]-(enh.yra[1]-nh_bound[1])/2.],'s',sym_size=4,sym_thick=2,/ov)
        perr = errorplot([-2.2],[enh.yra[1]-(enh.yra[1]-nh_bound[1])/2.],[mean(e_nhdet[where(iiqual_det)])],errorbar_capsize=0.1,linestyle="",/ov)
        perr = plot([-2.2],[enh.yra[1]-(enh.yra[1]-nh_bound[1])/2.],"S",sym_filled=1,sym_size=1.,col='white',/ov)
        perr = plot([-2.2],[enh.yra[1]-(enh.yra[1]-nh_bound[1])/2.],"S",sym_filled=0,sym_size=1.,col=[12,123,220],/ov)
        ;; text
        t = text(enh.xra[1]-0.2,nh_bound[1]+0.12,'CT',col='black',target=pwnhnon,alignment=1.,font_size=14,font_name='Times',font_style='Bold',/data)

        ;; SECONDARY SOURCES
        pr = plot(logebv[where(iiqual_non and ~iiwac)],nhnon[where(iiqual_non and ~iiwac)],_extra=enh,current=1,/device,position=[dim[0]/2,mrg[1],(dim-mrg)[0],dim[1]/2],/nodata)
        pr.axes[1].showtext=0
        ;; shading
        pr = plot(enh.xra,nh_bound[1]*[1,1],linestyle='',fill_background=1,fill_level=enh.yra[1],fill_color='light grey',fill_transparency=0,/ov)
        ;; plot data
        prnhnon = arrow(transpose([[logebv[where(iiqual_non and ~iiwac)]],[logebv[where(iiqual_non and ~iiwac)]]]),transpose([[nhnon[where(iiqual_non and ~iiwac)]],[nhnon[where(iiqual_non and ~iiwac)]+0.1]]),color="orange",fill_transparency=75,target=[pr],/data,head_size=0.3,clip=1)
        prnhdet = plot(logebv[where(iiqual_det and ~iiwac)],nhdet[where(iiqual_det and ~iiwac)],"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov,name=" X-ray detected")
        p = plot(logebv[where(iiqual_det and ~iiwac)],nhdet[where(iiqual_det and ~iiwac)],"S",sym_filled=1,sym_size=1.5,transparency=85,/ov)
        prnhnon = plot(logebv[where(iiqual_non and ~iiwac)],nhnon[where(iiqual_non and ~iiwac)],"td",sym_filled=1,sym_size=1.,color="orange",/current,/nodata,position=[dim[0]/2,mrg[1],(dim-mrg)[0],dim[1]/2],_extra=ell,axis_style=0)
        a_rnhnon = axis('y',target=pr,location=[enh.xra[1],0,0],textpos=1,tickdir=1,title='$log !8N!7_{H,lim} [cm^{-2}]$',tickformat='(i0)',tickfont_name='Times',tickfont_size=16)    

        ;; legend
        l = legend(target=[prlldet,prllnon],position=[0.03,0.25],/relative,sample_width=0.,horizontal_spacing=0.06,horizontal_alignment=0.,font_name='Times')
        ;; axis
        xt = text(0.5,0.026,'$log !8E!7(!8B-V!7)$',alignment=0.5,font_size=14,font_name='Times')

        if keyword_set(sav) then begin
            print, '    SAVING PLOT'
            if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/'+savfile[i]+'.eps',/BITMAP else $
                                                         p.save,'figures/'+savfile[i]+'.png',resolution=res
        endif
    endfor    
endif



;;----------------------------------------------------------------------------------------
;; NH Distribution
;;----------------------------------------------------------------------------------------
if keyword_set(nh_dist) then begin
    print, '    PREPARING NH DISTRIBUTION'
    
    wpnorm = 0
    rpnorm = 0   
    wbnorm = 0
    rbnorm = 0    
    
    ;; BORUS source numbers
    ctwd_borus = commas(fix(total(yhist_wdet_borus)))
    ctwn_borus = commas(fix(total(yhist_wnon_borus)))
    ctrd_borus = commas(fix(total(yhist_rdet_borus)))
    ctrn_borus = commas(fix(total(yhist_rnon_borus)))
    ;; POWER LAW source numbers
    ctwd_power = commas(fix(total(yhist_wdet_power)))
    ctwn_power = commas(fix(total(yhist_wnon_power)))
    ctrd_power = commas(fix(total(yhist_rdet_power)))
    ctrn_power = commas(fix(total(yhist_rnon_power)))
   
    ;; normalize or not
    if (wpnorm eq 1) then begin
        ywnp = nm(yhist_wnon_power)
        ywdp = nm(yhist_wdet_power)
    endif else begin
        ywnp = yhist_wnon_power
        ywdp = yhist_wdet_power
    endelse
    if (rpnorm eq 1) then begin
        yrnp = nm(yhist_rnon_power)
        yrdp = nm(yhist_rdet_power)
    endif else begin
        yrnp = yhist_rnon_power
        yrdp = yhist_rdet_power
    endelse
    if (wbnorm eq 1) then begin
        ywnb = nm(yhist_wnon_borus)
        ywdb = nm(yhist_wdet_borus)
    endif else begin
        ywnb = yhist_wnon_borus
        ywdb = yhist_wdet_borus
    endelse
    if (rbnorm eq 1) then begin
        yrnb = nm(yhist_rnon_borus)
        yrdb = nm(yhist_rdet_borus)
    endif else begin
        yrnb = yhist_rnon_borus
        yrdb = yhist_rdet_borus
    endelse

    ;; lower limit arrows
    ;; arrow x/y span
    yarr = [1.:5.]
    yarr = transpose([[yarr/10.],[yarr/10.]])
    ;if (normalize eq 0) then yarr*=ew.yra[1]
    ;xarr = transpose([[intarr(5)+alog10(1.5e24)+0.1],[intarr(5)+alog10(1.5e24)+0.3]])
    xarr = transpose([[intarr(n_elements(yarr[0,*]))-0.1],[intarr(n_elements(yarr[0,*]))+0.1]])
    ;print, yarr
    
    ;; WISE AGN source detections 
    e = {xra:[20.8,25.2], $
         stairstep:1,fill_background:1, $
         font_name:'Times',font_size:16, $
         buffer:0}
    dim = [800,700]
    mrg = [80,65]
    if keyword_set(hide) then ew.buffer = 1    
    wpyra = [0.,ceil((max(ywdp)>max(ywnp))/100.)*110.]
    rpyra = [0.,ceil((max(yrdp)>max(yrnp))/100.)*100.]
    wbyra = [0.,ceil((max(ywdb)>max(ywnb))/100.)*100.]
    rbyra = [0.,ceil((max(yrdb)>max(yrnb))/100.)*100.]
        
    ;; POWER-LAW MODEL -- WISE AGN
    pw = plot(xhist_wnon_power,ywnp,_extra=e,yra=wpyra,dim=dim,position=[mrg[0],dim[1]/2.,dim[0]/2.,(dim-mrg)[1]],/device,/nodata)
    pw.axes[0].showtext=0
    ;; shading
    p = plot([alog10(1.5e24),e.xra[1]],pw.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=pw.yra[0],fill_color='light grey',fill_transparency=0,/ov)
    ;; data
    hwnon = plot(xhist_wnon_power,ywnp,_extra=e,col=[255,194,10],thick=2,linestyle='__',fill_color='orange',fill_transparency=75,/ov)
    re = gaussfit(xhist_wnon_power,ywnp,a,nterms=3)
    hwlim = arrow(xarr+a[1],yarr*pw.yra[1],color='orange',/ov,/data,thick=4,fill_transparency=75,head_size=0.8,target=hwnon)
    hwdet = plot(xhist_wdet_power,ywdp,_extra=e,col=[12,123,220],thick=2,fill_color='dodger blue',fill_transparency=75,/ov)
    ;; text
    t = text(0.5,0.90,'Power Law',target=pw,/relative,alignment=0.5,font_size=14,font_name='Times',font_style='Bold')
    t = text(0.05,0.90,'!16WISE!15 AGN',col='black',target=pw,/relative,font_size=14,font_name='Times',font_style='Bold')
    t = text(0.25,0.70,'X-ray det:',col=[12,123,220],target=pw,/relative,alignment=0.5,font_size=14,font_name='Times',font_style='Bold')
    t = text(0.25,0.65,strtrim(ctwd_borus,2)+' sources',col=[12,123,220],target=pw,/relative,alignment=0.5,font_size=14,font_name='Times',font_style='Bold')
    t = text(0.25,0.55,'X-ray non-det:',col=[255,194,10],target=pw,/relative,alignment=0.5,font_size=14,font_name='Times',font_style='Bold')
    t = text(0.25,0.50,strtrim(ctwn_borus,2)+' sources',col=[255,194,10],target=pw,/relative,alignment=0.5,font_size=14,font_name='Times',font_style='Bold')

    ;; POWER-LAW MODEL -- Secondary SOURCES
    pr = plot(xhist_rnon_power,yhist_rnon_power,_extra=e,yra=rpyra,current=1,position=[dim[0]/2.,dim[1]/2.,(dim-mrg)[0],(dim-mrg)[1]],/device,/nodata)
    pr.axes[1].showtext=0
    pr.axes[0].showtext=0
    ;; shading
    p = plot([alog10(1.5e24),pr.xra[1]],pr.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=pr.yra[0],fill_color='light grey',fill_transparency=0,/ov)
    ;; data
    hrnon = plot(xhist_rnon_power,yrnp,_extra=e,col=[255,194,10],thick=2,linestyle='__',fill_color='orange',fill_transparency=75,/ov)
    a_rnon = axis('y',target=hrnon,location=[pr.xra[1],0,0],textpos=1,tickdir=1,tickfont_name='Times',tickfont_size=14)    
    re = gaussfit(xhist_rnon_power,yrnp,a,nterms=3)
    hrlim = arrow(xarr+a[1],yarr*pr.yra[1],color='orange',/ov,/data,thick=4,fill_transparency=75,head_size=0.8,target=hrnon)
    hrdet = plot(xhist_rdet_power,yrdp,_extra=e,col=[12,123,220],thick=2,fill_color='dodger blue',fill_transparency=75,/ov)
    ;; text
    t = text(0.05,0.90,'Secondary',col='black',target=pr,/relative,font_size=14,font_name='Times',font_style='Bold')
    t = text(0.25,0.70,'X-ray det:',col=[12,123,220],target=pr,/relative,alignment=0.5,font_size=14,font_name='Times',font_style='Bold')
    t = text(0.25,0.65,strtrim(ctrd_borus,2)+' sources',col=[12,123,220],target=pr,/relative,alignment=0.5,font_size=14,font_name='Times',font_style='Bold')
    t = text(0.25,0.55,'X-ray non-det:',col=[255,194,10],target=pr,/relative,alignment=0.5,font_size=14,font_name='Times',font_style='Bold')
    t = text(0.25,0.50,strtrim(ctrn_borus,2)+' sources',col=[255,194,10],target=pr,/relative,alignment=0.5,font_size=14,font_name='Times',font_style='Bold')
    t = text(0.815,0.90,'CT',col='black',target=pr,font_size=14,font_name='Times',font_style='Bold',/relative)          
    
    ;; BORUS MODEL -- WISE AGN
    pw = plot(xhist_wnon_borus,ywnb,_extra=e,yra=wbyra,current=1,position=[mrg[0],mrg[1],dim[0]/2.,dim[1]/2.],/device,/nodata)
    pw.ytickvalues = pw.ytickvalues[0:-2]
    ;; shading
    p = plot([alog10(1.5e24),pw.xra[1]],pw.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=pw.yra[0],fill_color='light grey',fill_transparency=0,/ov)
    ;; data
    hwnon = plot(xhist_wnon_borus,ywnb,_extra=e,col=[255,194,10],thick=2,linestyle='__',fill_color='orange',fill_transparency=75,/ov)
    re = gaussfit(xhist_wnon_borus,ywnb,a,nterms=3)
    hwlim = arrow(xarr+a[1],yarr*pw.yra[1],color='orange',/ov,/data,thick=4,fill_transparency=75,head_size=0.8,target=hwnon)
    hwdet = plot(xhist_wdet_borus,ywdb,_extra=e,col=[12,123,220],thick=2,fill_color='dodger blue',fill_transparency=75,/ov)
    ;; text
    xmod = text(pw.xra[0]+diff(pw.xra)/2.,pw.yra[1]*0.89,'BORUS',target=pw,/data,alignment=0.5,font_size=14,font_name='Times',font_style='Bold')
    
    ;; BORUS MODEL -- Secondary SOURCES
    pr = plot(xhist_rnon_borus,yrnb,_extra=e,current=1,position=[dim[0]/2.,mrg[1],(dim-mrg)[0],dim[1]/2.],/device,/nodata)
    pr.axes[1].showtext=0
    ;pr.xtickvalues = pr.xtickvalues[1:-1]
    ;; shading
    p = plot([alog10(1.5e24),pr.xra[1]],pr.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=pr.yra[0],fill_color='light grey',fill_transparency=0,/ov)
    ;; data
    hrnon = plot(xhist_rnon_borus,yrnb,_extra=e,yra=rbyra,col=[255,194,10],thick=2,linestyle='__',fill_color='orange',fill_transparency=75,/ov,name=' X-ray non-det.')
    a_rnon = axis('y',target=hrnon,location=[pr.xra[1],0,0],textpos=1,tickdir=1,tickfont_name='Times',tickfont_size=14)    
    a_rnon.tickvalues = a_rnon.tickvalues[0:-2]
    re = gaussfit(xhist_rnon_borus,yrnb,a,nterms=3)
    hrlim = arrow(xarr+a[1],yarr*pr.yra[1],color='orange',/ov,/data,thick=4,fill_transparency=75,head_size=0.8,target=hrnon)
    hrdet = plot(xhist_rdet_borus,yrdb,_extra=e,col=[12,123,220],thick=2,fill_color='dodger blue',fill_transparency=75,/ov,name=' X-ray detected')

    ;; axes
    xt = text(0.5,0.035,'$log !8N!7_{H} [cm^{-2}]$',alignment=0.5,vertical_alignment=0.5,font_size=16, font_name='Times',/relative)            
    yt = text(0.030,0.5,'Frequency (!8WISE!7 AGN)',orientation=90.,alignment=0.5,vertical_alignment=0.5,font_size=16,font_name='Times')
    yt2 = text(0.970,0.5,'Frequency (Secondary)',orientation=90.,alignment=0.5,vertical_alignment=0.5,font_size=16,font_name='Times')
    ;; legend
    l = legend(target=[hrdet,hrnon],position=[0.510,0.48],/normal,horizontal_alignment=0.,sample_width=0.12,font_size=14,font_name='Times')
    
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/nh_dist.eps',/BITMAP else $
                                                     p.save,'figures/nh_dist.png',resolution=res
    endif
endif


;;----------------------------------------------------------------------------------------
;; Survival Analysis
;;----------------------------------------------------------------------------------------
if keyword_set(surv_anal) then begin
    print, '    PREPARING SURVIVAL ANALYSIS RESULTS'

    ;col = [[204,121,167],[213,94,0],[0,158,115],[0,114,178],[240,228,66]]
    rem_col = [65,182,196];[49,54,149]
    wise_col = [215,48,39]
   
    ;rem_col = [77,146,33]
    ;wise_col = [197,27,125]

    ;; KM ESTIMATOR
    ekm = {xra:[-3.,2.],yra:[0.,1.05], $
           sym_size:1.,sym_filled:0, $
           stairstep:1, $
           xtitle:'$!8R!7_{!8L!7_X}$',ytitle:'KM Estimator', $
           font_name:'Times',font_size:16, $
           buffer:0}
    if keyword_set(hide) then ekm.buffer = 1    
    p = plot(-abin,akm,'-',_extra=e,/nodata)
    pr = plot(-rbin,e_rkmhi,'-',col=rem_col,_extra=ekm,fill_background=1,fill_color=rem_col,fill_level=ekm.yra[0],fill_transparency=75,/ov)
    pr = plot(-rbin,e_rkmlo,'-',col=rem_col,_extra=ekm,fill_background=1,fill_color='white',fill_level=ekm.yra[0],/ov)
    pr = plot(-rbin,rkm,'--',col=rem_col,_extra=ekm,thick=2,/ov,name='Secondary')
    ;; white square to remove the extended lines
    pclear = plot([ekm.xra[0],min(-rbin)-0.015],[1.,1.],col='white',fill_background=1,fill_color='white',fill_level=ekm.yra[0],/ov)
    pw = plot(-wbin,e_wkmhi,'-',col=wise_col,_extra=ekm,fill_background=1,fill_color=wise_col,fill_level=ekm.yra[0],fill_transparency=75,/ov)
    pw = plot(-wbin,e_wkmlo,'-',col=wise_col,_extra=ekm,fill_background=1,fill_color='white',fill_level=ekm.yra[0],/ov)
    pw = plot(-wbin,wkm,'-',col=wise_col,_extra=ekm,thick=2,/ov,name='!8WISE!7 AGN')
    ;; white square to remove the extended lines
    pclear = plot([ekm.xra[0],min(-wbin)-0.054],[0.6,0.6],col='white',fill_background=1,fill_color='white',fill_level=ekm.yra[0],/ov)
    ;; legend
    l = legend(target=[pw,pr],position=[0.86,0.21],/relative,sample_width=0.15,horizontal_spacing=0.06,font_name='Times',font_size=14)
    
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/km_est.eps',/BITMAP else $
                                                     p.save,'figures/km_est.png',resolution=res
    endif

    ;; SURVIVAL ANALYSIS SIMULATED NH DISTRIBUTION
    enh = {xra:[20.6,25.6],yra:[0.,ceil(max(ywnh_borus)/100.)*100.], $
           stairstep:1,fill_background:1,thick:2, $
           xtitle:'$log !8N!7_{H} (simulated) [cm^{-2}]$',ytitle:'Frequency', $
           font_name:'Times',font_size:16, $
           buffer:0}
    if keyword_set(hide) then enh.buffer = 1    
    
    p = plot(xrnh_borus,yrnh_borus,_extra=enh,/nodata)
    p = plot([alog10(1.5e24),enh.xra[1]],enh.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=enh.yra[0],fill_color='light grey',fill_transparency=0,/ov)
    ;pr = plot(xrnh_borus,yrnh_borus,'-.',col=rem_col,_extra=enh,thick=2,fill_color=rem_col,fill_transparency=75,/ov,name='Secondary')
    pw = plot(xwnh_borus,ywnh_borus,'-',col=wise_col,_extra=enh,thick=2,fill_color=wise_col,fill_transparency=75,/ov,name='!8WISE!7 AGN')
    maxnh = (rl2nh(min(lldet[where(iiqual_det)]),model='BORUS'))[0]
    ibin = (where(xwnh_borus le maxnh))[-1]
    pdet = plot([xwnh_borus[ibin]],[ywnh_borus[ibin]],'|',sym_thick=2,sym_size=2,col=wise_col,/ov)
    ;; upper limit arrows based on Secondary bins    
    yarr = [1.:5.]
    yarr = transpose([[yarr/10.],[yarr/10.]])
    xarr = transpose([[intarr(n_elements(yarr[0,*]))-0.1],[intarr(n_elements(yarr[0,*]))+0.1]])
    hrlim = arrow(xarr+xwnh_borus[ibin+1],yarr*enh.yra[1],color=wise_col,/ov,/data,thick=4,fill_transparency=75,head_size=0.8,target=hrnon)
    ;; text
    ct = text(0.74,0.90,'CT',col='black',target=p,font_size=14,font_name='Times',font_style='Bold',/relative,fill_background=1,fill_color='light grey')          
    ;; legend
    l = legend(target=[pw],position=[0.32,0.94],/relative,sample_width=0.15,horizontal_spacing=0.06,font_name='Times',font_size=14)

    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/sim_nh.eps',/BITMAP else $
                                                     p.save,'figures/sim_nh.png',resolution=res
    endif    
endif


;;----------------------------------------------------------------------------------------
;; X-ray Stacking
;;----------------------------------------------------------------------------------------
if keyword_set(xstack) then begin
    print, '    PREPARING X-RAY STACKING RESULTS'

    file = file_search('stack_output/'+['det_wagn','det_ragn','nodet_wagn','nodet_ragn']+'_057_stack.png')
    if (file[0] eq '') then STOP
    
    ;; index shorthand
    group = ['WD','RD','WN','RN']
    ngroups = n_elements(group)
    iiwd = iiqual_det and xdet eq 'CHA' and iiwac
    iird = iiqual_det and xdet eq 'CHA' and ~iiwac
    iiwn = iiqual_non and xnon eq 'CHA' and iiwac
    iirn = iiqual_non and xnon eq 'CHA' and ~iiwac
    ;; x-ray stack flux estimates
    ;wdfx = 2.7119369e-13
    ;wdfx_str = '(9.53\pm2.53)\times10^{-14}'
    ;rdfx = 1.0583464e-13
    ;rdfx_str = '(2.43\pm0.57)\times10^{-14}'
    ;wnfx = 1.18e-15
    ;wnfx_str = '(1.18\pm0.31)\times10^{-15}'
    ;rnfx = 2.56e-16 
    ;rnfx_str = '(2.56\pm0.24)\times10^{-16}'

    dim = 800
    edg = 80
    pos = [[edg,dim/2,dim/2,dim-edg], $
           [dim/2,dim/2,dim-edg,dim-edg], $
           [edg,edg,dim/2,dim/2], $ 
           [dim/2,edg,dim-edg,dim/2]]
    
    off = [[0,0,0,0], $
           [-1,0,-1,0], $
           [0,1,0,1], $
           [-1,1,-1,1]]
    pos = pos+off
    xy = findgen(35,start=-17)
    ra = minmax(xy)

    e = {xra:ra,yra:ra, $
         font_name:'Times',font_size:16, $
         dim:[dim,dim],device:1, $
         buffer:0}
    if keyword_set(hide) then e.buffer = 1
             
    p = plot(xy,xy,/nodata,_extra=e,current=0,position=pos[*,0],xshowtext=0)
    p = plot(xy,xy,/nodata,_extra=e,current=1,position=pos[*,1],xshowtext=0,yshowtext=0)
    p = plot(xy,xy,/nodata,_extra=e,current=1,position=pos[*,2])
    p = plot(xy,xy,/nodata,_extra=e,current=1,position=pos[*,3],yshowtext=0)
    im1 = image(file[0],current=1,position=pos[*,0],/device)    
    im2 = image(file[1],current=1,position=pos[*,1],/device)
    im3 = image(file[2],current=1,position=pos[*,2],/device)
    im4 = image(file[3],current=1,position=pos[*,3],/device)    
    ;; catalog tags
    wiset = text(0.12,0.862,'!16WISE!15 AGN',/relative,font_name='Times',font_size=14,font_style='Bold')
    resit = text(0.52,0.862,'Secondary',/relative,font_name='Times',font_size=14,font_style='Bold')
    ;; detection tags
    dett = text(0.783,0.862,'X-ray det.',/relative,font_size=14,font_name='Times',font_style='Bold')
    nont = text(0.737,0.46,'X-ray non-det.',/relative,font_size=14,font_name='Times',font_style='Bold')
    ;; energy band
    engt = text(0.382,0.862,'0.5-7 keV',/relative,alignment=0.,font_size=14,font_name='Times',font_style='Bold')
    ;;; number of stacked sources
    ;npos = [[0.48,0.86],[0.88,0.86], $
    ;        [0.48,0.46],[0.88,0.46]]
    ;fpos = [[0.12,0.55],[0.52,0.55], $
    ;        [0.12,0.15],[0.52,0.15]]
    ;lpos = fpos
    ;lpos[1,*] -= 0.03
    ;zpos = lpos
    ;zpos[0,*] += 0.15
    ;for i = 0,ngroups-1 do begin
    ;    re = execute('ind = II'+group[i])
    ;    numt = text(npos[0,i],npos[1,i],'!16N!15$_{src}$ ('+commas(string(total(ind),format='(i0.0)'))+')',/relative,alignment=1.,font_name='Times',font_size=14,font_style='Bold')
    ;    re = execute('flx = '+group[i]+'fx')
    ;    re = execute('flx_str = '+group[i]+'fx_str')
    ;    flxt = text(fpos[0,i],fpos[1,i],'$'+flx_str+' erg s^{-1} cm^{-2}$',/relative,font_name='Times',font_size=14,font_style='Bold')
    ;    ;flxt = text(fpos[0,i],fpos[1,i],'$'+string(rnd((strsplit(flx,'e',/extract))[0],2),format='(d4.2)')+'x10^{'+(strsplit(flx,'e',/extract))[1]+'} erg s^{-1} cm^{-2}$',/relative,font_name='Times',font_size=14,font_style='Bold')
    ;    zz = z[where(ind)]
    ;    dl = mean(dlum(zz,/sq))
    ;    logl6 = alog10(mean(4.*!const.pi*dl*flx))
    ;    lumt = text(lpos[0,i],lpos[1,i],'log $\langle !16L!15_X\rangle=$'+string(rnd(logl6,2),format='(d5.2)'),/relative,font_size=14,font_name='Times',font_style='Bold')
    ;endfor

    xt = text(0.5,0.035,'offset in Dec. [arcsec.]',alignment=0.5,font_size=16,font_name='Times')
    yt = text(0.035,0.5,'offset in RA [arcsec.]',alignment=0.5,orientation=90.,font_size=16,font_name='Times')

    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/xstack.eps',/BITMAP else $
                                                     p.save,'figures/xstack.png',resolution=res
    endif
endif


if keyword_set(gstack) then begin
    print, '    PREPARING GALAXY STACKING RESULTS'

    soft = file_search('stack_output/*_052_stack.png')
    hard = file_Search('stack_output/*_27_stack.png')
    if (soft[0] eq '') then begin
        STOP
    endif

    dim = [800,1120]
    mrg = [80,80]
    psp = (dim[1]-2.*mrg[1])/3.
    pos = [[mrg[0],mrg[1]+psp*2,dim[0]/2,(dim-mrg)[1]], [dim[0]/2,mrg[1]+psp*2,(dim-mrg)[0],(dim-mrg)[1]], $
           [mrg[0],mrg[1]+psp,dim[0]/2,mrg[1]+psp*2],   [dim[0]/2,mrg[1]+psp,(dim-mrg)[0],mrg[1]+psp*2], $
           [mrg[0],mrg[1],dim[0]/2,mrg[1]+psp],         [dim[0]/2,mrg[1],(dim-mrg)[0],mrg[1]+psp]]       

    off = [[0,0,0,0], $
           [-1,0,-1,0], $
           [0,1,0,1], $
           [-1,1,-1,1], $
           [0,2,0,2], $
           [-1,2,-1,2]]
    pos = pos+off
    xy = findgen(35,start=-17)
    ra = minmax(xy)

    e = {xra:ra,yra:ra, $
         font_name:'Times',font_size:16, $
         dim:dim,device:1, $
         buffer:0}
    if keyword_set(hide) then e.buffer = 1

    p = plot(xy,xy,/nodata,_extra=e,current=0,position=pos[*,0],xshowtext=0)
    p = plot(xy,xy,/nodata,_extra=e,current=1,position=pos[*,1],xshowtext=0,yshowtext=0)
    p = plot(xy,xy,/nodata,_extra=e,current=1,position=pos[*,2],xshowtext=0)
    p = plot(xy,xy,/nodata,_extra=e,current=1,position=pos[*,3],xshowtext=0,yshowtext=0)
    p = plot(xy,xy,/nodata,_extra=e,current=1,position=pos[*,4])
    p = plot(xy,xy,/nodata,_extra=e,current=1,position=pos[*,5],yshowtext=0)
    im1 = image(soft[2],current=1,position=pos[*,0],/device)    
    im2 = image(hard[2],current=1,position=pos[*,1],/device)
    im3 = image(soft[0],current=1,position=pos[*,2],/device)
    im4 = image(hard[0],current=1,position=pos[*,3],/device)    
    im5 = image(soft[1],current=1,position=pos[*,4],/device)
    im6 = image(hard[1],current=1,position=pos[*,5],/device)    

    ;; catalog tags
    scdt = text(0.12,0.902,'Secondary',/relative,font_size=14,font_name='Times',font_style='Bold')
    remt = text(0.12,0.616,'Removed AGNs',/relative,font_size=14,font_name='Times',font_style='Bold')
    galt = text(0.12,0.332,'SED Galaxies',/relative,font_size=14,font_name='Times',font_style='Bold')
    ;; energy range
    softt = text(0.382,0.902,'0.5-2 keV',/relative,alignment=0.,font_size=14,font_name='Times',font_style='Bold')
    hardt = text(0.800,0.902,'2-7 keV',/relative,alignment=0.,font_size=14,font_name='Times',font_style='Bold')
        
    xt = text(0.5,0.028,'offset in Dec. [arcsec.]',alignment=0.5,font_size=16,font_name='Times')
    yt = text(0.035,0.5,'offset in RA [arcsec.]',alignment=0.5,orientation=90.,font_size=16,font_name='Times')

    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/gstack.eps',/BITMAP else $
                                                     p.save,'figures/gstack.png',resolution=res
    endif
endif



END












