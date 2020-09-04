PRO plot_xray_lack_agn, PHOT_SPEC = phot_spec, $
                        SEDS = seds, $
                        LX_LIR = lx_lir, $
                        FLUX_LIMIT = flux_limit, $
                        LUM_RATIO = lum_ratio, $
                        NH_DIST = nh_dist, $
                        SHOW = show, $
                        REDO = redo, $
                        SAV = sav, $
                        ALL = all       


;; load data
common _fits    
common _resamp
common _inf_cha 
common _inf_xmm 
common _inf_nst 
common _det_cha 
common _det_xmm 
common _det_nst 
common _det_wac 
common _xconv   
common _fxlim 
common _comp    
common _agnlum 
common _clean_cha
common _clean_xmm
common _clean_nst
common _quality  
common _combined
common _nhdist


file_mkdir,'figures'
;;----------------------------------------------------------------------------------------
;; SEDs
;;----------------------------------------------------------------------------------------
if keyword_set(phot_spec) then begin
    print, '    PREPARING REDSHIFT COMPARISON'
    
    file = file_search()
    ifile = where(strmatch(file,'phot_spec_srcs.sav'),nfile)
    if keyword_set(redo) then nfile = 0
    if (nfile eq 0) then begin
        ;; zCOSMOS
        ;dir_cosmos = '/Users/ccarroll/Research/surveys/COSMOS/zCOSMOS/cesam_zcosbrightspec20k_dr3_catalog_1513358585.fits'
        ;cosm = mrdfits(dir_cosmos,1)
        ;; SDSS DR12
        ;dir_sdss = '/Users/ccarroll/Research/surveys/SDSS/DR14/sdss-dr14-cat-part26.fits.gz'
        ;sdss = mrdfits(dir_sdss,1)
        ;; AGES
        dir_ages = '/Users/ccarroll/Research/surveys/NDWFS/AGES-survey-J_ApJS_200_8_sources.fits'
        ages = mrdfits(dir_ages,1)
        ;; SDSS DR12
        dir_sdss = '/Users/ccarroll/Research/surveys/SDSS/DR14/sdss-dr14-cat-part'+['37','38']+'.fits.gz'
        sdss = [mrdfits(dir_sdss[0],1),mrdfits(dir_sdss[1],1)]
        ;; quality sources
        ;cosm = cosm[where(cosm.zpec gt 0.,/null)] 
        ages = ages[where((ages.gal eq 1 or ages.qso eq 1 or ages.agn eq 1) and ages.z1 gt 0. and ages.s_n1 gt 3.,/null)]
        sdss = sdss[where(sdss.zp gt 0. and sdss.zp le 0.6,/null)]
        sdss = sdss[where(sdss.photoerrorclass ge -1 and sdss.photoerrorclass le 3,/null)]
    
        ;; matched to SDSS DR12
        ;spherematch,sdss.ra_sdss,sdss.dec_sdss,cosm.ra,cosm.dec,6./3600.,isdss,icosm,sep
        ;cosm = cosm[icosm]
        ;sdss = sdss[isdss]
        spherematch,sdss.ra_sdss,sdss.dec_sdss,ages._raj2000,ages._dej2000,6./3600.,isdss,iages,sep_sdss
        ages = ages[iages]
        sdss = sdss[isdss]
        sep_sdss *= 3600.
        
        isep = where(sep_sdss lt 1.,ct)
        if (ct gt 0.) then begin
            ages = ages[isep]
            sdss = sdss[isep]
            sep_sdss = sep_sdss[isep]
        endif
        
        save,ages,sdss,sep_sdss,/compress,file='phot_spec_srcs.sav'
    endif else restore, file[ifile]
    
    ;zs = cosm.zpec
    isep = where(sep_sdss lt 1.0,ct)
    if (ct gt 0.) then begin
        ages = ages[isep]
        sdss = sdss[isep]
        sep_sdss = sep_sdss[isep]
    endif    
    zs = ages.z1
    zp = sdss.zp
    ;; all matches
    delz = (zp-zs)/(1+zs)
    mad = medabsdev(delz)
    ;; matches in final sample
    ifin = where(iifinal and strmatch(ztype,'ZP'),/null)
    ra_fin = ra[ifin]
    dec_fin = dec[ifin]
    zp_fin = z[ifin]
    ;spherematch,ra[iagn],dec[iagn],cosm.ra,cosm.dec,6./3600.,isamp,icosm,sep
    spherematch,ra_fin,dec_fin,ages._raj2000,ages._dej2000,6./3600.,isamp,iages,sep_samp
    zs_in = ages[iages].z1
    zp_in = zp_fin[isamp]
    delz_in = (zp_in-zs_in)/(1+zs_in)
    mad_in = medabsdev(delz_in)
        
    e = {xra:[0.,1.0],yra:[0.,0.6],aspect_ratio:1, $
         sym_size:0.5,sym_filled:1,transparency:75, $
         xtitle:'$!8z!7_{spec} (VLT)$',ytitle:'$!8z!7_{phot} (SDSS)$', $
         font_name:'Times',font_size:14, $
         dimension:[700,800],buffer:1}
    if keyword_set(show) then e.buffer = 0
    col = [0,114,178]
    col_in = [[204,121,167],[213,94,0],[0,158,115],[0,114,178],[240,228,66]]
    p = plot(zs,zp,'o',color=col,_extra=e)
    p.position = [0.12735492,0.33393553,0.769520,0.83583007]
    p = plot([0.,1.],[0.,1.],'--',thick=2,/ov)
    p = plot(zs_in,zp_in,'o',color=col_in[*,4],sym_size=0.5,sym_filled=1,/ov)
    p = plot(zs_in,zp_in,'o',color=black,sym_size=0.75,sym_thick=1.,sym_filled=0,transparency=50,/ov)
    ioff = where(zs gt e.xra[1] and zp lt e.yra[1]-0.02,noff)
    if (noff gt 0) then poff = arrow(transpose([[make_array(noff,value=e.xra[1]-0.012)],[make_array(noff,value=e.xra[1])-0.002]]),transpose([[zp[ioff]],[zp[ioff]]]),color=col,/ov,/data,fill_transparency=100,head_size=0.8,clip=1,target=p)
    p.axes[0].showtext = 0
    p.ytickvalues = (p.ytickvalues)[0:-2]
    ;; residuals plot
    er = {xra:e.xra,yra:[-0.8,0.8],aspect_ratio:0, $
          sym_size:0.5,sym_filled:1,transparency:75, $
          xtitle:'$!8z!7_{spec} (VLT)$',ytitle:'$\Delta!8z!7 / (1+!8z!7_{spec})$', $
          font_name:'Times',font_size:14, $
          dimension:[700,800],buffer:1}
    pr = plot(zs,delz,'o',color=col,_extra=er,/current)
    pr.position = [p.position[0],0.13195312,p.position[2],p.position[1]]
    pr = plot(p.xra,[0,0],'--',thick=2,/ov)
    pr = plot(zs_in,delz_in,'o',color=col_in[*,4],sym_size=0.5,sym_filled=1,/ov)
    pr = plot(zs_in,delz_in,'o',color=black,sym_size=0.75,sym_thick=1.,sym_filled=0,transparency=60,/ov)
    iroff = where(zs gt er.xra[1] and delz gt er.yra[0],nroff)
    if (nroff gt 0) then proff = arrow(transpose([[make_array(nroff,value=e.xra[1]-0.012)],[make_array(nroff,value=e.xra[1])-0.002]]),transpose([[delz[iroff]],[delz[iroff]]]),color=col,/ov,/data,fill_transparency=100,head_size=0.8,clip=1,target=pr)
    if (er.yra[1] mod 2 ne 0) then pr.xtickvalues = (pr.xtickvalues)[0:-2]
    pr.ytickvalues = [-0.4,0.,0.4]
    pr.yminor = 3.
    t = text(target=pr,er.xra[1]*0.7,er.yra[1]*0.6,'$\sigma_{MAD} = '+string(rnd(mad,3),format='(d5.3)')+'$',/data,font_size=14,font_style='Bold',font_name='Times',fill_background=1,fill_color='white')
    ;; histogram of residuals
    yh = histogram(delz,bin=scott(delz),locations=xh)
    eh = {xra:[0.,ceil(max(yh)/95.)*100.],yra:[-0.3,0.3], $
          stairstep:1,fill_background:1,fill_color:[0,114,178],fill_transparency:10, $
          xtitle:'Frequency', $
          font_name:'Times',font_size:14}
    ph = plot(yh,xh,/current,_extra=eh)
    ph = plot(ph.xra,[0.,0.],'--',thick=2,/ov)
    ph.position = [p.position[2],pr.position[1],0.91303385,pr.position[3]]
    ph.axes[1].showtext = 0
    ph.xtickvalues = [0.:ph.xra[1]:ph.xra[1]/2.]
    ph.xminor = 3.
    ;; histogram of zs
    yh_zs = histogram(zs,bin=scott(zs),locations=xh_zs)
    eh_zs = {xra:e.xra,yra:eh.xra, $
            stairstep:1,fill_background:1,fill_color:[0,114,178],fill_transparency:10, $
            ytitle:'Frequency', $
            font_name:'Times',font_size:14}
    ph_zs = plot(xh_zs,yh_zs,/current,_extra=eh_zs)
    ph_zs.position = [p.position[0],p.position[3],p.position[2],0.95303385]
    ph_zs.axes[0].showtext = 0
    ph_zs.ytickvalues = ph.xtickvalues
    ph_zs.yminor = 3.
    ;; histogram of zp
    yh_zp = histogram(zp,bin=scott(zp),locations=xh_zp)
    eh_zp = {xra:eh.xra,yra:e.yra, $
            stairstep:1,fill_background:1,fill_color:[0,114,178],fill_transparency:10, $
            font_name:'Times',font_size:14}
    ph_zp = plot(yh_zp,xh_zp,/current,_extra=eh_zp)
    ph_zp.position = [p.position[2],p.position[1],ph.position[2],p.position[3]]
    ph_zp.axes[1].showtext = 0
    ph_zp.axes[0].showtext = 0
    ph_zp.xtickvalues = ph.xtickvalues
    ph_zp.xminor = 3.

    print, 'NUM. MATCHES: ' + strtrim(n_elements(zs),2)
    print, 'NUM. OFFSETS: ' + strtrim(noff,2)
    print, 'NUM. SAMPLE:   ' + strtrim(n_elements(zs_in),2)   
    print, 'SIGMA MAD:  ' + strtrim(mad,2)
    print, 'SIGMA MAD FINAL: ' + strtrim(mad_in,2)

    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/phot_v_spec.eps',/BITMAP else $
                                                     p.save,'figures/phot_v_spec.png';,resolution=20
    endif
endif


if keyword_set(seds) then begin
    print, '    PREPARING PAPER SED'
    ids = ['1237679254133735725','1237679437739524714','1237664291009921910','1237679321788317973']
    igal = ['1237652900229415056','1237655370892116222','1237653616933535892','1237651272441725119']
    inds = []
    for i = 0,3 do inds = [inds,where(objid eq ids[i])]
    print, iifinal_det[inds]
    print, iifinal_non[inds]
    print, iidet_wac[inds]
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
    nobj = n_elements(inds)
    plot_flux = flux[*,inds]
    plot_err = e_flux[*,inds]
    plot_bin = bin[*,inds]
    plot_id = objid[inds]
    plot_z = z[inds]
    plot_zerr = zerr[inds]
    plot_fits = param[*,inds]
    plot_ebvsig = ebv_sigm[1,inds]
    ;; extract model parameters
    plot_ebv = plot_fits[0,*]
    coeff = plot_fits[2:2+ntemps-1,*]
    plot_chi = plot_fits[-2:-1,*]
    obswav = wave
    objwav = rebin(obswav,n_elements(obswav),nobj)
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
    plot_z = strtrim(string(plot_z,format='(d5.3)'),2)
    plot_zerr = strtrim(string(plot_zerr,format='(d5.3)'),2)
    plot_ebv = strtrim(string(plot_ebv,format='(d5.2)'),2)
    plot_ebvsig = strtrim(string(plot_ebvsig,format='(d4.2)'),2)
    coeff = reform(strtrim(string(coeff,format='(e10.3)'),2),ntemps,nobj)
    plot_chi = strtrim(string(plot_chi[0,*],format='(d0.2)'),2)+' / '+strtrim(string(plot_chi[1,*],format='(i)'),2)
    ;; plot SEDs
    e = {xr:[0.05,30.],yra:[floor(min(plot_flux[where(finite(plot_flux))]))-1.5,ceil(max(plot_flux[where(finite(plot_flux))]))+2.0], $
         xlog:1, $
         font_name:'Times', $
         nodata:1,dimension:[1140,890], $
         buffer:1}
    if keyword_set(show) then e.buffer = 0
    pos = [[80,455,585,835],[585,455,1090,835],[80,75,585,455],[585,75,1090,455]]
    ;label = transpose([['         ID : '+strtrim(plot_id,2)],['$            !8z!7 : $'+plot_z],['$!8E!7(!8B-V!7)$ : '+plot_ebv],['$\chi^2 / DoF$ : '+chi]])
    ;label = transpose([['ID :          '+strtrim(plot_id,2)],['$!8z!7 :             $'+plot_z],['$!8E!7(!8B-V!7)$ : '+plot_ebv],['$\chi^2 / DoF$ : '+chi]])
    label = transpose([['ObjID: '+strtrim(plot_id,2)],['$!8z!7 = $'+plot_z+'$\pm$'+plot_zerr],['$!8E!7(!8B-V!7)$ = '+plot_ebv+'$\pm$'+plot_ebvsig],['$\chi^2 / DoF$ = '+plot_chi]])
    for i = 0,nobj-1 do begin
        if (i eq 0) then current = 0 else current = 1
        ig = where(plot_bin[*,i],/null)
        p = plot(objwav[ig,i],plot_flux[ig,i],_extra=e,position=pos[*,i],/DEVICE,current=current)
        if (i eq 0 or i eq 1) then p.axes[0].showtext = 0
        if (i eq 1 or i eq 3) then p.axes[1].showtext = 0
        ;for t = 0,ntemps-1 do re = execute('p = plot(tempwav[*,i],'+temps[t]+'[*,i],col=col[*,t],thick=2,linestyle=linestyle[t],/ov)')   ;; plot models
        p_agn = plot(tempwav[*,i],agn[*,i],col=col[*,0],thick=2,linestyle=linestyle[0],/ov,name=' AGN')
        p_ell = plot(tempwav[*,i],ell[*,i],col=col[*,1],thick=2,linestyle=linestyle[1],/ov,name=' ELL')
        p_sfg = plot(tempwav[*,i],sfg[*,i],col=col[*,2],thick=2,linestyle=linestyle[2],/ov,name=' SFG')
        p_irr = plot(tempwav[*,i],irr[*,i],col=col[*,3],thick=2,linestyle=linestyle[3],/ov,name=' IRR')
        p = plot(tempwav[*,i],model[*,i],col='dark slate grey',thick=2,/ov)                                                           ;; plot coadded models
        p = errorplot(objwav[ig,i],plot_flux[ig,i],plot_err[ig,i],'o',SYM_FILLED=1,LINESTYLE='',sym_size=1.5,errorbar_thick=2,/OV)          ;; plot photometry
        for t = 0,ntemps-1 do begin
            lab = text(0.07,e.yra[1]-0.32*(1.8+t),label[t,i],font_size=14,/DATA,target=p,font_name='Times')
            ;fit = text(1.95,yp-t*0.35,temps[t]+': '+coeff[t,i],col=col[*,t],font_size=12,/DATA,TARGET=p,font_name='Times')
        endfor
        if (i eq 1) then l = legend(target=[p_agn,p_ell,p_sfg,p_irr],position=[0.9,0.92],/auto_text_color,sample_width=0.1,horizontal_spacing=0.06,font_name='Times')
    endfor
    xt = text(0.5,0.02,'$\lambda (observed) [ \mum ]$',alignment=0.5,font_size=14,font_name='Times')
    yt = text(0.025,0.5,'$log  \nu!8F!7_\nu  [erg s^{-1}cm^{-2}]$',orientation=90,alignment=0.5,font_size=14,font_name='Times')
    
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/sed_models.eps',/BITMAP else $
                                                     p.save,'figures/sed_models.png';,resolution=20
    endif
endif



;;----------------------------------------------------------------------------------------
;; LX vs LIR
;;----------------------------------------------------------------------------------------
if keyword_set(lx_lir) then begin
    print, '    PREPARING PAPER LX vs LIR'
    
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
    ;lx[where(lx eq 0. and iifinal_det_cha)] = lx_cha[where(lx eq 0. and iifinal_det_cha)]
    ;lx[where(lx eq 0. and iifinal_det_xmm)] = lx_xmm[where(lx eq 0. and iifinal_det_xmm)]
    ;lx[where(lx eq 0. and iifinal_det_nst)] = lx_nst[where(lx eq 0. and iifinal_det_nst)]

    iqso1 = where(iifinal_det and ebv le 0.15)
    iqso2 = where(iifinal_det and ebv gt 0.15)
    
    ;; test (1+z) correction
    loglirz = alog10((10.^loglir)/(1+z))
    loglirz[where(loglir eq -9999)] = -9999
    ;loglir = loglirz

    binsz = 0.1
    hist1 = hist_2d(loglir[iqso1],loglx[iqso1],bin1=binsz,min1=41.,max1=47.5,bin2=binsz,min2=41.,max2=46.5)
    hist2 = hist_2d(loglir[iqso2],loglx[iqso2],bin1=binsz,min1=41.,max1=47.5,bin2=binsz,min2=41.,max2=46.5)
    ind1 = array_indices(hist1,where(hist1 eq max(hist1)))
    ind2 = array_indices(hist2,where(hist2 eq max(hist2)))
    xlum = [41.:47.5:binsz]
    ylum = [41.:46.5:binsz]
    
    e = {xra:[41.,47.5],yra:[41.,46.5], $
         font_name:'Times'}
    if keyword_set(show) then buff = 0 else buff = 1
    ;; [0,158,115]  [204,121,167]
    im1 = image('chen17_bw1.png',transparency=50,dimension=[640,877],position=[80,440,587,810],buffer=buff,/device)
    pnodata = plot(xrel_chen,yrel_chen,_extra=e,/current,position=im1.position)
    cqso1 = contour(hist1,xlum,ylum,c_thick=4,rgb_table=colortable(49,/reverse),/fill,/ov,c_label_show=0,name='$!8E(B-V)!7 \leq  0.15$')
    prel_chen = plot(xrel_chen,yrel_chen,'-',col=[0,158,115],thick=4,_extra=e,/current,position=im1.position,name='Chen+17')
    prel_fiore = plot(xrel_fiore,yrel_fiore,'-',col=[213,94,0],thick=4,/ov,name='Fiore+09')
    prel_chen.axes[0].showtext = 0
    chent = text(41.15,42.3,/data,target=prel_chen,'Chen+17',col=[0,158,115],font_size=14,font_style='Bold',font_name='Times',fill_background=1,fill_color='white')
    fioret = text(42.15,41.5,/data,target=prel_fiore,'Fiore+09',col=[213,94,0],font_size=14,font_style='Bold',font_name='Times',fill_background=1,fill_color='white')
    ;p = plot([median(xlum[ind1[0,*]])],[median(ylum[ind1[1,*]])],'o',col='black',sym_size=1.5,/ov)
    ;p = plot([median(xlum[ind1[0,*]])],[median(ylum[ind1[1,*]])],'S',col='black',sym_size=1.0,sym_filled=1,/ov)
    ;p = plot(xlum[ind1[0,*]],ylum[ind1[1,*]],'X',sym_size=1.5,sym_thick=1,/sym_filled,col=[0,158,115],/ov)
    im2 = image('chen17_bw2.png',transparency=50,/current,position=[80,75,587,445],/device)
    pnodata = plot(xrel_chen,yrel_chen,_extra=e,/current,position=im2.position)    
    cqso2 = contour(hist2,xlum,ylum,c_thick=4,rgb_table=colortable(62,/reverse),/fill,/ov,c_label_show=0,name='$!8E(B-V)!7 > 0.15$')
    prel_chen = plot(xrel_chen,yrel_chen,'-',col=[0,158,115],thick=4,_extra=e,/current,position=im2.position,name='Chen+17')    
    prel_fiore = plot(xrel_fiore,yrel_fiore,'-',col=[213,94,0],thick=4,/ov,name='Fiore+09')
    ;p = plot([median(xlum[ind2[0,*]])],[median(ylum[ind2[1,*]])],'o',col='black',sym_size=1.5,/ov)
    ;p = plot([median(xlum[ind2[0,*]]),0,0,0],[median(ylum[ind2[1,*]]),0,0,0],'S',col='black',sym_size=1.0,sym_filled=1,_extra=e,/ov)
    ;p = plot(xlum[ind2[0,*]],ylum[ind2[1,*]],'X',sym_size=1.5,sym_thick=3,/sym_filled,col=[0,158,115],/ov)
    leg = legend(target=[cqso1,cqso2],position=[0.15,0.41],/normal,vertical_alignment=0.,horizontal_alignment=0.,font_name='Times')   
    ;leg.position = [0.14,0.39]

    xt = text(0.5,0.04,'$log  !8L!7_{6 \mu m} [erg s^{-1}cm^{-2}]$',alignment=0.5,font_size=14,font_name='Times')
    yt = text(0.05,0.5,'$log  !8L!7_{2-10 keV} [erg s^{-1}cm^{-2}]$',orientation=90,alignment=0.5,font_size=14,font_name='Times')
        
    ;peak_circ = text(0.72,0.13,'$\U(25EF)$',col='dark slate grey',font_size=12,font_name='Times')
    ;peak_star = text(0.721,0.128,'$\U(2605)$',col='dark slate grey',font_size=14,font_name='Times')
    ;peak_text = text(0.76,0.128,'Peak Freq.',col='dark slate grey',font_size=12,font_name='Times')
    
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then im1.save,'figures/lx_v_lir.eps',/BITMAP else $
                                                     im1.save,'figures/lx_v_lir.png';,resolution=20
    endif
endif



;;----------------------------------------------------------------------------------------
;; Flux Limit
;;----------------------------------------------------------------------------------------
if keyword_set(flux_limit) then begin
    multi_sn = 1
    print, '    PREPARING PAPER FLUX LIMIT'
    dim = [1285,640]
    src_pos = [[80,75,465,585],[465,75,850,585],[850,75,1235,585]]
    axis_src = 0
    which_axis = 1
    leg_src = nfield-1
    leg_pos = [0.675,0.80]
    xpos = [0.5,0.035]
    ypos = [0.025,0.5]
    
    e = {yra:[-15.5,-9.], $
         xlog:1,ylog:0, $
         font_name:'Times', $
         buffer:1}
    if keyword_set(show) then e.buffer = 0
;    xra = [[8e2,max(texp_cha)],[8e2,max(texp_xmm)],[2e3,max(texp_nst)>3e5]]
    xra = [[8e2,3e5],[8e2,3e5],[2e3,3e5]]
    cat = ['   Chandra','XMM-Newton','    NuSTAR']
    
    for i = 0,nfield-1 do begin
        if (i eq 0) then current = 0 else current = 1
        re = execute('p = plot(CAT_EXP'+xfield[i]+',alog10(CAT_FLX'+xfield[i]+'),_extra=e,xra=xra[*,i],dimension=dim,position=src_pos[*,i],/device,current=current,/nodata)')
        pcat_dark = plot([1e-5,1e-5],[1e5,1e6],'s',col='dark grey',sym_filled=1,sym_size=0.5,transparency=75,/ov,name=' X-ray Catalog(s)')
        re = execute('ncat = n_elements(CAT_EXP'+xfield[i]+')')
        ;if (ncat gt 1e5) then re = execute('irand = round(randomu(seed,n_elements(CAT_EXP'+xfield[i]+')/10)*n_elements(CAT_EXP'+xfield[i]+'))') else $
                              re = execute('irand = lindgen(n_elements(CAT_EXP'+xfield[i]+'))')
        re = execute('pcat_lite = plot(CAT_EXP'+xfield[i]+'[irand],alog10(CAT_FLX'+xfield[i]+'[irand]),"s",col="light grey",sym_filled=1,sym_size=0.5,transparency=85,/ov)')
        ; plot open symbols first (makes for a cleaner plot)
        ;re = execute('inrm = where(IIAGN_NRM'+xfield[i]+',ctnrm)')
        ;if (ctnrm gt 0) then re = execute('pnrm = plot(TEXP'+xfield[i]+'[inrm],alog10(fxir[inrm]),"o",sym_filled=0,sym_size=0.5,color="orange",transparency=75,/ov)')
        re = execute('pnon = plot(TEXP'+xfield[i]+'[where(iifinal_non'+xfield[i]+')],logfxir[where(iifinal_non'+xfield[i]+')],"o",sym_filled=1,sym_size=0.5,color="orange",transparency=50,/ov,name=" X-ray non-det. [$!8F!7_{X,expected}$]")')
        ;; off plot range
        re = execute('ioff = where(iifinal_non'+xfield[i]+' and TEXP'+xfield[i]+' gt xra[1,i],noff)')
        if (noff gt 0) then pwoff = arrow(transpose([[make_array(noff,value=xra[1,i]*0.75)],[make_array(noff,value=xra[1,i])]]),transpose([[logfxir[ioff]],[logfxir[ioff]]]),color='orange',/ov,/data,thick=2,head_size=0.5,target=pnon)

        ;re = execute('idrm = where(IIAGN_DRM'+xfield[i]+',ctdrm)')
        ;if (ctdrm gt 0) then begin
        ;    re = execute('pdrm = plot(EXP'+xfield[i]+'[idrm],alog10(FLX'+xfield[i]+'[idrm]),"S",sym_filled=0,sym_size=1.,color="dodger blue",/ov)')
        ;    ;re = execute('pdro = plot(EXP'+xfield[i]+'[idrm],alog10(FLX'+xfield[i]+'[idrm]),"S",sym_filled=0,sym_size=1.5,transparency=85,/ov)')
        ;endif
        re = execute('pdet = plot(EXP'+xfield[i]+'[where(iifinal_det'+xfield[i]+')],alog10(FX'+xfield[i]+'[where(iifinal_det'+xfield[i]+')]),"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov,name="X-ray detected [$!8F!7_{X,observed}$]")')
        re = execute('pdo = plot(EXP'+xfield[i]+'[where(iifinal_det'+xfield[i]+')],alog10(FX'+xfield[i]+'[where(iifinal_det'+xfield[i]+')]),"S",sym_filled=0,sym_size=1.5,transparency=85,/ov)')
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
                snt = text(2.5e3,-12.5,'$!8S/N!7 \geq  5.0$',target=p,/DATA,font_size=12,font_name='Times')
                snt = text(2.5e3,-12.9,'$!8S/N!7 \geq  3.0$',target=p,/DATA,font_size=12,font_name='Times',col='red')
            endif
        endif        
        
        ;; ensure axes range
        p.xra = xra[*,i]
        
        ;; average error bar
        xpt = [1.5e3,1.5e3,3.8e3]
        ypt = [-15.]
        perr = plot([xpt[i]],ypt,'s',sym_size=4,sym_thick=2,/ov)
        re = execute('yer = [median(CAT_ERR'+xfield[i]+'/(alog(10.)*CAT_FLX'+xfield[i]+'))]')
        perr = errorplot([xpt[i]],ypt,yer,errorbar_capsize=0.1,linestyle='',/ov)
        perr = plot([xpt[i]],ypt,'S',col='white',sym_size=1.5,sym_filled=1,/ov)
        perr = plot([xpt[i]],ypt,'S',col='dodger blue',sym_size=1.5,sym_filled=0,/ov)
        if (i eq 2) then terr = text(xpt[i]*1.4,ypt,target=perr,/data,'$Median \sigma_{log F}$',vertical_alignment=0.5,font_size=12,font_name='Times')
        
        if (i ne axis_src) then p.axes[which_axis].showtext = 0
        t = text(0.94,0.9,cat[i],target=p,/RELATIVE,font_size=16,font_style='Bold italic',alignment=1.,font_name='Times')
        if (i eq leg_src) then l = legend(target=[pdet,pnon,pcat_dark],position=leg_pos,/normal,/auto_text_color,sample_width=0,horizontal_spacing=0.06,vertical_alignment=0.5,horizontal_alignment=0)
    endfor
    xt = text(xpos[0],xpos[1],'$exposure time [s]$',alignment=0.5,font_size=14,font_name='Times')
    yt = text(ypos[0],ypos[1],'$log  !8F!7_{2-10 keV} [erg s^{-1}cm^{-2}]$',orientation=90,alignment=0.5,font_size=14,font_name='Times')    
    
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/fx_limits.eps',/BITMAP else $
                                                     p.save,'figures/fx_limits.png';,resolution=20
    endif
endif



;;----------------------------------------------------------------------------------------
;; Luminosity Ratio
;;----------------------------------------------------------------------------------------
if keyword_set(lum_ratio) then begin
    print, '    PREPARING LUMINOSITY RATIOS'
    
    model = 'POWER'
    ;model = 'BORUS'
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; PLOT L/L VS E(B-V)
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ell = {xra:[-3.8,2.2],yra:[-3.,1.], $
           font_name:'Times', $
           dimension:[1130,880], $
           buffer:1}
    if keyword_set(show) then ell.buffer = 0

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
    pw = plot(lebv[where(iifinal_non)],llnon[where(iifinal_non)],_extra=ell,position=[75,445,585,825],/device,/nodata,ytickvalues=[-3.0,-2.0,-1.0,0.0,1.0],ytickformat='(d4.1)')
    ;; unobscured shading
    pw = plot(ell.xra,[1,1]*ell.yra[1],linestyle='',fill_background=1,fill_level=0.,fill_color='light grey',fill_transparency=80,/ov)
    ;; CT shading
    pw = plot(ell.xra,ll_bound[1]*[1,1],linestyle='',fill_background=1,fill_level=ell.yra[0],fill_color='light grey',fill_transparency=80,/ov)
    ;; plot data
    pwllnon = arrow(transpose([[lebv[where(iifinal_non and iidet_wac)]],[lebv[where(iifinal_non and iidet_wac)]]]),transpose([[llnon[where(iifinal_non and iidet_wac)]],[llnon[where(iifinal_non and iidet_wac)]-0.12]]),color="orange",fill_transparency=75,target=[pw],/data,head_size=0.3,clip=1,name=" X-ray non-det.")
    ;pwllnon = plot(lebv[where(iifinal_non and iidet_wac)],llnon[where(iifinal_non and iidet_wac)],"o",sym_filled=1,sym_size=0.5,color="orange",transparency=50,/ov,name=" X-ray non-det.")    
    pwlldet = plot(lebv[where(iifinal_det and iidet_wac)],lldet[where(iifinal_det and iidet_wac)],"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov,name=" X-ray detected")
    p = plot(lebv[where(iifinal_det and iidet_wac)],lldet[where(iifinal_det and iidet_wac)],"S",sym_filled=0,sym_size=1.5,transparency=85,/ov)
    ;; non-detections off plot range
    ioff = where(iifinal_non and iidet_wac and llnon ge ell.yra[1],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[1]-0.12)],[make_array(noff,value=ell.yra[1])-0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ioff = where(iifinal_non and iidet_wac and llnon le ell.yra[0],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[0]+0.12)],[make_array(noff,value=ell.yra[0])+0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    p = plot(lebv[where(iifinal_det and iidet_wac)],lldet[where(iifinal_det and iidet_wac)],"S",sym_filled=0,sym_size=1.5,transparency=85,/ov)
    ;; detections off plot range
    ioff = where(iifinal_det and iidet_wac and lldet gt ell.yra[1],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[1]-0.12)],[make_array(noff,value=ell.yra[1])-0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ioff = where(iifinal_det and iidet_wac and lldet lt ell.yra[0],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[0]+0.12)],[make_array(noff,value=ell.yra[0])+0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; median error bars
    perr = plot([-2.2],[ell.yra[0]+0.5],'s',sym_size=4,sym_thick=2,/ov)
    perr = errorplot([-2.2],[ell.yra[0]+0.5],[median(e_lldet[where(iifinal_det and iidet_wac)])],errorbar_capsize=0.1,linestyle="",/ov)
    perr = plot([-2.2],[ell.yra[0]+0.5],"S",sym_filled=1,sym_size=1.5,col='white',/ov)
    perr = plot([-2.2],[ell.yra[0]+0.5] ,"S",sym_filled=0,sym_size=1.5,col='dodger blue',/ov)
    ;perr = arrow([-1.9,-1.9],[-3.,-3.3],/data,target=pwllnon,color='orange',/ov,thick=2,head_size=0.5)
    ;perr = plot([-1.9],ell.yra[0]+0.4,"o",sym_filled=1,sym_size=1.,color="white",transparency=0,/ov)
    ;perr = plot([-1.9],ell.yra[0]+0.4,"o",sym_filled=0,sym_size=1.,color="orange",transparency=0,/ov)
    ;; CT lines
    p = plot(ell.xra[1]+[-0.7,0.],ll_bound[0]*[1,1],'-',thick=2,/ov)
    p = plot(ell.xra[1]+[-0.7,0.],ll_bound[1]*[1,1],'-',thick=2,/ov)
    ;; catalog sources text
    t = text(ell.xra[0]+0.22,0.42,'!16WISE!15 AGN',col='red',/data,target=[pwllnon],font_size=12,font_name='Times',font_style='Bold')
    t = text(ell.xra[0]+0.22,0.17,'sources',col='red',/data,target=[pwllnon],font_size=12,font_name='Times',font_style='Bold')
    t = text(ell.xra[1]-0.2,ll_bound[0]+0.15,'Unobscured',col='black',target=pwllnon,alignment=1.,font_size=12,font_name='Times',font_style='Bold',/data)
    t = text(ell.xra[1]-0.2,ll_bound[1]-0.3,'Compton',col='black',target=pwllnon,alignment=1.,font_size=12,font_name='Times',font_style='Bold',/data)
    t = text(ell.xra[1]-0.2,ll_bound[1]-0.55,'thick',col='black',target=pwllnon,alignment=1.,font_size=12,font_name='Times',font_style='Bold',/data)
    ;; NH lines + text
    t = text(ell.xra[0]+0.2,ell.yra[0]+0.4,'$!8N!7_H  [cm^{-2}]$',col='black',/data,target=[pwllnon],font_size=11,font_name='Times')
    for i = 0,n_elements(ll_lines)-1 do begin
        p = plot(nh_lines,ll_lines[i]*[1.,1.],'--',col='black',thick=2,/ov)
        t = text(ell.xra[0]+0.2,ypos[i],nhtext[i],col='black',/data,target=[pwllnon],font_size=11,font_name='Times')
    endfor
    ;; remove axes
    pwlldet.axes[0].showtext=0
    ;; model choice
    mt = text(-0.8,ell.yra[0]+0.5,'Power law model',target=pwlldet,alignment=0.5,vertical_alignment=0.5,font_size=14,font_name='Times',font_style='Bold',/data)

    ;; Remaining sources
    pr = plot(lebv[where(iifinal_non)],llnon[where(iifinal_non)],_extra=ell,position=[585,445,1095,825],current=1,/device,/nodata)
    ;; unobscured shading
    pr = plot(ell.xra,[1,1]*ell.yra[1],linestyle='',fill_background=1,fill_level=0.,fill_color='light grey',fill_transparency=80,/ov)
    ;; CT shading
    pr = plot(ell.xra,ll_bound[1]*[1,1],linestyle='',fill_background=1,fill_level=ell.yra[0],fill_color='light grey',fill_transparency=80,/ov)
    ;; plot data
    prllnon = arrow(transpose([[lebv[where(iifinal_non and ~iidet_wac)]],[lebv[where(iifinal_non and ~iidet_wac)]]]),transpose([[llnon[where(iifinal_non and ~iidet_wac)]],[llnon[where(iifinal_non and ~iidet_wac)]-0.12]]),color="orange",fill_transparency=75,target=[pr],/data,head_size=0.3,clip=1,name=" X-ray non-det.")
    ;prllnon = plot(lebv[where(iifinal_non and ~iidet_wac)],llnon[where(iifinal_non and ~iidet_wac)],"o",sym_filled=1,sym_size=0.5,color="orange",transparency=50,/ov,name=" X-ray non-det.")
    prlldet = plot(lebv[where(iifinal_det and ~iidet_wac)],lldet[where(iifinal_det and ~iidet_wac)],"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov,name=" X-ray detected")
    p = plot(lebv[where(iifinal_det and ~iidet_wac)],lldet[where(iifinal_det and ~iidet_wac)],"S",sym_filled=0,sym_size=1.5,transparency=85,/ov)
    ;; non-detections off plot range
    ioff = where(iifinal_non and ~iidet_wac and llnon ge ell.yra[1],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[1]-0.12)],[make_array(noff,value=ell.yra[1])-0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ioff = where(iifinal_non and ~iidet_wac and llnon le ell.yra[0],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[0]+0.12)],[make_array(noff,value=ell.yra[0])+0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; detections off plot range
    ioff = where(iifinal_det and ~iidet_wac and lldet gt ell.yra[1],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[1]-0.12)],[make_array(noff,value=ell.yra[1])-0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ioff = where(iifinal_det and ~iidet_wac and lldet lt ell.yra[0],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[0]+0.12)],[make_array(noff,value=ell.yra[0])+0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; median error bars
    perr = plot([-2.2],[ell.yra[0]+0.5],'s',sym_size=4,sym_thick=2,/ov)
    perr = errorplot([-2.2],[ell.yra[0]+0.5],[median(e_lldet[where(iifinal_det and ~iidet_wac)])],errorbar_capsize=0.1,linestyle="",/ov)
    perr = plot([-2.2],[ell.yra[0]+0.5],"S",sym_filled=1,sym_size=1.5,col='white',/ov)
    perr = plot([-2.2],[ell.yra[0]+0.5],"S",sym_filled=0,sym_size=1.5,col='dodger blue',/ov)
    terr = text([-2.2]+0.3,[ell.yra[0]+0.5],target=perr,/data,'$Median \sigma_{!8R_{L}!7}$',vertical_alignment=0.5,font_size=12,font_name='Times')
    ;perr = arrow([-1.9,-1.9],[-3.,-3.3],/data,target=prllnon,color='orange',/ov,thick=2,head_size=0.5)
    ;perr = plot([-1.9],ell.yra[0]+0.4,"o",sym_filled=1,sym_size=1.,color="white",transparency=0,/ov)
    ;perr = plot([-1.9],ell.yra[0]+0.4,"o",sym_filled=0,sym_size=1.,color="orange",transparency=0,/ov)
    ;; CT lines
    p = plot(ell.xra[1]+[-0.7,0.],ll_bound[0]*[1,1],'-',thick=2,/ov)
    p = plot(ell.xra[1]+[-0.7,0.],ll_bound[1]*[1,1],'-',thick=2,/ov)
    ;; Catalog sources text
    t = text(ell.xra[0]+0.22,0.42,'Remaining',col='red',/data,target=[prllnon],font_size=12,font_name='Times',font_style='Bold')
    t = text(ell.xra[0]+0.22,0.17,'sources',col='red',/data,target=[prllnon],font_size=12,font_name='Times',font_style='Bold')
    ;; NH lines + text
    t = text(ell.xra[0]+0.2,ell.yra[0]+0.4,'$!8N!7_H  [cm^{-2}]$',col='black',/data,target=[prllnon],font_size=11,font_name='Times')
    for i = 0,n_elements(ll_lines)-1 do begin
        p = plot(nh_lines,ll_lines[i]*[1.,1.],'--',col='black',thick=2,/ov)
        t = text(ell.xra[0]+0.2,ypos[i],nhtext[i],col='black',/data,target=[prllnon],alignment=0.,font_size=11,font_name='Times')
    endfor
    
    ;; remove axes
    prlldet.axes[1].showtext = 0
    prlldet.axes[0].showtext = 0
    
    ;; legend
    ;l = legend(target=[prdet,prnon],position=[0.865,0.81],/normal,/auto_text_color,sample_width=0.,horizontal_spacing=0.06,font_name='Times')
    ;l.position = [0.956,0.874]
    
    ;xt = text(0.5,0.04,'$log  !8E!7(!8B-V!7)$',alignment=0.5,font_size=14,font_name='Times')
    yt = text(0.02,0.75,'$log  !8L!7_X / !8L!7_X(!8L!7_{IR})  (2$-$10 keV)$',orientation=90,alignment=0.5,font_size=14,font_name='Times')        
    yt.position = [0.007524362,0.59069969,0.026571130,0.85248213]

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; PLOT NH VS E(B-V)
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    enh = {xra:[-3.8,2.2],yra:[20.5,25.5], $
           font_name:'Times', $
           dimension:[1130,880], $
           buffer:1}
    if keyword_set(show) then enh.buffer = 0
    
    nhxdet = rl2nh(lldet,model=model)
    nhxnon = rl2nh(llnon,model=model)
    
    ;; [unobscured,Compton-thick]     
    nh_bound = alog10([1e21,1.5e24])

    ;; WISE AGN Catalog sources
    pw = plot(lebv[where(iifinal_non and iidet_wac)],llnon[where(iifinal_non and iidet_wac)],_extra=enh,position=[75,65,585,445],current=1,/device,/nodata)
    ;; CT shading
    pw = plot(enh.xra,nh_bound[1]*[1,1],linestyle='',fill_background=1,fill_level=enh.yra[1],fill_color='light grey',fill_transparency=80,/ov)
    ;; plot data
    pwnhnon = arrow(transpose([[lebv[where(iifinal_non and iidet_wac)]],[lebv[where(iifinal_non and iidet_wac)]]]),transpose([[nhxnon[where(iifinal_non and iidet_wac)]],[nhxnon[where(iifinal_non and iidet_wac)]+0.12]]),color="orange",fill_transparency=75,target=[pw],/data,head_size=0.3,clip=1,name=" X-ray non-det.")
    ;pwnhnon = plot(lebv[where(iifinal_non and iidet_wac)],nhxnon[where(iifinal_non and iidet_wac)],"o",sym_filled=1,sym_size=0.5,color="orange",transparency=50,/ov,name=" X-ray non-det.")
    pwnhdet = plot(lebv[where(iifinal_det and iidet_wac)],nhxdet[where(iifinal_det and iidet_wac)],"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov,name=" X-ray detected")
    p = plot(lebv[where(iifinal_det and iidet_wac)],nhxdet[where(iifinal_det and iidet_wac)],"S",sym_filled=0,sym_size=1.5,transparency=85,/ov)
    ;; non-detections off plot range
    ioff = where(iifinal_non and iidet_wac and nhxnon ge enh.yra[1],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[1]-0.12)],[make_array(noff,value=enh.yra[1])-0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ioff = where(iifinal_non and iidet_wac and nhxnon le enh.yra[0],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[0]+0.12)],[make_array(noff,value=enh.yra[0])+0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; detections off plot range-0.2
    ioff = where(iifinal_det and iidet_wac and nhxdet gt enh.yra[1],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[1]-0.12)],[make_array(noff,value=enh.yra[1])-0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ioff = where(iifinal_det and iidet_wac and nhxdet lt enh.yra[0],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[0]+0.12)],[make_array(noff,value=enh.yra[0])+0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; CT lines
    p = plot(enh.xra[0]+[0.,0.7],nh_bound[1]*[1,1],'-',thick=2,/ov)
    ;; catalog sources text
    t = text(enh.xra[0]+0.22,nh_bound[1]+0.24,'Compton',col='black',target=pwnhnon,font_size=12,font_name='Times',font_style='Bold',/data)
    t = text(enh.xra[0]+0.22,nh_bound[1]+0.09,'thick',col='black',target=pwnhnon,font_size=12,font_name='Times',font_style='Bold',/data)

    ;; Remaining sources
    pr = plot(lebv[where(iifinal_non and ~iidet_wac)],nhxnon[where(iifinal_non and ~iidet_wac)],_extra=enh,position=[585,65,1095,445],current=1,/device,/nodata)
    ;; CT shading
    pr = plot(enh.xra,nh_bound[1]*[1,1],linestyle='',fill_background=1,fill_level=enh.yra[1],fill_color='light grey',fill_transparency=80,/ov)
    ;; plot data
    prnhnon = arrow(transpose([[lebv[where(iifinal_non and ~iidet_wac)]],[lebv[where(iifinal_non and ~iidet_wac)]]]),transpose([[nhxnon[where(iifinal_non and ~iidet_wac)]],[nhxnon[where(iifinal_non and ~iidet_wac)]+0.12]]),color="orange",fill_transparency=75,target=[pr],/data,head_size=0.3,clip=1,name=" X-ray non-det.")
    ;prnhnon = plot(lebv[where(iifinal_non and ~iidet_wac)],nhxnon[where(iifinal_non and ~iidet_wac)],"o",sym_filled=1,sym_size=0.5,color="orange",transparency=50,/ov,name=" X-ray non-det.")
    prnhdet = plot(lebv[where(iifinal_det and ~iidet_wac)],nhxdet[where(iifinal_det and ~iidet_wac)],"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov,name=" X-ray detected")
    p = plot(lebv[where(iifinal_det and ~iidet_wac)],nhxdet[where(iifinal_det and ~iidet_wac)],"S",sym_filled=1,sym_size=1.5,transparency=85,/ov)
    ;; non-detections off plot range
    ioff = where(iifinal_non and ~iidet_wac and nhxnon ge enh.yra[1],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[1]-0.12)],[make_array(noff,value=enh.yra[1])-0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ioff = where(iifinal_non and ~iidet_wac and nhxnon le enh.yra[0],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[0]+0.12)],[make_array(noff,value=enh.yra[0])+0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; detections off plot range-0.2
    ioff = where(iifinal_det and ~iidet_wac and nhxdet gt enh.yra[1],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[1]-0.12)],[make_array(noff,value=enh.yra[1])-0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ioff = where(iifinal_det and ~iidet_wac and nhxdet lt enh.yra[0],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[0]+0.12)],[make_array(noff,value=enh.yra[0])+0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; CT lines
    p = plot(enh.xra[0]+[0.,0.7],nh_bound[1]*[1,1],'-',thick=2,/ov)
    ;; catalog sources text
    ;t = text(e.xra[0]+0.2,nh_bound[1]+0.2,'Compton',col='black',target=pwnon,font_size=12,font_name='Times',font_style='Bold',/data)
    ;t = text(e.xra[0]+0.2,nh_bound[1]+0.05,'thick',col='black',target=pwnon,font_size=12,font_name='Times',font_style='Bold',/data)
    ;; remove axes
    prnhdet.axes[1].showtext=0    
    ;; legend
    l = legend(target=[prnhdet,prnhnon],/normal,/auto_text_color,sample_width=0.,horizontal_spacing=0.06,font_name='Times')
    l.position = [0.69,0.49]
    
    xt = text(0.52,0.02,'$log  !8E!7(!8B-V!7)$',alignment=0.5,font_size=14,font_name='Times')
    yt = text(0.02,0.25,'$log  !8N!7_H [cm^{-2}]$',orientation=90,alignment=0.5,font_size=14,font_name='Times')        
    yt.position = [0.007524362,0.22634560,0.026571130,0.35319986]
    
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/rlum_b.eps',/BITMAP else $
                                                     p.save,'figures/rlum_b.png';,resolution=20                             
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
    pw = plot(lebv[where(iifinal_non)],llnon[where(iifinal_non)],_extra=ell,position=[75,445,585,825],/device,/nodata,ytickvalues=[-3.0,-2.0,-1.0,0.0,1.0],ytickformat='(d4.1)')
    ;; unobscured shading
    pw = plot(ell.xra,[1,1],linestyle='',fill_background=1,fill_level=0.,fill_color='light grey',fill_transparency=80,/ov)
    ;; CT shading
    pw = plot(ell.xra,ll_bound[1]*[1,1],linestyle='',fill_background=1,fill_level=ell.yra[0],fill_color='light grey',fill_transparency=80,/ov)
    ;; plot data
    pwllnon = arrow(transpose([[lebv[where(iifinal_non and iidet_wac)]],[lebv[where(iifinal_non and iidet_wac)]]]),transpose([[llnon[where(iifinal_non and iidet_wac)]],[llnon[where(iifinal_non and iidet_wac)]-0.12]]),color="orange",fill_transparency=75,target=[pw],/data,head_size=0.3,clip=1,name=" X-ray non-det.")
    ;pwllnon = plot(lebv[where(iifinal_non and iidet_wac)],llnon[where(iifinal_non and iidet_wac)],"o",sym_filled=1,sym_size=0.5,color="orange",transparency=50,/ov,name=" X-ray non-det.")
    pwlldet = plot(lebv[where(iifinal_det and iidet_wac)],lldet[where(iifinal_det and iidet_wac)],"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov,name=" X-ray detected")
    p = plot(lebv[where(iifinal_det and iidet_wac)],lldet[where(iifinal_det and iidet_wac)],"S",sym_filled=0,sym_size=1.5,transparency=85,/ov)
    ;; non-detections off plot range
    ioff = where(iifinal_non and iidet_wac and llnon ge ell.yra[1],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[1]-0.12)],[make_array(noff,value=ell.yra[1])-0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ioff = where(iifinal_non and iidet_wac and llnon le ell.yra[0],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[0]+0.12)],[make_array(noff,value=ell.yra[0])+0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    p = plot(lebv[where(iifinal_det and iidet_wac)],lldet[where(iifinal_det and iidet_wac)],"S",sym_filled=0,sym_size=1.5,transparency=85,/ov)
    ;; detections off plot range
    ioff = where(iifinal_det and iidet_wac and lldet gt ell.yra[1],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[1]-0.12)],[make_array(noff,value=ell.yra[1])-0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ioff = where(iifinal_det and iidet_wac and lldet lt ell.yra[0],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[0]+0.12)],[make_array(noff,value=ell.yra[0])+0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; median error bars
    perr = plot([-2.2],[ell.yra[0]+0.5],'s',sym_size=4,sym_thick=2,/ov)
    perr = errorplot([-2.2],[ell.yra[0]+0.5],[median(e_lldet[where(iifinal_det and iidet_wac)])],errorbar_capsize=0.1,linestyle="",/ov)
    perr = plot([-2.2],[ell.yra[0]+0.5],"S",sym_filled=1,sym_size=1.5,col='white',/ov)
    perr = plot([-2.2],[ell.yra[0]+0.5],"S",sym_filled=0,sym_size=1.5,col='dodger blue',/ov)
    ;perr = arrow([-1.9,-1.9],[-3.,-3.3],/data,target=pwllnon,color='orange',/ov,thick=2,head_size=0.5)
    ;perr = plot([-1.9],ell.yra[0]+0.4,"o",sym_filled=1,sym_size=1.,color="white",transparency=0,/ov)
    ;perr = plot([-1.9],ell.yra[0]+0.4,"o",sym_filled=0,sym_size=1.,color="orange",transparency=0,/ov)
    ;; CT lines
    p = plot(ell.xra[1]+[-0.7,0.],ll_bound[0]*[1,1],'-',thick=2,/ov)
    p = plot(ell.xra[1]+[-0.7,0.],ll_bound[1]*[1,1],'-',thick=2,/ov)
    ;; catalog sources text
    t = text(ell.xra[0]+0.22,0.42,'!16WISE!15 AGN',col='red',/data,target=[pwllnon],font_size=12,font_name='Times',font_style='Bold')
    t = text(ell.xra[0]+0.22,0.17,'sources',col='red',/data,target=[pwllnon],font_size=12,font_name='Times',font_style='Bold')
    t = text(ell.xra[1]-0.2,ll_bound[0]+0.15,'Unobscured',col='black',target=pwllnon,alignment=1.,font_size=12,font_name='Times',font_style='Bold',/data)
    t = text(ell.xra[1]-0.2,ll_bound[1]-0.3,'Compton',col='black',target=pwllnon,alignment=1.,font_size=12,font_name='Times',font_style='Bold',/data)
    t = text(ell.xra[1]-0.2,ll_bound[1]-0.55,'thick',col='black',target=pwllnon,alignment=1.,font_size=12,font_name='Times',font_style='Bold',/data)
    ;; NH lines + text
    t = text(ell.xra[0]+0.2,ell.yra[0]+0.4,'$!8N!7_H  [cm^{-2}]$',col='black',/data,target=[pwllnon],font_size=11,font_name='Times')
    for i = 0,n_elements(ll_lines)-1 do begin
        p = plot(nh_lines,ll_lines[i]*[1.,1.],'--',col='black',thick=2,/ov)
        t = text(ell.xra[0]+0.2,ypos[i],nhtext[i],col='black',/data,target=[pwllnon],font_size=11,font_name='Times')
    endfor
    ;; remove axes
    pwlldet.axes[0].showtext=0
    ;; model choice
    mt = text(-0.8,ell.yra[0]+0.5,'BORUS model',target=pwlldet,alignment=0.5,vertical_alignment=0.5,font_size=14,font_name='Times',font_style='Bold',/data)
    
    ;; Remaining sources
    pr = plot(lebv[where(iifinal_non)],llnon[where(iifinal_non)],_extra=ell,position=[585,445,1095,825],current=1,/device,/nodata)
    ;; unobscured shading
    pr = plot(ell.xra,[1,1],linestyle='',fill_background=1,fill_level=0.,fill_color='light grey',fill_transparency=80,/ov)
    ;; CT shading
    pr = plot(ell.xra,ll_bound[1]*[1,1],linestyle='',fill_background=1,fill_level=ell.yra[0],fill_color='light grey',fill_transparency=80,/ov)
    ;; data
    prllnon = arrow(transpose([[lebv[where(iifinal_non and ~iidet_wac)]],[lebv[where(iifinal_non and ~iidet_wac)]]]),transpose([[llnon[where(iifinal_non and ~iidet_wac)]],[llnon[where(iifinal_non and ~iidet_wac)]-0.12]]),color="orange",fill_transparency=75,target=[pr],/data,head_size=0.3,clip=1,name=" X-ray non-det.")
    ;prllnon = plot(lebv[where(iifinal_non and ~iidet_wac)],llnon[where(iifinal_non and ~iidet_wac)],"o",sym_filled=1,sym_size=0.5,color="orange",transparency=50,/ov,name=" X-ray non-det.")
    prlldet = plot(lebv[where(iifinal_det and ~iidet_wac)],lldet[where(iifinal_det and ~iidet_wac)],"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov,name=" X-ray detected")
    p = plot(lebv[where(iifinal_det and ~iidet_wac)],lldet[where(iifinal_det and ~iidet_wac)],"S",sym_filled=0,sym_size=1.5,transparency=85,/ov)
    ;; non-detections off plot range
    ioff = where(iifinal_non and ~iidet_wac and llnon ge ell.yra[1],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[1]-0.12)],[make_array(noff,value=ell.yra[1])-0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ioff = where(iifinal_non and ~iidet_wac and llnon le ell.yra[0],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[0]+0.12)],[make_array(noff,value=ell.yra[0])+0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; detections off plot range
    ioff = where(iifinal_det and ~iidet_wac and lldet gt ell.yra[1],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[1]-0.12)],[make_array(noff,value=ell.yra[1])-0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ioff = where(iifinal_det and ~iidet_wac and lldet lt ell.yra[0],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=ell.yra[0]+0.12)],[make_array(noff,value=ell.yra[0])+0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; median error bars
    perr = plot([-2.2],[ell.yra[0]+0.5],'s',sym_size=4,sym_thick=2,/ov)
    perr = errorplot([-2.2],[ell.yra[0]+0.5],[median(e_lldet[where(iifinal_det and ~iidet_wac)])],errorbar_capsize=0.1,linestyle="",/ov)
    perr = plot([-2.2],[ell.yra[0]+0.5],"S",sym_filled=1,sym_size=1.5,col='white',/ov)
    perr = plot([-2.2],[ell.yra[0]+0.5],"S",sym_filled=0,sym_size=1.5,col='dodger blue',/ov)
    terr = text([-2.2]+0.3,[ell.yra[0]+0.5],target=perr,/data,'$Median \sigma_{!8R_{L}!7}$',vertical_alignment=0.5,font_size=12,font_name='Times')
    ;perr = arrow([-1.9,-1.9],[-3.,-3.3],/data,target=prllnon,color='orange',/ov,thick=2,head_size=0.5)
    ;perr = plot([-1.9],ell.yra[0]+0.4,"o",sym_filled=1,sym_size=1.,color="white",transparency=0,/ov)
    ;perr = plot([-1.9],ell.yra[0]+0.4,"o",sym_filled=0,sym_size=1.,color="orange",transparency=0,/ov)
    ;; CT lines
    p = plot(ell.xra[1]+[-0.7,0.],ll_bound[0]*[1,1],'-',thick=2,/ov)
    p = plot(ell.xra[1]+[-0.7,0.],ll_bound[1]*[1,1],'-',thick=2,/ov)
    ;; Catalog sources text
    t = text(ell.xra[0]+0.22,0.42,'Remaining',col='red',/data,target=[prllnon],font_size=12,font_name='Times',font_style='Bold')
    t = text(ell.xra[0]+0.22,0.17,'sources',col='red',/data,target=[prllnon],font_size=12,font_name='Times',font_style='Bold')
    ;; NH lines + text
    t = text(ell.xra[0]+0.2,ell.yra[0]+0.4,'$!8N!7_H  [cm^{-2}]$',col='black',/data,target=[prllnon],font_size=11,font_name='Times')
    for i = 0,n_elements(ll_lines)-1 do begin
        p = plot(nh_lines,ll_lines[i]*[1.,1.],'--',col='black',thick=2,/ov)
        t = text(ell.xra[0]+0.2,ypos[i],nhtext[i],col='black',/data,target=[prllnon],alignment=0.,font_size=11,font_name='Times')
    endfor
    
    ;; remove axes
    prlldet.axes[1].showtext=0
    prlldet.axes[0].showtext=0
    
    ;; legend
    ;l = legend(target=[prdet,prnon],position=[0.865,0.81],/normal,/auto_text_color,sample_width=0.,horizontal_spacing=0.06,font_name='Times')
    ;l.position = [0.956,0.874]
    
    ;xt = text(0.5,0.04,'$log  !8E!7(!8B-V!7)$',alignment=0.5,font_size=14,font_name='Times')
    yt = text(0.02,0.75,'$log  !8L!7_X / !8L!7_X(!8L!7_{IR})  (2$-$10 keV)$',orientation=90,alignment=0.5,font_size=14,font_name='Times')        
    yt.position = [0.007524362,0.59069969,0.026571130,0.85248213]

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; PLOT NH VS E(B-V)
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    nhxdet = rl2nh(lldet,model=model)
    nhxnon = rl2nh(llnon,model=model)
    
    ;; [unobscured,Compton-thick]     
    nh_bound = alog10([1e21,1.5e24])

    ;; WISE AGN Catalog sources
    pw = plot(lebv[where(iifinal_non and iidet_wac)],llnon[where(iifinal_non and iidet_wac)],_extra=enh,position=[75,65,585,445],current=1,/device,/nodata)
    ;; CT shading
    pw = plot(enh.xra,nh_bound[1]*[1,1],linestyle='',fill_background=1,fill_level=enh.yra[1],fill_color='light grey',fill_transparency=80,/ov)
    ;; plot data
    pwnhnon = arrow(transpose([[lebv[where(iifinal_non and iidet_wac)]],[lebv[where(iifinal_non and iidet_wac)]]]),transpose([[nhxnon[where(iifinal_non and iidet_wac)]],[nhxnon[where(iifinal_non and iidet_wac)]+0.12]]),color="orange",fill_transparency=75,target=[pw],/data,head_size=0.3,clip=1,name=" X-ray non-det.")
    ;pwnhnon = plot(lebv[where(iifinal_non and iidet_wac)],nhxnon[where(iifinal_non and iidet_wac)],"o",sym_filled=1,sym_size=0.5,color="orange",transparency=50,/ov,name=" X-ray non-det.")
    pwnhdet = plot(lebv[where(iifinal_det and iidet_wac)],nhxdet[where(iifinal_det and iidet_wac)],"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov,name=" X-ray detected")
    p = plot(lebv[where(iifinal_det and iidet_wac)],nhxdet[where(iifinal_det and iidet_wac)],"S",sym_filled=0,sym_size=1.5,transparency=85,/ov)
    ;; non-detections off plot range
    ioff = where(iifinal_non and iidet_wac and nhxnon ge enh.yra[1],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[1]-0.12)],[make_array(noff,value=enh.yra[1])-0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ioff = where(iifinal_non and iidet_wac and nhxnon le enh.yra[0],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[0]+0.12)],[make_array(noff,value=enh.yra[0])+0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; detections off plot range
    ioff = where(iifinal_det and iidet_wac and nhxdet gt enh.yra[1],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[1]-0.12)],[make_array(noff,value=enh.yra[1])-0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ioff = where(iifinal_det and iidet_wac and nhxdet lt enh.yra[0],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[0]+0.12)],[make_array(noff,value=enh.yra[0])+0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; CT lines
    p = plot(enh.xra[0]+[0.,0.7],nh_bound[1]*[1,1],'-',thick=2,/ov)
    ;; catalog sources text
    t = text(enh.xra[0]+0.22,nh_bound[1]+0.24,'Compton',col='black',target=pwnhnon,font_size=12,font_name='Times',font_style='Bold',/data)
    t = text(enh.xra[0]+0.22,nh_bound[1]+0.09,'thick',col='black',target=pwnhnon,font_size=12,font_name='Times',font_style='Bold',/data)

    ;; Remaining sources
    pr = plot(lebv[where(iifinal_non and ~iidet_wac)],nhxnon[where(iifinal_non and ~iidet_wac)],_extra=enh,position=[585,65,1095,445],current=1,/device,/nodata)
    ;; CT shading
    pr = plot(enh.xra,nh_bound[1]*[1,1],linestyle='',fill_background=1,fill_level=enh.yra[1],fill_color='light grey',fill_transparency=80,/ov)
    ;; plot data
    prnhnon = arrow(transpose([[lebv[where(iifinal_non and ~iidet_wac)]],[lebv[where(iifinal_non and ~iidet_wac)]]]),transpose([[nhxnon[where(iifinal_non and ~iidet_wac)]],[nhxnon[where(iifinal_non and ~iidet_wac)]+0.12]]),color="orange",fill_transparency=75,target=[pr],/data,head_size=0.3,clip=2,name=" X-ray non-det.")
    ;prnhnon = plot(lebv[where(iifinal_non and ~iidet_wac)],nhxnon[where(iifinal_non and ~iidet_wac)],"o",sym_filled=1,sym_size=0.5,color="orange",transparency=50,/ov,name=" X-ray non-det.")
    prnhdet = plot(lebv[where(iifinal_det and ~iidet_wac)],nhxdet[where(iifinal_det and ~iidet_wac)],"S",sym_filled=1,sym_size=1.,color="dodger blue",/ov,name=" X-ray detected")
    p = plot(lebv[where(iifinal_det and ~iidet_wac)],nhxdet[where(iifinal_det and ~iidet_wac)],"S",sym_filled=1,sym_size=1.5,transparency=85,/ov)
    ;; non-detections off plot range
    ioff = where(iifinal_non and ~iidet_wac and nhxnon ge enh.yra[1],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[1]-0.12)],[make_array(noff,value=enh.yra[1])-0.02]]),color='orange',/ov,thick=2,/data,fill_transparency=95,head_size=0.8,target=p)
    ioff = where(iifinal_non and ~iidet_wac and nhxnon le enh.yra[0],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[0]+0.12)],[make_array(noff,value=enh.yra[0])+0.02]]),color='orange',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; detections off plot range
    ioff = where(iifinal_det and ~iidet_wac and nhxdet gt enh.yra[1],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[1]-0.12)],[make_array(noff,value=enh.yra[1])-0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ioff = where(iifinal_det and ~iidet_wac and nhxdet lt enh.yra[0],noff)
    if (noff gt 0) then pwoff = arrow(transpose([[lebv[ioff]],[lebv[ioff]]]),transpose([[make_array(noff,value=enh.yra[0]+0.12)],[make_array(noff,value=enh.yra[0])+0.02]]),color='dodger blue',/ov,/data,fill_transparency=100,head_size=0.8,target=p)
    ;; CT lines
    p = plot(enh.xra[0]+[0.,0.7],nh_bound[1]*[1,1],'-',thick=2,/ov)
    ;; catalog sources text
    ;t = text(e.xra[0]+0.2,nh_bound[1]+0.2,'Compton',col='black',target=pwnon,font_size=12,font_name='Times',font_style='Bold',/data)
    ;t = text(e.xra[0]+0.2,nh_bound[1]+0.05,'thick',col='black',target=pwnon,font_size=12,font_name='Times',font_style='Bold',/data)
    ;; remove axes
    prnhdet.axes[1].showtext=0    
    ;; legend
    l = legend(target=[prnhdet,prnhnon],/normal,/auto_text_color,sample_width=0.,horizontal_spacing=0.06,font_name='Times')
    l.position = [0.69,0.49]
    
    xt = text(0.52,0.02,'$log  !8E!7(!8B-V!7)$',alignment=0.5,font_size=14,font_name='Times')
    yt = text(0.02,0.25,'$log  !8N!7_H [cm^{-2}]$',orientation=90,alignment=0.5,font_size=14,font_name='Times')        
    yt.position = [0.007524362,0.22634560,0.026571130,0.35319986]
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/rlum_a.eps',/BITMAP else $
                                                     p.save,'figures/rlum_a.png';,resolution=20                             
    endif
endif



;;----------------------------------------------------------------------------------------
;; NH Distribution
;;----------------------------------------------------------------------------------------
if keyword_set(nh_dist) then begin
    print, '    PREPARING NH DISTRIBUTION'
    
    normalize = 1   
    e = {xra:[21.,25.],yra:[0.,100.], $
         stairstep:1,fill_background:1, $
         dimension:[1130,880], $
         font_name:'Times', $
         buffer:1}
    if keyword_set(show) then e.buffer = 0    
    if (normalize eq 1) then e.yra = [0.,1.2]
    
    ;; BORUS MODEL
    ;; shaded plot boundaries
    ;; [unobscured,Compton-thick]
    model = 'BORUS'     
    nh_bound = alog10([1e21,1.5e24])
    ll_bound = rl2nh(nh_bound,model=model,/lum_out)
    
    ;; WISE AGN Catalog sources
    ;; sources in plot
    ;iiwd_borus = xhist_wdet_borus gt e.xra[0] and xhist_wdet_borus lt e.xra[1]
    ;iiwn_borus = xhist_wnon_borus gt e.xra[0] and xhist_wnon_borus lt e.xra[1]
    ;iwd_borus = where(iiwd_borus)
    ;iwn_borus = where(iiwn_borus)
    ;; number of sources in plot
    ctwd_borus = commas(fix(total(yhist_wdet_borus)))
    ctwn_borus = commas(fix(total(yhist_wnon_borus)))
    ;ctwdlo_borus = fix(total(yhist_wdet_borus[where(~iiwd_borus and xhist_wdet_borus le e.xra[0])]))
    ;ctwnlo_borus = fix(total(yhist_wnon_borus[where(~iiwn_borus and xhist_wnon_borus le e.xra[0])]))
    ;ctwdhi_borus = fix(total(yhist_wdet_borus[where(~iiwd_borus and xhist_wdet_borus ge e.xra[1])]))
    ;ctwnhi_borus = fix(total(yhist_wnon_borus[where(~iiwn_borus and xhist_wnon_borus ge e.xra[1])]))
    ;ctwlo_borus = ctwdlo_borus+ctwnlo_borus
    ;ctwhi_borus = ctwdhi_borus+ctwnhi_borus

    ;; make plot
    pw = plot(xhist_wnon_borus,yhist_wnon_borus,_extra=e,position=[75,445,585,825],/device,/nodata)
    ;pw.xtickvalues = [ceil(e.xra[0]):floor(e.xra[1])]
    pw.axes[0].showtext=0
    if (normalize eq 1) then pw.ytickvalues = [0.0:1.2:0.2]
    ;; CT shading
    p = plot([alog10(1.5e24),e.xra[1]],e.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=e.yra[0],fill_color='light grey',fill_transparency=80,/ov)
    ;; data
    if (normalize eq 1) then begin
        hwnon = plot(xhist_wnon_borus,nm(yhist_wnon_borus),_extra=e,col='orange',thick=2,linestyle='__',fill_color='orange',fill_transparency=75,/ov,name=' X-ray non-det.')
        hwdet = plot(xhist_wdet_borus,nm(yhist_wdet_borus),_extra=e,col='dodger blue',thick=2,fill_color='dodger blue',fill_transparency=75,/ov,name=' X-ray detected')
        hwlim = arrow(transpose([[intarr(4)+alog10(1.5e24)+0.1],[intarr(4)+alog10(1.5e24)+0.3]]),transpose([[findgen(4,start=5)/10.],[findgen(4,start=5)/10.]]),color='orange',/ov,/data,thick=4,fill_transparency=75,head_size=0.8,target=hwnon)
    endif else begin
        hwnon = plot(xhist_wnon_borus,yhist_wnon_borus,_extra=e,col='orange',thick=2,linestyle='__',fill_color='orange',fill_transparency=75,/ov,name=' X-ray non-det.')
        hwdet = plot(xhist_wdet_borus,yhist_wdet_borus,_extra=e,col='dodger blue',thick=2,fill_color='dodger blue',fill_transparency=75,/ov,name=' X-ray detected')
        hwlim = arrow(transpose([[intarr(4)+alog10(1.5e24)+0.1],[intarr(4)+alog10(1.5e24)+0.3]]),transpose([[e.yra[1]/findgen(4,start=2)],[e.yra[1]/findgen(4,start=2)]]),color='orange',/ov,/data,thick=4,fill_transparency=75,head_size=0.8,target=hwnon)
    endelse
    ;; CT lines
    p = plot([1,1]*alog10(1.5e24),[e.yra[1]*0.875,e.yra[1]],'-',thick=2,/ov)
    p = plot([1,1]*alog10(1.5e24),[e.yra[0],e.yra[1]*0.125],'-',thick=2,/ov)
    ;ct = text(24.22,e.yra[1]*0.902,'Compton',col='black',target=pw,font_name='Times',font_style='Bold',/data)            
    ;ct = text(24.22,e.yra[1]*0.86,'thick',col='black',target=pw,font_name='Times',font_style='Bold',/data)            

    ;; print sources in plot range
    t = text(e.xra[0]+1.0,0.75,'X-ray det:',col='dodger blue',target=pw,/data,alignment=0.5,font_name='Times',font_style='Bold')
    t = text(e.xra[0]+1.0,0.70,strtrim(ctwd_borus,2)+' sources',col='dodger blue',target=pw,/data,alignment=0.5,font_name='Times',font_style='Bold')
    t = text(e.xra[0]+1.0,0.60,'X-ray non-det:',col='orange',target=pw,/data,alignment=0.5,font_name='Times',font_style='Bold')
    t = text(e.xra[0]+1.0,0.55,strtrim(ctwn_borus,2)+' sources',col='orange',target=pw,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;; print s outside plode plot range
    ;; exit stage left
    ;t = text(e.xra[0]+0.4,00.27,strtrim(ctwlo_borus,2),col='dark grey',target=pw,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[0]+0.4,00.22,'sources',col='dark grey',target=pw,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[0]+0.4,00.12,'$\Leftarrow$',target=pw,/data,alignment=0.5,col='dark grey',font_size=32)
    ;; exit stage right
    ;t = text(e.xra[1]-0.4,0.5,strtrim(ctwhi,2),col='dark grey',target=pw,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[1]-0.4,0.45,'sources',col='dark grey',target=pw,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[1]-0.4,0.35,'$\Rightarrow$',target=pw,/data,alignment=0.5,col='dark grey',font_size=32)
    ;; Catalog sources
    t = text(e.xra[0]+0.14,e.yra[1]*0.902,'!16WISE!15 AGN',col='red',target=pw,/data,font_size=12,font_name='Times',font_style='Bold')
    t = text(e.xra[0]+0.14,e.yra[1]*0.86,'sources',col='red',target=pw,/data,font_size=12,font_name='Times',font_style='Bold')
    xmod = text(e.xra[0]+diff(e.xra)/3.,e.yra[1]*0.89,'BORUS model',target=pw,/data,font_size=14,font_name='Times',font_style='Bold')

    ;; Remaining sources
    ;; sources in plot
    ;iird_borus = xhist_rdet_borus gt 21. and xhist_rdet_borus lt 25.8
    ;iirn_borus = xhist_rnon_borus gt 21. and xhist_rnon_borus lt 25.8
    ;ird_borus = where(iird_borus)
    ;irn_borus = where(iirn_borus)
    ;; number of sources in plot
    ctrd_borus = commas(fix(total(yhist_rdet_borus)))
    ctrn_borus = commas(fix(total(yhist_rnon_borus)))
    ;ctrdlo_borus = fix(total(yhist_rdet_borus[where(~iird_borus and xhist_rdet_borus le e.xra[0])]))
    ;ctrnlo_borus = fix(total(yhist_rnon_borus[where(~iirn_borus and xhist_rnon_borus le e.xra[0])]))
    ;ctrdhi_borus = fix(total(yhist_rdet_borus[where(~iird_borus and xhist_rdet_borus ge e.xra[1])]))
    ;ctrnhi_borus = fix(total(yhist_rnon_borus[where(~iirn_borus and xhist_rnon_borus ge e.xra[1])]))
    ;ctrlo_borus = ctrdlo_borus+ctrnlo_borus
    ;ctrhi_borus = ctrdhi_borus+ctrnhi_borus
    
    ;; make plot
    pr = plot(xhist_rnon_borus,yhist_rnon_borus,_extra=e,current=1,position=[585,445,1095,825],/device,/nodata)
    pr.axes[0].showtext=0
    pr.axes[1].showtext=0
    ;pr.xtickvalues = [22.:24.]
    ;; CT shading
    p = plot([alog10(1.5e24),e.xra[1]],e.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=e.yra[0],fill_color='light grey',fill_transparency=80,/ov)
    ;; data
    if (normalize eq 1) then begin
        hrnon = plot(xhist_rnon_borus,nm(yhist_rnon_borus),_extra=e,col='orange',thick=2,linestyle='__',fill_color='orange',fill_transparency=75,/ov,name=' X-ray non-det.')
        hrdet = plot(xhist_rdet_borus,nm(yhist_rdet_borus),_extra=e,col='dodger blue',thick=2,fill_color='dodger blue',fill_transparency=75,/ov,name=' X-ray detected')
        hrlim = arrow(transpose([[intarr(4)+alog10(1.5e24)+0.1],[intarr(4)+alog10(1.5e24)+0.3]]),transpose([[findgen(4,start=5)/10.],[findgen(4,start=5)/10.]]),color='orange',/ov,/data,thick=4,fill_transparency=75,head_size=0.8,target=hrnon)
    endif else begin
        hrnon = plot(xhist_rnon_borus,yhist_rnon_borus,_extra=e,col='orange',thick=2,linestyle='__',fill_color='orange',fill_transparency=75,/ov,name=' X-ray non-det.')
        hrdet = plot(xhist_rdet_borus,yhist_rdet_borus,_extra=e,col='dodger blue',thick=2,fill_color='dodger blue',fill_transparency=75,/ov,name=' X-ray detected')
        hrlim = arrow(transpose([[intarr(4)+alog10(1.5e24)+0.1],[intarr(4)+alog10(1.5e24)+0.3]]),transpose([[e.yra[1]/findgen(4,start=2)],[e.yra[1]/findgen(4,start=2)]]),color='orange',/ov,/data,thick=4,fill_transparency=75,head_size=0.8,target=hrnon)
    endelse
    ;; CT lines
    p = plot([1,1]*alog10(1.5e24),[e.yra[1]*0.875,e.yra[1]],'-',thick=2,/ov)
    p = plot([1,1]*alog10(1.5e24),[e.yra[0],e.yra[1]*0.125],'-',thick=2,/ov)
    ;ct = text(24.25,1.10,'Compton',col='black',target=pr,font_name='Times',font_style='Bold',/data)            
    ;ct = text(24.25,1.05,'thick',col='black',target=pr,font_name='Times',font_style='Bold',/data)            

    ;; sources in plot range
    t = text(e.xra[0]+1.0,0.75,'X-ray det:',col='dodger blue',target=pr,/data,alignment=0.5,font_name='Times',font_style='Bold')
    t = text(e.xra[0]+1.0,0.7,strtrim(ctrd_borus,2)+' sources',col='dodger blue',target=pr,/data,alignment=0.5,font_name='Times',font_style='Bold')
    t = text(e.xra[0]+1.0,0.60,'X-ray non-det:',col='orange',target=pr,/data,alignment=0.5,font_name='Times',font_style='Bold')
    t = text(e.xra[0]+1.0,0.55,strtrim(ctrn_borus,2)+' sources',col='orange',target=pr,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;; sources outside plot range
    ;; exit stage left
    ;t = text(e.xra[0]+0.4,0.27,strtrim(ctrlo_borus,2),col='dark grey',target=pr,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[0]+0.4,0.22,'sources',col='dark grey',target=pr,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[0]+0.4,0.12,'$\Leftarrow$',target=pr,/data,col='dark grey',alignment=0.5,font_size=32)
    ;; exit stage right
    ;t = text(e.xra[1]-0.3,0.5,strtrim(ctrhi,2),col='dark grey',target=pr,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[1]-0.3,0.45,'sources',col='dark grey',target=pr,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[1]-0.3,0.35,'$\Rightarrow$',target=pr,/data,alignment=0.5,col='dark grey',font_size=32)
    ;; Catalog sources
    t = text(e.xra[0]+0.14,e.yra[1]*0.902,'Remaining',col='red',target=pr,/data,font_size=12,font_name='Times',font_style='Bold')
    t = text(e.xra[0]+0.14,e.yra[1]*0.86,'sources',col='red',target=pr,/data,font_size=12,font_name='Times',font_style='Bold')

    xt = text(0.52,0.03,'$log  !8N!7_{H} [cm^{-2}]$',alignment=0.5,vertical_alignment=0.5,font_size=14,font_name='Times')            
    yt = text(0.015,0.5,'Frequency [normalized]',orientation=90.,alignment=0.5,vertical_alignment=0.5,font_size=14,font_name='Times')
    
    ;; POWER-LAW MODEL
    ;; shaded plot boundaries
    ;; [unobscured,Compton-thick]
    model = 'POWER'     
    nh_bound = alog10([1e21,1.5e24])
    ll_bound = rl2nh(nh_bound,model=model,/lum_out)

    ;; WISE AGN Catalog sources
    ;; sources in plot
    ;iiwd_power = xhist_wdet_power gt e.xra[0] and xhist_wdet_power lt e.xra[1]
    ;iiwn_power = xhist_wnon_power gt e.xra[0] and xhist_wnon_power lt e.xra[1]
    ;iwd_power = where(iiwd_power)
    ;iwn_power = where(iiwn_power)
    ;; number of sources in plot
    ;ctwd_power = commas(fix(total(yhist_wdet_power)))
    ;ctwn_power = commas(fix(total(yhist_wnon_power)))
    ;ctwdlo_power = fix(total(yhist_wdet_power[where(~iiwd_power and xhist_wdet_power le e.xra[0])]))
    ;ctwnlo_power = fix(total(yhist_wnon_power[where(~iiwn_power and xhist_wnon_power le e.xra[0])]))
    ;ctwdhi_power = fix(total(yhist_wdet_power[where(~iiwd_power and xhist_wdet_power ge e.xra[1])]))
    ;ctwnhi_power = fix(total(yhist_wnon_power[where(~iiwn_power and xhist_wnon_power ge e.xra[1])]))
    ;ctwlo_power = ctwdlo_power+ctwnlo_power
    ;ctwhi_power = ctwdhi_power+ctwnhi_power

    ;; make plot
    pw = plot(xhist_wnon_power,yhist_wnon_power,_extra=e,current=1,position=[75,65,585,445],/device,/nodata)
    pw.xtickvalues = [21.:24.:1.]
    if (normalize eq 1) then pw.ytickvalues = [0.0:1.:0.2]
    ;; CT shading
    p = plot([alog10(1.5e24),e.xra[1]],e.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=e.yra[0],fill_color='light grey',fill_transparency=80,/ov)
    ;; data
    if (normalize eq 1) then begin
        hwnon = plot(xhist_wnon_power,nm(yhist_wnon_power),_extra=e,col='orange',thick=2,linestyle='__',fill_color='orange',fill_transparency=75,/ov,name=' X-ray non-det.')
        hwdet = plot(xhist_wdet_power,nm(yhist_wdet_power),_extra=e,col='dodger blue',thick=2,fill_color='dodger blue',fill_transparency=75,/ov,name=' X-ray detected')
        hwlim = arrow(transpose([[intarr(4)+alog10(1.5e24)+0.1],[intarr(4)+alog10(1.5e24)+0.3]]),transpose([[findgen(4,start=5)/10.],[findgen(4,start=5)/10.]]),color='orange',/ov,/data,thick=4,fill_transparency=75,head_size=0.8,target=hwnon)
    endif else begin
        hwnon = plot(xhist_wnon_power,yhist_wnon_power,_extra=e,col='orange',thick=2,linestyle='__',fill_color='orange',fill_transparency=75,/ov,name=' X-ray non-det.')
        hwdet = plot(xhist_wdet_power,yhist_wdet_power,_extra=e,col='dodger blue',thick=2,fill_color='dodger blue',fill_transparency=75,/ov,name=' X-ray detected')
        hwlim = arrow(transpose([[intarr(4)+alog10(1.5e24)+0.1],[intarr(4)+alog10(1.5e24)+0.3]]),transpose([[e.yra[1]/findgen(4,start=2)],[e.yra[1]/findgen(4,start=2)]]),color='orange',/ov,/data,thick=4,fill_transparency=75,head_size=0.8,target=hwnon)
    endelse 
    ;; CT lines
    p = plot([1,1]*alog10(1.5e24),[e.yra[1]*0.875,e.yra[1]],'-',thick=2,/ov)
    p = plot([1,1]*alog10(1.5e24),[e.yra[0],e.yra[1]*0.125],'-',thick=2,/ov)
    ;ct = text(24.22,e.yra[1]*0.922,'Compton',col='black',target=pw,font_name='Times',font_style='Bold',/data)            
    ;ct = text(24.22,e.yra[1]*0.88,'thick',col='black',target=pw,font_name='Times',font_style='Bold',/data)            

    ;; print sources in plot range
    ;t = text(e.xra[0]+1.0,0.75,'X-ray det:',col='dodger blue',target=pw,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[0]+1.0,0.7,strtrim(ctwd_power,2)+' sources',col='dodger blue',target=pw,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[0]+1.0,0.60,'X-ray non-det:',col='orange',target=pw,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[0]+1.0,0.55,strtrim(ctwn_power,2)+' sources',col='orange',target=pw,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;; print s outside plode plot range
    ;; exit stage left
    ;t = text(e.xra[0]+0.4,00.27,strtrim(ctwlo_power,2),col='dark grey',target=pw,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[0]+0.4,00.22,'sources',col='dark grey',target=pw,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[0]+0.4,00.12,'$\Leftarrow$',target=pw,/data,alignment=0.5,col='dark grey',font_size=32)
    ;; exit stage right
    ;t = text(e.xra[1]-0.4,0.5,strtrim(ctwhi,2),col='dark grey',target=pw,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[1]-0.4,0.45,'sources',col='dark grey',target=pw,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[1]-0.4,0.35,'$\Rightarrow$',target=pw,/data,alignment=0.5,col='dark grey',font_size=32)
    ;; Catalog sources
    ;t = text(e.xra[0]+0.15,0.24,'!16WISE!15 AGN',col='red',target=pw,/data,font_size=12,font_name='Times',font_style='Bold')
    ;t = text(e.xra[0]+0.15,0.18,'sources',col='red',target=pw,/data,font_size=12,font_name='Times',font_style='Bold')
    xmod = text(e.xra[0]+diff(e.xra)/3.,e.yra[1]*0.89,'Power law model',target=pw,/data,font_size=14,font_name='Times',font_style='Bold')
    
    ;; Remaining sources
    ;; sources in plot
    ;iird_power = xhist_rdet_power gt e.xra[0] and xhist_rdet_power lt e.xra[1]
    ;iirn_power = xhist_rnon_power gt e.xra[0] and xhist_rnon_power lt e.xra[1]
    ;ird_power = where(iird_power)
    ;irn_power = where(iirn_power)
    ;; number of sources in plot
    ctrd_power = commas(fix(total(yhist_rdet_power)))
    ctrn_power = commas(fix(total(yhist_rnon_power)))
    ;ctrdlo_power = fix(total(yhist_rdet_power[where(~iird_power and xhist_rdet_power le e.xra[0])]))
    ;ctrnlo_power = fix(total(yhist_rnon_power[where(~iirn_power and xhist_rnon_power le e.xra[0])]))
    ;ctrdhi_power = fix(total(yhist_rdet_power[where(~iird_power and xhist_rdet_power ge e.xra[1])]))
    ;ctrnhi_power = fix(total(yhist_rnon_power[where(~iirn_power and xhist_rnon_power ge e.xra[1])]))
    ;ctrlo_power = ctrdlo_power+ctrnlo_power
    ;ctrhi_power = ctrdhi_power+ctrnhi_power
    
    ;; make plot
    pr = plot(xhist_rnon_power,yhist_rnon_power,_extra=e,current=1,position=[585,65,1095,445],/device,/nodata)
    pr.xtickvalues = [22.:25.:1.]
    pr.axes[1].showtext=0
    ;; CT shading
    p = plot([alog10(1.5e24),e.xra[1]],e.yra[1]*[1.,1.],linestyle='',fill_background=1,fill_level=e.yra[0],fill_color='light grey',fill_transparency=80,/ov)
    ;; data
    if (normalize eq 1) then begin
        hrnon = plot(xhist_rnon_power,nm(yhist_rnon_power),_extra=e,col='orange',thick=2,linestyle='__',fill_color='orange',fill_transparency=75,/ov,name=' X-ray non-det.')
        hrdet = plot(xhist_rdet_power,nm(yhist_rdet_power),_extra=e,col='dodger blue',thick=2,fill_color='dodger blue',fill_transparency=75,/ov,name=' X-ray detected')
        hrlim = arrow(transpose([[intarr(4)+alog10(1.5e24)+0.1],[intarr(4)+alog10(1.5e24)+0.3]]),transpose([[findgen(4,start=5)/10.],[findgen(4,start=5)/10.]]),color='orange',/ov,/data,thick=4,fill_transparency=75,head_size=0.8,target=hrnon)
    endif else begin
        hrnon = plot(xhist_rnon_power,yhist_rnon_power,_extra=e,col='orange',thick=2,linestyle='__',fill_color='orange',fill_transparency=75,/ov,name=' X-ray non-det.')
        hrdet = plot(xhist_rdet_power,yhist_rdet_power,_extra=e,col='dodger blue',thick=2,fill_color='dodger blue',fill_transparency=75,/ov,name=' X-ray detected')
        hrlim = arrow(transpose([[intarr(4)+alog10(1.5e24)+0.1],[intarr(4)+alog10(1.5e24)+0.3]]),transpose([[e.yra[1]/findgen(4,start=2)],[e.yra[1]/findgen(4,start=2)]]),color='orange',/ov,/data,thick=4,fill_transparency=75,head_size=0.8,target=hrnon)
    endelse
    ;; CT lines
    p = plot([1,1]*alog10(1.5e24),[e.yra[1]*0.875,e.yra[1]],'-',thick=2,/ov)
    p = plot([1,1]*alog10(1.5e24),[e.yra[0],e.yra[1]*0.125],'-',thick=2,/ov)
    ct = text(24.22,e.yra[1]*0.095,'Compton',col='black',target=pr,font_name='Times',font_style='Bold',/data)          
    ct = text(24.22,e.yra[1]*0.053,'thick',col='black',target=pr,font_name='Times',font_style='Bold',/data)            

    ;;sources in plot range
    ;t = text(e.xra[0]+1.0,0.75,'X-ray det:',col='dodger blue',target=pr,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[0]+1.0,0.7,strtrim(ctrd_power,2)+' sources',col='dodger blue',target=pr,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[0]+1.0,0.60,'X-ray non-det:',col='orange',target=pr,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[0]+1.0,0.55,strtrim(ctrn_power,2)+' sources',col='orange',target=pr,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;; sources outside plot range
    ;; exit stage left
    ;t = text(e.xra[0]+0.4,0.27,strtrim(ctrlo_power,2),col='dark grey',target=pr,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[0]+0.4,0.22,'sources',col='dark grey',target=pr,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[0]+0.4,0.12,'$\Leftarrow$',target=pr,/data,col='dark grey',alignment=0.5,font_size=32)
    ;; exit stage right
    ;t = text(e.xra[1]-0.3,0.5,strtrim(ctrhi,2),col='dark grey',target=pr,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[1]-0.3,0.45,'sources',col='dark grey',target=pr,/data,alignment=0.5,font_name='Times',font_style='Bold')
    ;t = text(e.xra[1]-0.3,0.35,'$\Rightarrow$',target=pr,/data,alignment=0.5,col='dark grey',font_size=32)
    ;; Catalog sources
    ;t = text(e.xra[0]+0.15,0.24,'Remaining',col='red',target=pr,/data,font_size=12,font_name='Times',font_style='Bold')
    ;t = text(e.xra[0]+0.15,0.18,'sources',col='red',target=pr,/data,font_size=12,font_name='Times',font_style='Bold')

    ;xt = text(0.5,0.05,'$log  !8N!7_{H} [cm^{-2}]$',alignment=0.5,vertical_alignment=0.5,font_size=14,font_name='Times')            
    ;yt = text(0.015,0.5,'Frequency',orientation=90.,alignment=0.5,vertical_alignment=0.5,font_size=14,font_name='Times')

    l = legend(target=[hwdet,hwnon],/normal,/auto_text_color,sample_width=0.1,horizontal_spacing=0.06,font_name='Times')
    l.position = [0.745,0.49]

    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/nh_dist.eps',/BITMAP else $
                                                     p.save,'figures/nh_dist.png';,resolution=20
                             
    endif
endif


END




