PRO surv_analysis, FMT = fmt

;common _fits    
;common _resamp
;common _comp    
;common _inf_cha 
;common _inf_xmm 
;common _inf_nst 
;common _det_cha 
;common _det_xmm 
;common _det_nst 
common _wac 
common _xconv   
;common _fxlim 
common _agnlum 
common _quality  
common _combined
common _nhdist


surv_dir = 'surv_analysis'
if keyword_set(asurv) then $
    if (file_search(surv_dir) ne '') then file_delete,surv_dir,/recursive

file_mkdir,surv_dir
pushd,surv_dir


case strupcase(fmt) of
    'ASURV': begin
        spawn,'ln -s ~/Research/statistics/asurv-master/asurv .'
        spawn,'ln -s /Users/ccarroll/IDLWorkspace/libraries/cmc/asurv_idl/all.com .'
        spawn,'ln -s /Users/ccarroll/IDLWorkspace/libraries/cmc/asurv_idl/wagn.com .'
        spawn,'ln -s /Users/ccarroll/IDLWorkspace/libraries/cmc/asurv_idl/ragn.com .'
        spawn,'ln -s /Users/ccarroll/IDLWorkspace/libraries/cmc/asurv_idl/all2.com .'
        spawn,'ln -s /Users/ccarroll/IDLWorkspace/libraries/cmc/asurv_idl/wagn2.com .'
        spawn,'ln -s /Users/ccarroll/IDLWorkspace/libraries/cmc/asurv_idl/ragn2.com .'
        spawn,'ln -s /Users/ccarroll/IDLWorkspace/libraries/cmc/asurv_idl/all3.com'
        ;; X-ray det/non-det, luminosity ratio, SED obscured/unobscured/borderline
        xd = intarr(nsrc)+99
        xd[where(iiqual_non)] = 1
        xd[where(iiqual_det)] = 0
        rl = dblarr(nsrc)-9999.
        rl[where(iiqual_non)] = llnon[where(iiqual_non)]
        rl[where(iiqual_det)] = lldet[where(iiqual_det)]
        obsc = intarr(nsrc)+99
        obsc[where(ebv+e_ebv lt 0.15)] = 0
        obsc[where(ebv-e_ebv gt 0.15)] = 1
        obsc[where((ebv lt 0.15 and ebv+e_ebv gt 0.15) or (ebv gt 0.15 and ebv-e_ebv lt 0.15))] = 2
        ;; CUT TO RUN SURVIVAL ANALYSIS ON ONLY ANALYSIS SET
        xd = xd[where(iiqual)]
        rl = rl[where(iiqual)]
        obsc = obsc[where(iiqual)]
        iwac = iiwac[where(iiqual)]
        
        lir6 = loglir[where(iiqual)]
        lx210 = lir6*0.
        lx210[where(iiqual_det)] = loglx[where(iiqual_det)]
        lx210[where(iiqual_non)] = loglxir[where(iiqual_non)]
        
        ;; limit RL boundary values
        ;rl[where(rl gt 0.,/null)] = 0.
        ;rl[where(rl lt -2.,/null)] = -2.

        ;; remove outliers in detections
        ;iikeep = bytarr(n_elements(rl))
        ;iikeep[where(xd eq -1,/null)] = 1
        ;resistant_mean,rl[where(xd)],3.,mn,sigmn,nrej,goodvec=ig
        ;iikeep[(where(xd eq 0))[ig]] = 1

        ;; RL starts with a detection
        ;iikeep = rl gt min(rl[where(xd eq 0)])

        ;; RL ³ -1
        ;iikeep = rl ge -1

        ;; sparse non-detections
        ;iikeep = bytarr(n_elements(rl))
        ;ixd = where(xd eq 0,nxd,complement=ixn,ncomplement=nxn,/null)
        ;if (nxd gt 0) then iikeep[ixd] = 1
        ;if (nxn gt 0) then iikeep[ixn[randomi(nxd/2.,nxn)]] = 1

        ;; cut on above subsets
        ;ikeep = where(iikeep,ct)
        ;if (ct eq 0) then stop
        ;xd = xd[ikeep]
        ;rl = rl[ikeep]
        ;obsc = obsc[ikeep]
        ;iwac = iwac[ikeep]

        ;; WISE AGN and remaning
        xdw = xd[where(iwac)]
        rlw = rl[where(iwac)]
        obscw = obsc[where(iwac)]
        xdr = xd[where(~iwac)]
        rlr = rl[where(~iwac)]
        obscr = obsc[where(~iwac)]
        
        ;; save input sources; WISE AGN/remaining may be subset
        save,xd,rl,obsc,iwac, $
             xdw,rlw,obscw, $
             xdr,rlr,obscr, $
             file='../asurv_input.sav'

        ;; univariate, KM test
        forprint,xd,(-1.)*rl,format='2x,i2,2x,d7.4',textout='all.dat',/nocomment
        forprint,xdw,(-1.)*rlw,format='2x,i2,2x,d7.4',textout='wagn.dat',/nocomment
        forprint,xdr,(-1.)*rlr,format='2x,i2,2x,d7.4',textout='ragn.dat',/nocomment

        ;; univariate, two-sample test [0 remaining, 1 WISEAGN]
        forprint,iwac,xd,(-1.)*rl,format='2x,i2,2x,i2,2x,d7.4',textout='all2.dat',/nocomment

        ;; univariate, two-sample test [0 unobscured, 1 obscured, 2 unsure]
        ;forprint,obsc,xd,(-1.)*rl,format='2x,a2,2x,a2,2x,d7.4',textout='all2.dat',/nocomment
        forprint,obscw,xdw,(-1.)*rlw,format='2x,i2,2x,i2,2x,d7.4',textout='wagn2.dat',/nocomment
        forprint,obscr,xdr,(-1.)*rlr,format='2x,i2,2x,i2,2x,d7.4',textout='ragn2.dat',/nocomment

        ;; bivariate tests
        forprint,xd,lir6,(-1.)*rl,format='2x,i2,2x,d4.1,2x,d5.2',textout='all3.dat',/nocomment

        spawn,'./asurv'
    
        ;; pull bin size
        readcol,'wagn.com',col,format='a'
        svbinsz = double(col[9])

        ;; univariate, KM
        ;; ALL SOURCES
        readcol,'all.out',abinlo,abinhi,akm,aerr,format='x,d,x,d,d',/silent
        ;; add first line
        abinhi = [abinlo[0],abinhi]
        abinlo = [min(rl),abinlo]
        akm = [1.,akm]
        aerr = [0.,aerr]
        abin = abinlo-(abinlo-abinhi)/2.
        readcol,'all.out',abinc,ad,format='d,d',/silent
        ;; remove bad line
        abinc = abinc[1:-1]
        ad = ad[1:-1]
        ;; WISE AGN SOURCES
        readcol,'wagn.out',wbinlo,wbinhi,wkm,werr,format='x,d,x,d,d',/silent
        ;; add first line
        wbinhi = [wbinlo[0],wbinhi]
        wbinlo = [min(rlw),wbinlo]
        wkm = [1.,wkm]
        werr = [0.,werr]
        wbin = wbinlo-(wbinlo-wbinhi)/2.
        readcol,'wagn.out',wbinc,wd,format='d,d',/silent
        ;; remove bad line
        wbinc = wbinc[1:-1]
        wd = wd[1:-1]
        ;; REMAINING SOURCES
        readcol,'ragn.out',rbinlo,rbinhi,rkm,rerr,format='x,d,x,d,d',/silent
        ;; add first line
        rbinhi = [rbinlo[0],rbinhi]
        rbinlo = [min(rlr),rbinlo]
        rkm = [1.,rkm]
        rerr = [0.,rerr]
        rbin = rbinlo-(rbinlo-rbinhi)/2.
        readcol,'ragn.out',rbinc,rd,format='d,d',/silent
        ;; remove bad line
        rbinc = rbinc[1:-1]
        rd = rd[1:-1]

        sav_vars = ['SVBINSZ','ABINLO','ABINHI','AKM','AERR','ABIN','ABINC','AD', $
                              'WBINLO','WBINHI','WKM','WERR','WBIN','WBINC','WD', $
                              'RBINLO','RBINHI','RKM','RERR','RBIN','RBINC','RD']
        sav_inds = []


        svnhbinsz = nhbinsz
        ;svnhbinsz = 0.5
        ;svnhbinsz = svbinsz

        ;; recreate luminosity ratio arrays
        ;; use central luminosity ratio value and distribute randomly around it using histogram frequency
        simrla = []
        simrlw = []
        simrlr = []
        for i = 0,n_elements(wbinc)-1 do begin
            if (ad[i] gt 0.) then simrla = [simrla,abinc[i]+(randomu(seed,round(ad[i]))-0.5)*svbinsz]
            if (wd[i] gt 0.) then simrlw = [simrlw,wbinc[i]+(randomu(seed,round(wd[i]))-0.5)*svbinsz]
            if (rd[i] gt 0.) then simrlr = [simrlr,rbinc[i]+(randomu(seed,round(rd[i]))-0.5)*svbinsz]
        endfor

        anh_power = rl2nh(simrla,model='POWER')
        wnh_power = rl2nh(simrlw,model='POWER')
        rnh_power = rl2nh(simrlr,model='POWER')
        anh_borus = rl2nh(simrla,model='BORUS')
        wnh_borus = rl2nh(simrlw,model='BORUS')
        rnh_borus = rl2nh(simrlr,model='BORUS')

        yanh_power = histogram(anh_power,locations=xanh_power,bin=svnhbinsz)
        ywnh_power = histogram(wnh_power,locations=xwnh_power,bin=svnhbinsz)
        yrnh_power = histogram(rnh_power,locations=xrnh_power,bin=svnhbinsz)
        yanh_borus = histogram(anh_borus,locations=xanh_borus,bin=svnhbinsz)
        ywnh_borus = histogram(wnh_borus,locations=xwnh_borus,bin=svnhbinsz)
        yrnh_borus = histogram(rnh_borus,locations=xrnh_borus,bin=svnhbinsz)

        xanh_power = [xanh_power[0]-svnhbinsz,xanh_power,xanh_power[-1]+svnhbinsz]
        xwnh_power = [xwnh_power[0]-svnhbinsz,xwnh_power,xwnh_power[-1]+svnhbinsz]
        xrnh_power = [xrnh_power[0]-svnhbinsz,xrnh_power,xrnh_power[-1]+svnhbinsz]
        xanh_borus = [xanh_borus[0]-svnhbinsz,xanh_borus,xanh_borus[-1]+svnhbinsz]
        xwnh_borus = [xwnh_borus[0]-svnhbinsz,xwnh_borus,xwnh_borus[-1]+svnhbinsz]
        xrnh_borus = [xrnh_borus[0]-svnhbinsz,xrnh_borus,xrnh_borus[-1]+svnhbinsz]
        yanh_power = [0.,yanh_power,0.]
        ywnh_power = [0.,ywnh_power,0.]
        yrnh_power = [0.,yrnh_power,0.]
        yanh_borus = [0.,yanh_borus,0.]
        ywnh_borus = [0.,ywnh_borus,0.]
        yrnh_borus = [0.,yrnh_borus,0.]

        sav_vars = [sav_vars,'SIMRLA','SIMRLW','SIMRLR', $
                             'ANH_POWER','XANH_POWER','YANH_POWER', $
                             'WNH_POWER','XWNH_POWER','YWNH_POWER', $
                             'RNH_POWER','XRNH_POWER','YRNH_POWER', $
                             'ANH_BORUS','XANH_BORUS','YANH_BORUS', $
                             'WNH_BORUS','XWNH_BORUS','YWNH_BORUS', $
                             'RNH_BORUS','XRNH_BORUS','YRNH_BORUS', $
                             'SVNHBINSZ']
        sav_inds = [sav_inds]


        ;; pull bin size
        readcol,'wagn2.com',col,format='a'
        svbinsz = double(col[11])

        ;; ALL SOURCES
        readcol,'all2.out',abinlo1,abinhi1,akm1,aerr1,format='x,d,x,d,d,d';,/silent,skipline=wskip11,numline=(wnum11-wskip11)
        ;; split samples
        ii = where(diff(abinlo1) lt 0.)
        abinlo2 = abinlo1[ii+1:-1]
        abinlo1 = abinlo1[0:ii]
        abinhi2 = abinhi1[ii+1:-1]
        abinhi1 = abinhi1[0:ii]
        akm2 = akm1[ii+1:-1]
        akm1 = akm1[0:ii]
        aerr2 = aerr1[ii+1:-1]
        aerr1 = aerr1[0:ii]
        ;; correct first lines
        abinhi1 = [abinlo1[0],abinhi1]
        abinlo1 = [min(rl),abinlo1]
        akm1 = [1.,akm1]
        aerr1 = [0.,aerr1]
        abinhi2 = [abinlo2[0],abinhi2]
        abinlo2 = [min(rl),abinlo2]
        akm2 = [1.,akm2]
        aerr2 = [0.,aerr2]
        abin1 = abinlo1-(abinlo1-abinhi1)/2.
        abin2 = abinlo2-(abinlo2-abinhi2)/2.

        readcol,'all2.out',abinc1,ad1,format='d,d',/silent
        ;; remove bad lines
        abinc1 = abinc1[1:-1]
        ad1 = ad1[1:-1]
        ii = (where(diff(abinc1) lt 0))[-1]
        abinc2 = abinc1[ii+1:-1]
        ad2 = ad1[ii+1:-1]
        abinc1 = abinc1[0:ii-1]
        ad1 = ad1[0:ii-1]
        if (total(abinc1-abinc2) ne 0.) then stop

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
        wbinlo1 = [min(rlw),wbinlo1]
        wkm1 = [1.,wkm1]
        werr1 = [0.,werr1]
        wbinhi2 = [wbinlo2[0],wbinhi2]
        wbinlo2 = [min(rlw),wbinlo2]
        wkm2 = [1.,wkm2]
        werr2 = [0.,werr2]
        wbin1 = wbinlo1-(wbinlo1-wbinhi1)/2.
        wbin2 = wbinlo2-(wbinlo2-wbinhi2)/2.

        readcol,'wagn2.out',wbinc1,wd1,format='d,d',/silent
        ;; remove bad lines
        wbinc1 = wbinc1[1:-1]
        wd1 = wd1[1:-1]
        ii = (where(diff(wbinc1) lt 0))[-1]
        wbinc2 = wbinc1[ii+1:-1]
        wd2 = wd1[ii+1:-1]
        wbinc1 = wbinc1[0:ii-1]
        wd1 = wd1[0:ii-1]
        if (total(wbinc1-wbinc2) ne 0.) then stop

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
        rbinlo1 = [min(rlr),rbinlo1]
        rkm1 = [1.,rkm1]
        rerr1 = [0.,rerr1]
        rbinhi2 = [rbinlo2[0],rbinhi2]
        rbinlo2 = [min(rlr),rbinlo2]
        rkm2 = [1.,rkm2]
        rerr2 = [0.,rerr2]
        rbin1 = rbinlo1-(rbinlo1-rbinhi1)/2.
        rbin2 = rbinlo2-(rbinlo2-rbinhi2)/2.

        readcol,'ragn2.out',rbinc1,rd1,format='d,d',/silent
        ;; remove bad lines
        rbinc1 = rbinc1[1:-1]
        rd1 = rd1[1:-1]
        ii = (where(diff(rbinc1) lt 0))[-1]
        rbinc2 = rbinc1[ii+1:-1]
        rd2 = rd1[ii+1:-1]
        rbinc1 = rbinc1[0:ii-1]
        rd1 = rd1[0:ii-1]
        if (total(rbinc1-rbinc2) ne 0.) then stop

        sav_vars = [sav_vars,'ABINLO1','ABINHI1','AKM1','AERR1','ABIN1','ABINC1','AD1', $
                             'ABINLO2','ABINHI2','AKM2','AERR2','ABIN2','ABINC2','AD2', $
                             'WBINLO1','WBINHI1','WKM1','WERR1','WBIN1','WBINC1','WD1', $
                             'WBINLO2','WBINHI2','WKM2','WERR2','WBIN2','WBINC2','WD2', $
                             'RBINLO1','RBINHI1','RKM1','RERR1','RBIN1','RBINC1','RD1', $
                             'RBINLO2','RBINHI2','RKM2','RERR2','RBIN2','RBINC2','RD2']
        sav_inds = [sav_inds]


        ;; recreate luminosity ratio arrays
        ;; use central luminosity ratio value and distribute randomly around it using histogram frequency
        simrla1 = []
        simrla2 = []
        simrlw1 = []
        simrlw2 = []
        simrlr1 = []
        simrlr2 = []
        for i = 0,n_elements(wbinc1)-1 do begin
            if (ad1[i] gt 0.) then simrla1 = [simrla1,abinc1[i]+(randomu(seed,round(ad1[i]))-0.5)*svbinsz]
            if (ad2[i] gt 0.) then simrla2 = [simrla2,abinc2[i]+(randomu(seed,round(ad2[i]))-0.5)*svbinsz]
            if (wd1[i] gt 0.) then simrlw1 = [simrlw1,wbinc1[i]+(randomu(seed,round(wd1[i]))-0.5)*svbinsz]
            if (wd2[i] gt 0.) then simrlw2 = [simrlw2,wbinc2[i]+(randomu(seed,round(wd2[i]))-0.5)*svbinsz]
            if (rd1[i] gt 0.) then simrlr1 = [simrlr1,rbinc1[i]+(randomu(seed,round(rd1[i]))-0.5)*svbinsz]
            if (rd2[i] gt 0.) then simrlr2 = [simrlr2,rbinc2[i]+(randomu(seed,round(rd2[i]))-0.5)*svbinsz]
        endfor

        ;; ALL SOURCES
        anh1_power = rl2nh(simrla1,model='POWER')
        anh2_power = rl2nh(simrla2,model='POWER')
        anh1_borus = rl2nh(simrla1,model='BORUS')
        anh2_borus = rl2nh(simrla2,model='BORUS')

        yanh1_power = histogram(anh1_power,locations=xanh1_power,bin=svnhbinsz)
        yanh2_power = histogram(anh2_power,locations=xanh2_power,bin=svnhbinsz)
        yanh1_borus = histogram(anh1_borus,locations=xanh1_borus,bin=svnhbinsz)
        yanh2_borus = histogram(anh2_borus,locations=xanh2_borus,bin=svnhbinsz)

        xanh1_power = [xanh1_power[0]-svnhbinsz,xanh1_power,xanh1_power[-1]+svnhbinsz]
        xanh2_power = [xanh2_power[0]-svnhbinsz,xanh2_power,xanh2_power[-1]+svnhbinsz]
        yanh1_power = [0.,yanh1_power,0.]
        yanh2_power = [0.,yanh2_power,0.]
        xanh1_borus = [xanh1_borus[0]-svnhbinsz,xanh1_borus,xanh1_borus[-1]+svnhbinsz]
        xanh2_borus = [xanh2_borus[0]-svnhbinsz,xanh2_borus,xanh2_borus[-1]+svnhbinsz]
        yanh1_borus = [0.,yanh1_borus,0.]
        yanh2_borus = [0.,yanh2_borus,0.]

        ;; WISE AGN SOURCES
        wnh1_power = rl2nh(simrlw1,model='POWER')
        wnh2_power = rl2nh(simrlw2,model='POWER')
        wnh1_borus = rl2nh(simrlw1,model='BORUS')
        wnh2_borus = rl2nh(simrlw2,model='BORUS')

        ywnh1_power = histogram(wnh1_power,locations=xwnh1_power,bin=svnhbinsz)
        ywnh2_power = histogram(wnh2_power,locations=xwnh2_power,bin=svnhbinsz)
        ywnh1_borus = histogram(wnh1_borus,locations=xwnh1_borus,bin=svnhbinsz)
        ywnh2_borus = histogram(wnh2_borus,locations=xwnh2_borus,bin=svnhbinsz)

        xwnh1_power = [xwnh1_power[0]-svnhbinsz,xwnh1_power,xwnh1_power[-1]+svnhbinsz]
        xwnh2_power = [xwnh2_power[0]-svnhbinsz,xwnh2_power,xwnh2_power[-1]+svnhbinsz]
        ywnh1_power = [0.,ywnh1_power,0.]
        ywnh2_power = [0.,ywnh2_power,0.]
        xwnh1_borus = [xwnh1_borus[0]-svnhbinsz,xwnh1_borus,xwnh1_borus[-1]+svnhbinsz]
        xwnh2_borus = [xwnh2_borus[0]-svnhbinsz,xwnh2_borus,xwnh2_borus[-1]+svnhbinsz]
        ywnh1_borus = [0.,ywnh1_borus,0.]
        ywnh2_borus = [0.,ywnh2_borus,0.]

        ;; REMAINING SOURCES
        rnh1_power = rl2nh(simrlr1,model='POWER')
        rnh2_power = rl2nh(simrlr2,model='POWER')
        rnh1_borus = rl2nh(simrlr1,model='BORUS')
        rnh2_borus = rl2nh(simrlr2,model='BORUS')

        yrnh1_power = histogram(rnh1_power,locations=xrnh1_power,bin=svnhbinsz)
        yrnh2_power = histogram(rnh2_power,locations=xrnh2_power,bin=svnhbinsz)
        yrnh1_borus = histogram(rnh1_borus,locations=xrnh1_borus,bin=svnhbinsz)
        yrnh2_borus = histogram(rnh2_borus,locations=xrnh2_borus,bin=svnhbinsz)

        xrnh1_power = [xrnh1_power[0]-svnhbinsz,xrnh1_power,xrnh1_power[-1]+svnhbinsz]
        xrnh2_power = [xrnh2_power[0]-svnhbinsz,xrnh2_power,xrnh2_power[-1]+svnhbinsz]
        yrnh1_power = [0.,yrnh1_power,0.]
        yrnh2_power = [0.,yrnh2_power,0.]
        xrnh1_borus = [xrnh1_borus[0]-svnhbinsz,xrnh1_borus,xrnh1_borus[-1]+svnhbinsz]
        xrnh2_borus = [xrnh2_borus[0]-svnhbinsz,xrnh2_borus,xrnh2_borus[-1]+svnhbinsz]
        yrnh1_borus = [0.,yrnh1_borus,0.]
        yrnh2_borus = [0.,yrnh2_borus,0.]

        sav_vars = [sav_vars,'SIMRLA1','ANH1_POWER','XANH1_POWER','YANH1_POWER', $
                             'SIMRLA2','ANH2_POWER','XANH2_POWER','YANH2_POWER', $
                             'SIMRLW1','WNH1_POWER','XWNH1_POWER','YWNH1_POWER', $
                             'SIMRLW2','WNH2_POWER','XWNH2_POWER','YWNH2_POWER', $
                             'SIMRLR1','RNH1_POWER','XRNH1_POWER','YRNH1_POWER', $
                             'SIMRLR2','RNH2_POWER','XRNH2_POWER','YRNH2_POWER', $
                                       'ANH1_BORUS','XANH1_BORUS','YANH1_BORUS', $
                                       'ANH2_BORUS','XANH2_BORUS','YANH2_BORUS', $
                                       'WNH1_BORUS','XWNH1_BORUS','YWNH1_BORUS', $
                                       'WNH2_BORUS','XWNH2_BORUS','YWNH2_BORUS', $
                                       'RNH1_BORUS','XRNH1_BORUS','YRNH1_BORUS', $
                                       'RNH2_BORUS','XRNH2_BORUS','YRNH2_BORUS']
        sav_inds = [sav_inds]
        end
    'R': begin
        ;; CUT TO RUN SURVIVAL ANALYSIS ON ONLY ANALYSIS SET
        iq = where(iiqual,nq)
        xd = iiqual_det[iq]
        rl = (llnon>lldet)[iq]
        ;obsc = intarr(nq)+99
        ;obsc[where((ebv+e_ebv)[iq] le 0.15)] = 1
        ;obsc[where((ebv-e_ebv)[iq] gt 0.15)] = 2
        ;obsc[where((ebv[iq] lt 0.15 and (ebv+e_ebv)[iq] gt 0.15) or (ebv[iq] gt 0.15 and (ebv-e_ebv)[iq] lt 0.15))] = 3
        obsc = (ebv[iq] ge 0.15)+1
        wagn = (iiwac[iq] eq 0)+1
        ;; WISE AGN
        xdw = xd[where(wagn eq 1)]
        rlw = rl[where(wagn eq 1)]
        obscw = obsc[where(wagn eq 1)]
        agnw = wagn[where(wagn eq 1)]
        ;; Remaining sources
        xdr = xd[where(wagn eq 2)]
        rlr = rl[where(wagn eq 2)]
        obscr = obsc[where(wagn eq 2)]
        agnr = wagn[where(wagn eq 2)]
        ;; save survival analysis input
        save,xd,rl,obsc,wagn, $
             xdw,rlw,obscw,agnw, $
             xdr,rlr,obscr,agnr, $
             file='../rsurv_input.sav'
        
        ;; test equal sample size
        ;re = equal_cenc(xd,rl,obsc,wagn,xdeq=xd,rleq=rl,obsceq=obsc,wagneq=wagn)
        ;re = equal_cenc(xdw,rlw,obscw,agnw,xdeq=xdw,rleq=rlw,obsceq=obscw,wagneq=agnw)
        ;re = equal_cenc(xdr,rlr,obscr,agnr,xdeq=xdr,rleq=rlr,obsceq=obscr,wagneq=agnr)
                
        ;; univariate, KM test [0 right-censored, 1 detection, 2 left-censored]
        ;; FACTOR OF (-1) APPLIED WITHIN THE R CODE, NOT HERE
        forprint,xd,rl,obsc,wagn,format='i1,2x,d7.4,2x,i1,2x,i1',textout='rsurv_all.dat',/nocomment
        forprint,xdw,rlw,obscw,agnw,format='i1,2x,d7.4,2x,i1,2x,i1',textout='rsurv_wagn.dat',/nocomment
        forprint,xdr,rlr,obscr,agnr,format='i1,2x,d7.4,2x,i1,2x,i1',textout='rsurv_ragn.dat',/nocomment
        
        print, ''
        print, '*********************************************'
        print, '****    WAITING FOR FURTHER INPUT.       ****'
        print, '****    RUN R SCRIPT...                  ****'
        print, '*********************************************'
        stop
        print, '****    ...READING SURVIVAL OUTPUT       ****'
        print, '****    CONTINUING SURVIVAL ANALYSIS.    ****'
        print, '*********************************************'
        print, ''

        ;; read output -- single
        readcol,'km-output-all.txt',abin,akm,e_akmhi,e_akmlo,format='x,d,d,d,d'
        readcol,'km-output-wagn.txt',wbin,wkm,e_wkmhi,e_wkmlo,format='x,d,d,d,d'
        readcol,'km-output-ragn.txt',rbin,rkm,e_rkmhi,e_rkmlo,format='x,d,d,d,d'

        svbinsz = 0.1
        ad = diff_km(akm,abin,svbinsz,n_elements(xd),binc=abinc)
        wd = diff_km(wkm,wbin,svbinsz,n_elements(xdw),binc=wbinc)
        rd = diff_km(rkm,rbin,svbinsz,n_elements(xdr),binc=rbinc)

        sav_vars = ['ABIN','AKM','E_AKMHI','E_AKMLO','ABIN','ABINC','AD', $
                    'WBIN','WKM','E_WKMHI','E_WKMLO','WBIN','WBINC','WD', $
                    'RBIN','RKM','E_RKMHI','E_RKMLO','RBIN','RBINC','RD', $
                    'SVBINSZ']
                              
        sav_inds = []

        ;; BUILD THE SURVIVAL ANALYSIS SIMULATED RL DISTRIBUTION
        simrla = []
        simrlw = []
        simrlr = []
        for i = 0,n_elements(abinc)-1 do begin
            if (ad[i] gt 0.) then simrla = [simrla,abinc[i]+(randomu(seed,round(ad[i]))-0.5)*svbinsz]
            if (wd[i] gt 0.) then simrlw = [simrlw,wbinc[i]+(randomu(seed,round(wd[i]))-0.5)*svbinsz]
            if (rd[i] gt 0.) then simrlr = [simrlr,rbinc[i]+(randomu(seed,round(rd[i]))-0.5)*svbinsz]
        endfor

        ;; CONVERT SIMULATED RL DIST INTO NH DIST        
        anh_power = rl2nh(-simrla,model='POWER')+(randomu(seed,n_elements(simrla))-0.5)*mean(e_nhdet_power[where(iiqual_det)]) > 21.
        wnh_power = rl2nh(-simrlw,model='POWER')+(randomu(seed,n_elements(simrlw))-0.5)*mean(e_nhdet_power[where(iiqual_det)]) > 21.
        rnh_power = rl2nh(-simrlr,model='POWER')+(randomu(seed,n_elements(simrlr))-0.5)*mean(e_nhdet_power[where(iiqual_det)]) > 21.
        anh_borus = rl2nh(-simrla,model='BORUS')+(randomu(seed,n_elements(simrla))-0.5)*mean(e_nhdet_borus[where(iiqual_det)]) > 21.
        wnh_borus = rl2nh(-simrlw,model='BORUS')+(randomu(seed,n_elements(simrlw))-0.5)*mean(e_nhdet_borus[where(iiqual_det)]) > 21.
        rnh_borus = rl2nh(-simrlr,model='BORUS')+(randomu(seed,n_elements(simrlr))-0.5)*mean(e_nhdet_borus[where(iiqual_det)]) > 21.

        svnhbinsz = scott(wnh_borus);nhbinsz
        yanh_power = histogram(anh_power,locations=xanh_power,bin=svnhbinsz,min=21.-svnhbinsz,max=max(anh_power)+svnhbinsz)
        ywnh_power = histogram(wnh_power,locations=xwnh_power,bin=svnhbinsz,min=21.-svnhbinsz,max=max(wnh_power)+svnhbinsz)
        yrnh_power = histogram(rnh_power,locations=xrnh_power,bin=svnhbinsz,min=21.-svnhbinsz,max=max(rnh_power)+svnhbinsz)
        yanh_borus = histogram(anh_borus,locations=xanh_borus,bin=svnhbinsz,min=21.-svnhbinsz,max=max(anh_borus)+svnhbinsz)
        ywnh_borus = histogram(wnh_borus,locations=xwnh_borus,bin=svnhbinsz,min=21.-svnhbinsz,max=max(wnh_borus)+svnhbinsz)
        yrnh_borus = histogram(rnh_borus,locations=xrnh_borus,bin=svnhbinsz,min=21.-svnhbinsz,max=max(rnh_borus)+svnhbinsz)

        sav_vars = [sav_vars,'SIMRLA','SIMRLW','SIMRLR', $
                             'ANH_POWER','XANH_POWER','YANH_POWER', $
                             'WNH_POWER','XWNH_POWER','YWNH_POWER', $
                             'RNH_POWER','XRNH_POWER','YRNH_POWER', $
                             'ANH_BORUS','XANH_BORUS','YANH_BORUS', $
                             'WNH_BORUS','XWNH_BORUS','YWNH_BORUS', $
                             'RNH_BORUS','XRNH_BORUS','YRNH_BORUS', $
                             'SVNHBINSZ']
        sav_inds = [sav_inds]
        
        ;; read output -- split by obscuration
        readcol,'km-output-all-obsc.txt',abin1,akm1,e_akmhi1,e_akmlo1,format='x,d,d,d,d'
        readcol,'km-output-wagn-obsc.txt',wbin1,wkm1,e_wkmhi1,e_wkmlo1,format='x,d,d,d,d'
        readcol,'km-output-ragn-obsc.txt',rbin1,rkm1,e_rkmhi1,e_rkmlo1,format='x,d,d,d,d'

        i2 = where(diff(abin1) lt 0.)+1
        abin2 = abin1[i2:-1]
        akm2 = akm1[i2:-1]
        e_akmhi2 = e_akmhi1[i2:-1]
        e_akmlo2 = e_akmlo1[i2:-1]
        abin1 = abin1[0:i2-1]
        akm1 = akm1[0:i2-1]
        e_akmhi1 = e_akmhi1[0:i2-1]
        e_akmlo1 = e_akmlo1[0:i2-1]
        
        i2 = where(diff(wbin1) lt 0.)+1
        wbin2 = wbin1[i2:-1]
        wkm2 = wkm1[i2:-1]
        e_wkmhi2 = e_wkmhi1[i2:-1]
        e_wkmlo2 = e_wkmlo1[i2:-1]
        wbin1 = wbin1[0:i2-1]
        wkm1 = wkm1[0:i2-1]
        e_wkmhi1 = e_wkmhi1[0:i2-1]
        e_wkmlo1 = e_wkmlo1[0:i2-1]

        i2 = where(diff(rbin1) lt 0.)+1
        rbin2 = rbin1[i2:-1]
        rkm2 = rkm1[i2:-1]
        e_rkmhi2 = e_rkmhi1[i2:-1]
        e_rkmlo2 = e_rkmlo1[i2:-1]
        rbin1 = rbin1[0:i2-1]
        rkm1 = rkm1[0:i2-1]
        e_rkmhi1 = e_rkmhi1[0:i2-1]
        e_rkmlo1 = e_rkmlo1[0:i2-1]

        ;ad1 = diff_km(akm1,abin1,svbinsz,n_elements(xd),binc=abinc1)
        ;wd1 = diff_km(wkm1,wbin1,svbinsz,n_elements(xdw),binc=wbinc1)
        ;rd1 = diff_km(rkm1,rbin1,svbinsz,n_elements(xdr),binc=rbinc1)
        ;ad2 = diff_km(akm2,abin2,svbinsz,n_elements(xd),binc=abinc2)
        ;wd2 = diff_km(wkm2,wbin2,svbinsz,n_elements(xdw),binc=wbinc2)
        ;rd2 = diff_km(rkm2,rbin2,svbinsz,n_elements(xdr),binc=rbinc2)
        ad1 = diff_km(akm1,abin1,svbinsz,total(obsc eq 1),binc=abinc1)
        wd1 = diff_km(wkm1,wbin1,svbinsz,total(obscw eq 1),binc=wbinc1)
        rd1 = diff_km(rkm1,rbin1,svbinsz,total(obscr eq 1),binc=rbinc1)
        ad2 = diff_km(akm2,abin2,svbinsz,total(obsc eq 2),binc=abinc2)
        wd2 = diff_km(wkm2,wbin2,svbinsz,total(obscw eq 2),binc=wbinc2)
        rd2 = diff_km(rkm2,rbin2,svbinsz,total(obscr eq 2),binc=rbinc2)

        sav_vars = [sav_vars,'ABIN1','AKM1','E_AKMHI1','E_AKMLO1','ABIN1','ABINC1','AD1', $
                             'ABIN2','AKM2','E_AKMHI2','E_AKMLO2','ABIN2','ABINC2','AD2', $
                             'WBIN1','WKM1','E_WKMHI1','E_WKMLO1','WBIN1','WBINC1','WD1', $
                             'WBIN2','WKM2','E_WKMHI2','E_WKMLO2','WBIN2','WBINC2','WD2', $
                             'RBIN1','RKM1','E_RKMHI1','E_RKMLO1','RBIN1','RBINC1','RD1', $
                             'RBIN2','RKM2','E_RKMHI2','E_RKMLO2','RBIN2','RBINC2','RD2']
                              
        sav_inds = []

        ;; BUILD THE SURVIVAL ANALYSIS SIMULATED RL DISTRIBUTION
        simrla1 = []
        simrla2 = []
        simrlw1 = []
        simrlw2 = []
        simrlr1 = []
        simrlr2 = []

        for i = 0,n_elements(wbinc1)-1 do begin
            if (ad1[i] gt 0.) then simrla1 = [simrla1,abinc1[i]+(randomu(seed,round(ad1[i]))-0.5)*svbinsz]
            if (ad2[i] gt 0.) then simrla2 = [simrla2,abinc2[i]+(randomu(seed,round(ad2[i]))-0.5)*svbinsz]
            if (wd1[i] gt 0.) then simrlw1 = [simrlw1,wbinc1[i]+(randomu(seed,round(wd1[i]))-0.5)*svbinsz]
            if (wd2[i] gt 0.) then simrlw2 = [simrlw2,wbinc2[i]+(randomu(seed,round(wd2[i]))-0.5)*svbinsz]
            if (rd1[i] gt 0.) then simrlr1 = [simrlr1,rbinc1[i]+(randomu(seed,round(rd1[i]))-0.5)*svbinsz]
            if (rd2[i] gt 0.) then simrlr2 = [simrlr2,rbinc2[i]+(randomu(seed,round(rd2[i]))-0.5)*svbinsz]
        endfor

        ;; ALL SOURCES
        anh1_power = rl2nh(-simrla1,model='POWER')+(randomu(seed,n_elements(simrla1))-0.5)*mean(e_nhdet_power[where(iiqual_det)]) > 21.
        anh2_power = rl2nh(-simrla2,model='POWER')+(randomu(seed,n_elements(simrla2))-0.5)*mean(e_nhdet_power[where(iiqual_det)]) > 21.
        anh1_borus = rl2nh(-simrla1,model='BORUS')+(randomu(seed,n_elements(simrla1))-0.5)*mean(e_nhdet_borus[where(iiqual_det)]) > 21.
        anh2_borus = rl2nh(-simrla2,model='BORUS')+(randomu(seed,n_elements(simrla2))-0.5)*mean(e_nhdet_borus[where(iiqual_det)]) > 21.

        yanh1_power = histogram(anh1_power,locations=xanh1_power,bin=svnhbinsz)
        yanh2_power = histogram(anh2_power,locations=xanh2_power,bin=svnhbinsz)
        yanh1_borus = histogram(anh1_borus,locations=xanh1_borus,bin=svnhbinsz)
        yanh2_borus = histogram(anh2_borus,locations=xanh2_borus,bin=svnhbinsz)

        xanh1_power = [xanh1_power[0]-svnhbinsz,xanh1_power,xanh1_power[-1]+svnhbinsz]
        xanh2_power = [xanh2_power[0]-svnhbinsz,xanh2_power,xanh2_power[-1]+svnhbinsz]
        yanh1_power = [0.,yanh1_power,0.]
        yanh2_power = [0.,yanh2_power,0.]
        xanh1_borus = [xanh1_borus[0]-svnhbinsz,xanh1_borus,xanh1_borus[-1]+svnhbinsz]
        xanh2_borus = [xanh2_borus[0]-svnhbinsz,xanh2_borus,xanh2_borus[-1]+svnhbinsz]
        yanh1_borus = [0.,yanh1_borus,0.]
        yanh2_borus = [0.,yanh2_borus,0.]

        ;; WISE AGN SOURCES
        wnh1_power = rl2nh(-simrlw1,model='POWER')+(randomu(seed,n_elements(simrlw1))-0.5)*mean(e_nhdet_power[where(iiqual_det)]) > 21.
        wnh2_power = rl2nh(-simrlw2,model='POWER')+(randomu(seed,n_elements(simrlw2))-0.5)*mean(e_nhdet_power[where(iiqual_det)]) > 21.
        wnh1_borus = rl2nh(-simrlw1,model='BORUS')+(randomu(seed,n_elements(simrlw1))-0.5)*mean(e_nhdet_borus[where(iiqual_det)]) > 21.
        wnh2_borus = rl2nh(-simrlw2,model='BORUS')+(randomu(seed,n_elements(simrlw2))-0.5)*mean(e_nhdet_borus[where(iiqual_det)]) > 21.

        ywnh1_power = histogram(wnh1_power,locations=xwnh1_power,bin=svnhbinsz)
        ywnh2_power = histogram(wnh2_power,locations=xwnh2_power,bin=svnhbinsz)
        ywnh1_borus = histogram(wnh1_borus,locations=xwnh1_borus,bin=svnhbinsz)
        ywnh2_borus = histogram(wnh2_borus,locations=xwnh2_borus,bin=svnhbinsz)

        xwnh1_power = [xwnh1_power[0]-svnhbinsz,xwnh1_power,xwnh1_power[-1]+svnhbinsz]
        xwnh2_power = [xwnh2_power[0]-svnhbinsz,xwnh2_power,xwnh2_power[-1]+svnhbinsz]
        ywnh1_power = [0.,ywnh1_power,0.]
        ywnh2_power = [0.,ywnh2_power,0.]
        xwnh1_borus = [xwnh1_borus[0]-svnhbinsz,xwnh1_borus,xwnh1_borus[-1]+svnhbinsz]
        xwnh2_borus = [xwnh2_borus[0]-svnhbinsz,xwnh2_borus,xwnh2_borus[-1]+svnhbinsz]
        ywnh1_borus = [0.,ywnh1_borus,0.]
        ywnh2_borus = [0.,ywnh2_borus,0.]

        ;; REMAINING SOURCES
        rnh1_power = rl2nh(-simrlr1,model='POWER')+(randomu(seed,n_elements(simrlr1))-0.5)*mean(e_nhdet_power[where(iiqual_det)]) > 21.
        rnh2_power = rl2nh(-simrlr2,model='POWER')+(randomu(seed,n_elements(simrlr2))-0.5)*mean(e_nhdet_power[where(iiqual_det)]) > 21.
        rnh1_borus = rl2nh(-simrlr1,model='BORUS')+(randomu(seed,n_elements(simrlr1))-0.5)*mean(e_nhdet_borus[where(iiqual_det)]) > 21.
        rnh2_borus = rl2nh(-simrlr2,model='BORUS')+(randomu(seed,n_elements(simrlr2))-0.5)*mean(e_nhdet_borus[where(iiqual_det)]) > 21.

        yrnh1_power = histogram(rnh1_power,locations=xrnh1_power,bin=svnhbinsz)
        yrnh2_power = histogram(rnh2_power,locations=xrnh2_power,bin=svnhbinsz)
        yrnh1_borus = histogram(rnh1_borus,locations=xrnh1_borus,bin=svnhbinsz)
        yrnh2_borus = histogram(rnh2_borus,locations=xrnh2_borus,bin=svnhbinsz)

        xrnh1_power = [xrnh1_power[0]-svnhbinsz,xrnh1_power,xrnh1_power[-1]+svnhbinsz]
        xrnh2_power = [xrnh2_power[0]-svnhbinsz,xrnh2_power,xrnh2_power[-1]+svnhbinsz]
        yrnh1_power = [0.,yrnh1_power,0.]
        yrnh2_power = [0.,yrnh2_power,0.]
        xrnh1_borus = [xrnh1_borus[0]-svnhbinsz,xrnh1_borus,xrnh1_borus[-1]+svnhbinsz]
        xrnh2_borus = [xrnh2_borus[0]-svnhbinsz,xrnh2_borus,xrnh2_borus[-1]+svnhbinsz]
        yrnh1_borus = [0.,yrnh1_borus,0.]
        yrnh2_borus = [0.,yrnh2_borus,0.]

        sav_vars = [sav_vars,'SIMRLA1','SIMRLA2','SIMRLW1','SIMRLW2','SIMRLR1','SIMRLR2', $
                             'ANH1_POWER','XANH1_POWER','YANH1_POWER', $ 
                             'ANH2_POWER','XANH2_POWER','YANH2_POWER', $ 
                             'WNH1_POWER','XWNH1_POWER','YWNH1_POWER', $ 
                             'WNH2_POWER','XWNH2_POWER','YWNH2_POWER', $ 
                             'RNH1_POWER','XRNH1_POWER','YRNH1_POWER', $ 
                             'RNH2_POWER','XRNH2_POWER','YRNH2_POWER', $ 
                             'ANH1_BORUS','XANH1_BORUS','YANH1_BORUS', $ 
                             'ANH2_BORUS','XANH2_BORUS','YANH2_BORUS', $ 
                             'WNH1_BORUS','XWNH1_BORUS','YWNH1_BORUS', $ 
                             'WNH2_BORUS','XWNH2_BORUS','YWNH2_BORUS', $ 
                             'RNH1_BORUS','XRNH1_BORUS','YRNH1_BORUS', $ 
                             'RNH2_BORUS','XRNH2_BORUS','YRNH2_BORUS'] 
        sav_inds = [sav_inds]
        end
    'IDL': begin
        sfa = rl_surv(rl,xd,0.1)
        sfw = rl_surv(rlw,xdw,0.1)
        sfr = rl_surv(rlr,xdr,0.1)
        end
    else: begin
        print, 'NO VALID SURVIVAL ANALYSIS FORMAT SPECIFIED.'
        PRINT, 'RETURNING.'
        return
        end
endcase

popd
sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="surv_anal.sav"')


END













