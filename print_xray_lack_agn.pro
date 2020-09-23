PRO print_xray_lack_agn, COMMANDS = commands, $
                         RECORD = record, $
                         FITS = fits, $
                         MRT = mrt
                         
                         
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
common _clean_cha
common _clean_xmm
common _clean_nst
common _quality  
common _combined
common _nhdist

restore,'phot_spec_errs.sav'

file_mkdir,'tables'

if keyword_set(commands) then begin
    print, '    PRINTING NEWCOMMANDS LIST'
    table_name='tables/list_newcommands.tex' 
    OPENW,1,table_name
    printf,1,'%% Sample numbers'
    printf,1,'\newcommand{\rawsdss}{187,176,677\xspace}'
    printf,1,'\newcommand{\rawxdqso}{5,537,436\xspace}'
    printf,1,'\newcommand{\rawwise}{747,634,026\xspace}'
    printf,1,'\newcommand{\rawunw}{469,037,483\xspace}'
    printf,1,'\newcommand{\rawuk}{84,440,794\xspace}'
    printf,1,'\newcommand{\raw2m}{470,992,970\xspace}'
    printf,1,'\newcommand{\rawgalex}{\xspace}'

    printf,1,'\newcommand{\ninit}{80,236,261\xspace}'
    printf,1,'\newcommand{\nisdss}{78,794,318\xspace}'
    printf,1,'\newcommand{\nizsupp}{3,672\xspace}'
    printf,1,'\newcommand{\nixdqso}{1,441,943\xspace}'
    printf,1,'\newcommand{\niwise}{80,236,261\xspace}'
    printf,1,'\newcommand{\niunw}{80,230,855\xspace}'
    printf,1,'\newcommand{\niuk}{18,671,693\xspace}'
    printf,1,'\newcommand{\ni2m}{4,271,320\xspace}'
    printf,1,'\newcommand{\nigalex}{7,775,179\xspace}'

    printf,1,'\newcommand{\nfin}{80,236,261\xspace}'
    printf,1,'\newcommand{\nfsdss}{78,794,318\xspace}'
    printf,1,'\newcommand{\nfzsupp}{3,672\xspace}'
    printf,1,'\newcommand{\nfxdqso}{1,441,943\xspace}'
    printf,1,'\newcommand{\nfwise}{80,236,261\xspace}'
    printf,1,'\newcommand{\nfunw}{80,230,855\xspace}'
    printf,1,'\newcommand{\nfuk}{18,671,693\xspace}'
    printf,1,'\newcommand{\nf2m}{0\xspace}'
    printf,1,'\newcommand{\nfgalex}{7,775,179\xspace}'

    
    printf,1,'\newcommand{\ngal}{'+commas(ngal)+'\xspace}'
    printf,1,'\newcommand{\nagn}{'+commas(nagn)+'\xspace}'
    printf,1,'\newcommand{\ndet}{'+commas(n_elements(where(iidet,/null)))+'\xspace}'
    printf,1,'\newcommand{\nnon}{'+commas(n_elements(where(~iidet,/null)))+'\xspace}'
    printf,1,'\newcommand{\nfgal}{'+commas(0)+'\xspace}'
    printf,1,'\newcommand{\nfagn}{'+commas(n_elements(where(iifinal,/null)))+'\xspace}'
    printf,1,'\newcommand{\nfdet}{'+commas(n_elements(where(iifinal_det,/null)))+'\xspace}'
    printf,1,'\newcommand{\nfnon}{'+commas(n_elements(where(iifinal_non,/null)))+'\xspace}'
    printf,1,'%% WISE AGNs'
    printf,1,'\newcommand{\ninfwac}{'+commas(n_elements(where(iidet_wac,/null)))+'\xspace}'
    printf,1,'\newcommand{\ndetwac}{'+commas(n_elements(where(iidet_wac and iidet,/null)))+'\xspace}'
    printf,1,'\newcommand{\nnonwac}{'+commas(n_elements(where(iidet_wac and ~iidet,/null)))+'\xspace}'
    printf,1,'\newcommand{\nfinfwac}{'+commas(n_elements(where(iidet_wac and iifinal,/null)))+'\xspace}'
    printf,1,'\newcommand{\nfdetwac}{'+commas(n_elements(where(iidet_wac and iifinal_det,/null)))+'\xspace}'
    printf,1,'\newcommand{\nfnonwac}{'+commas(n_elements(where(iidet_wac and iifinal_non,/null)))+'\xspace}'
    printf,1,'%% Chandra'
    printf,1,'\newcommand{\ninfcha}{'+commas(n_elements(where(iiinf_cha,/null)))+'\xspace}'
    printf,1,'\newcommand{\ndetcha}{'+commas(n_elements(where(iiinf_cha and iidet_cha,/null)))+'\xspace}'
    printf,1,'\newcommand{\nnoncha}{'+commas(n_elements(where(iiinf_cha and ~iidet_cha,/null)))+'\xspace}'
    printf,1,'\newcommand{\nfinfcha}{'+commas(n_elements(where(iiinf_cha and (iifinal_det_cha or iifinal_non_cha),/null)))+'\xspace}'
    printf,1,'\newcommand{\nfdetcha}{'+commas(n_elements(where(iifinal_det_cha,/null)))+'\xspace}'
    printf,1,'\newcommand{\nfnoncha}{'+commas(n_elements(where(iifinal_non_cha,/null)))+'\xspace}'
    printf,1,'%% XMM'
    printf,1,'\newcommand{\ninfxmm}{'+commas(n_elements(where(iiinf_xmm,/null)))+'\xspace}'
    printf,1,'\newcommand{\ndetxmm}{'+commas(n_elements(where(iiinf_xmm and iidet_xmm,/null)))+'\xspace}'
    printf,1,'\newcommand{\nnonxmm}{'+commas(n_elements(where(iiinf_xmm and ~iidet_xmm,/null)))+'\xspace}'
    printf,1,'\newcommand{\nfinfxmm}{'+commas(n_elements(where(iiinf_xmm and (iifinal_det_xmm or iifinal_non_xmm),/null)))+'\xspace}'
    printf,1,'\newcommand{\nfdetxmm}{'+commas(n_elements(where(iifinal_det_xmm,/null)))+'\xspace}'
    printf,1,'\newcommand{\nfnonxmm}{'+commas(n_elements(where(iifinal_non_xmm,/null)))+'\xspace}'
    printf,1,'%% NuSTAR'
    printf,1,'\newcommand{\ninfnst}{'+commas(n_elements(where(iiinf_nst,/null)))+'\xspace}'
    printf,1,'\newcommand{\ndetnst}{'+commas(n_elements(where(iiinf_nst and iidet_nst,/null)))+'\xspace}'
    printf,1,'\newcommand{\nnonnst}{'+commas(n_elements(where(iiinf_nst and ~iidet_nst,/null)))+'\xspace}'
    printf,1,'\newcommand{\nfinfnst}{'+commas(n_elements(where(iiinf_nst and (iifinal_det_nst or iifinal_non_nst),/null)))+'\xspace}'
    printf,1,'\newcommand{\nfdetnst}{'+commas(n_elements(where(iifinal_det_nst,/null)))+'\xspace}'
    printf,1,'\newcommand{\nfnonnst}{'+commas(n_elements(where(iifinal_non_nst,/null)))+'\xspace}'
    printf,1,'% Extraneous numbers'
    printf,1,'\newcommand{\nsrc}{'+commas(nsrc)+'\xspace}'
    printf,1,'\newcommand{\ninf}{'+commas(n_elements(where(iiinf,/null)))+'\xspace}'
    printf,1,'\newcommand{\prctwagn}{'+strtrim(ceil(total(iidet_wac and iifinal)/total(iifinal)*100.),2)+'\%\xspace}'
    printf,1,'\newcommand{\prctukidss}{'+strtrim(ceil(total(total(bin[9:12,*],1) gt 0)/nsrc*100),2)+'\%\xspace}'
    printf,1,'\newcommand{\prcttwom}{'+strtrim(ceil(total(total(bin[13:15,*],1) gt 0)/nsrc*100),2)+'\%\xspace}'
    printf,1,'\newcommand{\prctgalex}{'+strtrim(ceil(total(total(bin[16:17,*],1) gt 0)/nsrc*100),2)+'\%\xspace}'
    printf,1,'\newcommand{\prctspec}{${\sim}$'+strtrim(round(total(strmatch(ztype,'ZS*'))/nsrc*100),2)+'\%\xspace}'
    printf,1,'\newcommand{\prctphot}{${\sim}$'+strtrim(round(total(strmatch(ztype,'ZP'))/nsrc*100),2)+'\%\xspace}'
    printf,1,'\newcommand{\prctcorr}{${\sim}$'+strtrim(ceil(total(iicorr and iifinal)/total(iifinal)*100),2)+'\%\xspace}'
    printf,1,'\newcommand{\prctdetwagn}{'+strtrim(round(total(iidet_wac and iifinal_det)/total(iifinal)*100),2)+'\%\xspace}'
    printf,1,'\newcommand{\prctnonwagn}{'+strtrim(round(total(iidet_wac and iifinal_non)/total(iifinal)*100),2)+'\%\xspace}'
    printf,1,'\newcommand{\prctctwagn}{${\sim}$'+strtrim(round(total(iidet_wac and iifinal_non and llnon le (rl2nh(alog10(1.5e24),model='BORUS',/lum_out))[0])/total(iidet_wac and iifinal_non)*100),2)+'\%\xspace}'
    printf,1,'\newcommand{\mad}{'+string(rnd(mad,3),format="(d5.3)")+'\xspace}'
    printf,1,'\newcommand{\madfin}{'+string(rnd(madfin,3),format="(d5.3)")+'\xspace}'
    CLOSE,1
endif

if keyword_set(record) then begin
    print, '    PRINTING RECORDS TABLE'
    table_name='tables/table1.tex'
    OPENW,2,table_name
    printf,2,'\begin{deluxetable}{lrr}'
    printf,2,'\tablecolumns{3}'
    printf,2,'\tablecaption{Number of sources with various selections.\label{sample}}'
    printf,2,'\tablehead{\colhead{Selection} & \colhead{Initial sample} & \colhead{Final sample}\\ '
    printf,2,'\colhead{} & \colhead{\#} & \colhead{\#}} '
    printf,2,'\startdata'
    printf,2,'SED Galaxy                  &  \ngal     &  \nfgal  \\'
    printf,2,'SED Galaxy+AGN              &  \nagn     &  \nfagn  \\'
    printf,2,'X-ray detected              &  \ndet     &  \nfdet  \\'
    printf,2,'X-ray non-det.              &  \nnon     &  \nfnon  \\'
    printf,2,'\cutinhead{\wise AGNs}'
    printf,2,'in \emph{WISE} AGN Catalog  &  \ninfwac  &  \nfinfwac  \\'
    printf,2,'X-ray detected              &  \ndetwac  &  \nfdetwac  \\'
    printf,2,'X-ray non-det.              &  \nnonwac  &  \nfnonwac  \\'
    printf,2,'\cutinhead{\chandra}'
    printf,2,'in X-ray field              &  \ninfcha  &  \nfinfcha  \\'
    printf,2,'X-ray detected              &  \ndetcha  &  \nfdetcha  \\'
    printf,2,'X-ray non-det.              &  \nnoncha  &  \nfnoncha  \\'
    printf,2,'\cutinhead{\xmm}'
    printf,2,'in X-ray field              &  \ninfxmm  &  \nfinfxmm  \\'
    printf,2,'X-ray detected              &  \ndetxmm  &  \nfdetxmm  \\'
    printf,2,'X-ray non-det.              &  \nnonxmm  &  \nfnonxmm  \\'
    printf,2,'\cutinhead{\nustar}'
    printf,2,'in X-ray field              &  \ninfnst  &  \nfinfnst  \\'
    printf,2,'X-ray detected              &  \ndetnst  &  \nfdetnst  \\'
    printf,2,'X-ray non-det.              &  \nnonnst  &  \nfnonnst  \\'
    printf,2,'\enddata'
    printf,2,'\tablecomments{Total initial sample: \nsrc sources, which exist within \chandra, \xmm, and/or \nustar observations.}'
    printf,2,'%\tablerefs{}'
    printf,2,'\end{deluxetable}'
    CLOSE,2
endif



if keyword_set(fits) then begin
    print, '    PRINTING FIT OUTPUT TABLE'
    rows = 10
    ;; shorthand 
    ind = where(iifinal,ct)
    if (ct eq 0.) then stop
    inst = [['-','c'],['-','x'],['-','n']]

    ;; LX field
    lx_str = strarr(nsrc)
    lx_str[where(iifinal_non)] = '\textless \,'+string(loglxlim[where(iifinal_non)],format='(f5.2)')
    lx_str[where(iifinal_det)] = ' '+string(loglx[where(iifinal_det)],format='(f5.2)')+'$\pm$'+string(e_loglx[where(iifinal_det)],format='(f5.2)')+' '
    
    ;; table output strings
    id_str = strtrim(objid[ind],2)
    ra_str = strtrim(ra[ind],2)
    dec_str = strtrim(dec[ind],2)
    z_str = string(rnd(z[ind],3),format='(f5.3)')+'$\pm$'+strtrim(string(rnd(zerr[ind],3),format='(f5.3)'),2)
    ;ebv_str = string(rnd(ebv[ind],2),format='(f5.2)')+'\pm'+string(rnd(ebv_sigm[1,ind],2),format='(f5.2)')
    ebv_str = string(rnd(ebv[ind]<50.0,2),format='(f5.2)')+'$\pm$'+strtrim(string(rnd(e_ebv[ind],2),format='(f5.2)'),2)
    chi_str = string(rnd(chi[ind],2),format='(f6.2)')+' / '+strtrim(fix(dof[ind]),2)
    perc_str = string(perc_agn[ind],format='(i3)')
    ;lir_str = string(rnd(loglir[ind],2),format='(f5.2)')+'\pm'+string(rnd(lir_sigm[1,ind]/(alog(10.)*lir[ind]),2),format='(f5.2)')
    lir_str = string(rnd(loglir[ind],2),format='(f5.2)')+'$\pm$'+strtrim(string(rnd(e_loglir[ind],2),format='(f5.2)'),2)
    ;lx_str = string(rnd(loglx[ind],2),format='(f5.2)')+'\pm'+string(rnd(e_loglx[ind],2),format='(f5.2)')
    lx_str = lx_str[ind]
    wagn_str = strtrim(fix(iidet_wac[ind]),2)
    field_str = inst[iiinf_cha[ind],0]+'/'+inst[iiinf_xmm[ind],1]+'/'+inst[iiinf_nst[ind],2]
    
    ;; all output lines
    ;line = id_str+'  &  '+ra_str+'  &  '+dec_str+'  &  '+z_str+'  &  '+ebv_str+'  &  '+chi_str+'  &  '+perc_str+'  &  '+lir_str+'  &  '+lx_str+'  &  '+wagn_str+'  &  '+field_str+'  \\ '
    line = ra_str+'  &  '+dec_str+'  &  '+z_str+'  &  '+ebv_str+'  &  '+chi_str+'  &  '+perc_str+'  &  '+lir_str+'  &  '+lx_str+'  &  '+wagn_str+'  &  '+field_str+'  \\ '
    isort = sort(ra[ind])
    line = line[reverse(isort)]

    table_name='tables/table2.tex' 
    OPENW,3,table_name
    ;printf,3,'\begin{deluxetable*}{ccccccccccc}'
    printf,3,'\begin{deluxetable*}{rrrrrrrrrr}'
    printf,3,'\tablecolumns{12}'
    printf,3,'\tablecaption{Source and SED modeling output parameters.\label{fits}}'
    ;printf,3,'\tablehead{\colhead{SDSS ObjID} & \colhead{RA} & \colhead{Dec} & \colhead{z} & \colhead{$E(\bv)$} & \colhead{$\chi^2 /$ df.} & \colhead{$C_{\text{AGN}}$} & \colhead{$\log L_{\text{IR}}$} & \colhead{$\log L_{\text{X}}$} & \colhead{\wise} & \colhead{X-ray}\\ '
    printf,3,'\tablehead{\colhead{RA} & \colhead{Dec} & \colhead{z} & \colhead{$E(\bv)$} & \colhead{$\chi^2 /$ df.} & \colhead{$C_{\text{AGN}}$} & \colhead{$\log L_{\text{IR}}$} & \colhead{$\log L_{\text{X}}$} & \colhead{\wise} & \colhead{X-ray}\\ '
    ;printf,3,'\colhead{} & \colhead{[deg]} & \colhead{[deg]} & \colhead{} & \colhead{} & \colhead{} & \colhead{[\%]} & \colhead{[erg\,s$^{-1}$]} & \colhead{[erg\,s$^{-1}$]} & \colhead{AGN} & \colhead{Field}} '
    printf,3,'\colhead{[deg]} & \colhead{[deg]} & \colhead{} & \colhead{} & \colhead{} & \colhead{[\%]} & \colhead{[erg\,s$^{-1}$]} & \colhead{[erg\,s$^{-1}$]} & \colhead{AGN} & \colhead{Field}} '
    printf,3,'\startdata'
    for i = 0,rows-1 do printf,3,line[i]
    printf,3,'\enddata'
    printf,3,'\tablecomments{Source data and modeling output parameters: right ascension, declination, redshift, color excess \ebv, chi-squared--degrees of freedom, percent of realizations with AGN contribution, IR luminosity, X-ray luminosity, matched to the \wise AGN Catalog, and observations by X-ray field (\chandra, \xmm, \nustar). Uncertainties on \ebv, $\log L_{\text{IR}}$, and $\log L_{\text{X}}$ were estimated via the median aboolute deviation of all source realizations containing AGN contribution.}'
    printf,3,'%\tablerefs{}'
    printf,3,'\end{deluxetable*}'
    CLOSE,3
endif


if keyword_set(mrt) then begin
    print, '    PRINTING FIT OUTPUT TABLE'
    ;rows = 10
    ;; shorthand 
    ind = where(iifinal,ct)
    if (ct eq 0.) then stop
    inst = [['-','c'],['-','x'],['-','n']]

    ;; LX field
    lt_str = strarr(nsrc)
    lt_str[where(iifinal_non)] = '<'
    lx_str = strarr(nsrc)
    lx_str[where(iifinal_non)] = string(lxlim[where(iifinal_non)]*1e-7,format='(e10.4)')
    lx_str[where(iifinal_det)] = string(lx[where(iifinal_det)]*1e-7,format='(e10.4)')
    lxerr_str = strarr(nsrc)
    lxerr_str[where(iifinal_det)] = string(e_lx[where(iifinal_det)]*1e-7,format='(e10.4)')
    
    ;; table output strings
    id_str = strtrim(objid[ind],2)
    ra_str = strtrim(ra[ind],2)
    dec_str = strtrim(dec[ind],2)
    z_str = string(rnd(z[ind],3),format='(f5.3)')
    zerr_str = strtrim(string(rnd(zerr[ind],3),format='(f5.3)'),2)
    ;ebv_str = string(rnd(ebv[ind],2),format='(f5.2)')+'\pm'+string(rnd(ebv_sigm[1,ind],2),format='(f5.2)')
    ebv_str = string(rnd(ebv[ind]<50.0,2),format='(f5.2)')
    ebverr_str = strtrim(string(rnd(e_ebv[ind],2),format='(f5.2)'),2)
    chi_str = string(rnd(chi[ind],2),format='(f6.2)')
    dof_str = strtrim(fix(dof[ind]),2)
    perc_str = string(agn_perc[ind],format='(i3)')
    ;lir_str = string(rnd(loglir[ind],2),format='(f5.2)')+'\pm'+string(rnd(lir_sigm[1,ind]/(alog(10.)*lir[ind]),2),format='(f5.2)')
    lir_str = string(lir[ind]*1e-7,format='(e10.4)')
    lirerr_str = string(e_lir[ind]*1e-7,format='(e10.4)')
    ;lx_str = string(rnd(loglx[ind],2),format='(f5.2)')+'\pm'+string(rnd(e_loglx[ind],2),format='(f5.2)')
    lt_str = lt_str[ind]
    lx_str = lx_str[ind]
    lxerr_str = lxerr_str[ind]
    wagn_str = strtrim(fix(iidet_wac[ind]),2)
    field_str = inst[iiinf_cha[ind],0]+'/'+inst[iiinf_xmm[ind],1]+'/'+inst[iiinf_nst[ind],2]
    
    ;; all output lines
    ;line = id_str+'  &  '+ra_str+'  &  '+dec_str+'  &  '+z_str+'  &  '+ebv_str+'  &  '+chi_str+'  &  '+perc_str+'  &  '+lir_str+'  &  '+lx_str+'  &  '+wagn_str+'  &  '+field_str+'  \\ '
    line = ra_str+'  &  '+dec_str+'  &  '+z_str+'  &  '+zerr_str+'  &  '+ebv_str+'  &  '+ebverr_str+'  &  '+chi_str+'  &  '+dof_str+'  &  '+perc_str+'  &  '+lir_str+'  &  '+lirerr_str+'  &  '+lt_str+'  &  '+lx_str+'  &  '+lxerr_str+'  &  '+wagn_str+'  &  '+field_str
    isort = sort(ra[ind])
    line = line[reverse(isort)]

    rows = n_elements(line)
    
    table_name='tables/table2.txt' 
    OPENW,4,table_name
    ;printf,3,'\begin{deluxetable*}{ccccccccccc}'
    ;printf,4,'\begin{deluxetable*}{rrrrrrrrrr}'
    ;printf,4,'\tablecolumns{12}'
    ;printf,4,'\tablecaption{Source and SED modeling output parameters.\label{fits}}'
    ;printf,3,'\tablehead{\colhead{SDSS ObjID} & \colhead{RA} & \colhead{Dec} & \colhead{z} & \colhead{$E(\bv)$} & \colhead{$\chi^2 /$ df.} & \colhead{$C_{\text{AGN}}$} & \colhead{$\log L_{\text{IR}}$} & \colhead{$\log L_{\text{X}}$} & \colhead{\wise} & \colhead{X-ray}\\ '
    ;printf,4,'\tablehead{\colhead{RA} & \colhead{Dec} & \colhead{z} & \colhead{$E(\bv)$} & \colhead{$\chi^2 /$ df.} & \colhead{$C_{\text{AGN}}$} & \colhead{$\log L_{\text{IR}}$} & \colhead{$\log L_{\text{X}}$} & \colhead{\wise} & \colhead{X-ray}\\ '
    ;printf,3,'\colhead{} & \colhead{[deg]} & \colhead{[deg]} & \colhead{} & \colhead{} & \colhead{} & \colhead{[\%]} & \colhead{[erg\,s$^{-1}$]} & \colhead{[erg\,s$^{-1}$]} & \colhead{AGN} & \colhead{Field}} '
    ;printf,4,'\colhead{[deg]} & \colhead{[deg]} & \colhead{} & \colhead{} & \colhead{} & \colhead{[\%]} & \colhead{[erg\,s$^{-1}$]} & \colhead{[erg\,s$^{-1}$]} & \colhead{AGN} & \colhead{Field}} '
    ;printf,4,'\startdata'
    for i = 0,rows-1 do printf,4,line[i]
    ;printf,4,'\enddata'
    ;printf,4,'\tablecomments{Source data and modeling output parameters: right ascension, declination, redshift, color excess \ebv, chi-squared--degrees of freedom, percent of realizations with AGN contribution, IR luminosity, X-ray luminosity, matched to the \wise AGN Catalog, and observations by X-ray field (\chandra, \xmm, \nustar). Uncertainties on \ebv, $\log L_{\text{IR}}$, and $\log L_{\text{X}}$ were estimated via the median aboolute deviation of all source realizations containing AGN contribution.}'
    ;printf,4,'%\tablerefs{}'
    ;printf,4,'\end{deluxetable*}'
    CLOSE,4
endif



END


