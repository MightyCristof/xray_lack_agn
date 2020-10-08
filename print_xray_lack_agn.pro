PRO print_xray_lack_agn, COMMANDS = commands, $
                         RECORD = record, $
                         PROPERTY_CUTS = property_cuts, $
                         FITS = fits, $
                         QUALITY_CUTS = quality_cuts, $
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
;common _clean_cha
;common _clean_xmm
;common _clean_nst
common _quality  
common _combined
common _nhdist

restore,'phot_spec_errs.sav'

file_mkdir,'tables'

if keyword_set(commands) then begin
    print, '    PRINTING NEWCOMMANDS LIST'

    OPENW,99,'tables/my_commands.tex' 
    printf,99,'%% Instruments'
    printf,99,'\newcommand{\nustar}{\emph{NuSTAR}\xspace}'
    printf,99,'\newcommand{\xmm}{\emph{XMM-Newton}\xspace}'
    printf,99,'\newcommand{\chandra}{\emph{Chandra}\xspace}'
    printf,99,'\newcommand{\galex}{\emph{GALEX}\xspace}'
    printf,99,'\newcommand{\wise}{\emph{WISE}\xspace}'
    printf,99,'\newcommand{\xdqso}{\emph{XDQSOz}\xspace}'
    printf,99,'\newcommand{\numaster}{\textsc{numaster}\xspace}'
    printf,99,'\newcommand{\xmmmaster}{\textsc{xmmmaster}\xspace}'
    printf,99,'\newcommand{\chanmaster}{\textsc{chanmaster}\xspace}'
    printf,99,'%% Parameters'
    printf,99,'\newcommand{\ebv}{$E(\bv)_{\text{AGN}}$\xspace}'
    printf,99,'\newcommand{\w}[1]{\emph{W}#1\xspace}'
    printf,99,'\newcommand{\nh}{$N_{\text{H}}$\xspace}'
    printf,99,'\newcommand{\lx}{$L_{\text{X}}\xspace$}'
    printf,99,'\newcommand{\lir}{$L_{\text{MIR}}$\xspace}'
    printf,99,'\newcommand{\lxir}{$L_{\text{X}}(L_{\text{MIR}})$\xspace}'
    printf,99,'\newcommand{\rlum}{$R_{L_{\text{X}}}$\xspace}'
    printf,99,'\newcommand{\llim}{$L_{{\text{X}}_{\text{lim}}}$\xspace}'
    printf,99,'%% Units'
    printf,99,'\newcommand{\kev}{\,{\text{keV}}\xspace}'
    printf,99,'\newcommand{\mum}{\,\mu{\text{m}}\xspace}'
    printf,99,'\newcommand{\cmcm}{\,{\text{cm}}^{-2}\xspace}'
    printf,99,'%%'
    printf,99,'%%'
    printf,99,'%%'
    printf,99,'%% Quality cuts for analysis'    
    printf,99,'\newcommand{\ngal}{'+commas(ngal)+'\xspace}'
    printf,99,'\newcommand{\nagn}{'+commas(nagn)+'\xspace}'
    printf,99,'\newcommand{\ndet}{'+commas(long(total(iidet)))+'\xspace}'
    printf,99,'\newcommand{\nnon}{'+commas(long(total(~iidet)))+'\xspace}'
    printf,99,'\newcommand{\nqgal}{'+commas(0)+'\xspace}'
    printf,99,'\newcommand{\nqagn}{'+commas(long(total(iiqual)))+'\xspace}'
    printf,99,'\newcommand{\nqdet}{'+commas(long(total(iiqual_det)))+'\xspace}'
    printf,99,'\newcommand{\nqnon}{'+commas(long(total(iiqual_non)))+'\xspace}'
    printf,99,'%% WISE AGNs'
    printf,99,'\newcommand{\ninfwac}{'+commas(long(total(iidet_wac)))+'\xspace}'
    printf,99,'\newcommand{\ndetwac}{'+commas(long(total(iidet_wac and iidet)))+'\xspace}'
    printf,99,'\newcommand{\nnonwac}{'+commas(long(total(iidet_wac and ~iidet)))+'\xspace}'
    printf,99,'\newcommand{\nqinfwac}{'+commas(long(total(iidet_wac and iiqual)))+'\xspace}'
    printf,99,'\newcommand{\nqdetwac}{'+commas(long(total(iidet_wac and iiqual_det)))+'\xspace}'
    printf,99,'\newcommand{\nqnonwac}{'+commas(long(total(iidet_wac and iiqual_non)))+'\xspace}'
    printf,99,'%% Chandra'
    printf,99,'\newcommand{\ninfcha}{'+commas(long(total(iiinf_cha)))+'\xspace}'
    printf,99,'\newcommand{\ndetcha}{'+commas(long(total(iiinf_cha and iidet_cha)))+'\xspace}'
    printf,99,'\newcommand{\nnoncha}{'+commas(long(total(iiinf_cha and ~iidet_cha)))+'\xspace}'
    printf,99,'\newcommand{\nqinfcha}{'+commas(long(total(strmatch(xdet,'CHA') or strmatch(xnon,'CHA'))))+'\xspace}'
    printf,99,'\newcommand{\nqdetcha}{'+commas(long(total(strmatch(xdet,'CHA'))))+'\xspace}'
    printf,99,'\newcommand{\nqnoncha}{'+commas(long(total(strmatch(xnon,'CHA'))))+'\xspace}'
    printf,99,'%% XMM'
    printf,99,'\newcommand{\ninfxmm}{'+commas(long(total(iiinf_xmm)))+'\xspace}'
    printf,99,'\newcommand{\ndetxmm}{'+commas(long(total(iiinf_xmm and iidet_xmm)))+'\xspace}'
    printf,99,'\newcommand{\nnonxmm}{'+commas(long(total(iiinf_xmm and ~iidet_xmm)))+'\xspace}'
    printf,99,'\newcommand{\nqinfxmm}{'+commas(long(total(strmatch(xdet,'XMM') or strmatch(xnon,'XMM'))))+'\xspace}'
    printf,99,'\newcommand{\nqdetxmm}{'+commas(long(total(strmatch(xdet,'XMM'))))+'\xspace}'
    printf,99,'\newcommand{\nqnonxmm}{'+commas(long(total(strmatch(xnon,'XMM'))))+'\xspace}'
    printf,99,'%% NuSTAR'
    printf,99,'\newcommand{\ninfnst}{'+commas(long(total(iiinf_nst)))+'\xspace}'
    printf,99,'\newcommand{\ndetnst}{'+commas(long(total(iiinf_nst and iidet_nst)))+'\xspace}'
    printf,99,'\newcommand{\nnonnst}{'+commas(long(total(iiinf_nst and ~iidet_nst)))+'\xspace}'
    printf,99,'\newcommand{\nqinfnst}{'+commas(long(total(strmatch(xdet,'NST') or strmatch(xnon,'NST'))))+'\xspace}'
    printf,99,'\newcommand{\nqdetnst}{'+commas(long(total(strmatch(xdet,'NST'))))+'\xspace}'
    printf,99,'\newcommand{\nqnonnst}{'+commas(long(total(strmatch(xnon,'NST'))))+'\xspace}'
    printf,99,'% Extraneous numbers'
    printf,99,'\newcommand{\nsrc}{'+commas(nsrc)+'\xspace}'
    printf,99,'\newcommand{\ninf}{'+commas(long(total(iiinf)))+'\xspace}'
    printf,99,'\newcommand{\prctwagn}{'+strtrim(round(total(iidet_wac and iiqual)/total(iiqual)*100.),2)+'\%\xspace}'
    printf,99,'\newcommand{\prctukidss}{'+strtrim(round(total(total(bin[where(strmatch(band,'UK?')),*],1) gt 0)/nsrc*100),2)+'\%\xspace}'
    printf,99,'\newcommand{\prcttwom}{'+strtrim(round(total(total(bin[where(strmatch(band,'TWOM?')),*],1) gt 0)/nsrc*100),2)+'\%\xspace}'
    printf,99,'\newcommand{\prctgalex}{'+strtrim(round(total(total(bin[where(strmatch(band,'GALEX?')),*],1) gt 0)/nsrc*100),2)+'\%\xspace}'
    printf,99,'\newcommand{\prctspec}{${\sim}$'+strtrim(round(total(strmatch(ztype,'ZS*'))/nsrc*100),2)+'\%\xspace}'
    printf,99,'\newcommand{\prctphot}{${\sim}$'+strtrim(round(total(~strmatch(ztype,'ZS*'))/nsrc*100),2)+'\%\xspace}'
    printf,99,'\newcommand{\prctcorr}{${\sim}$'+strtrim(round(total(iicorr ne 0)/nsrc*100),2)+'\%\xspace}'
    printf,99,'\newcommand{\prctdetwagn}{'+strtrim(round(total(iidet_wac and iiqual_det)/total(iiqual)*100),2)+'\%\xspace}'
    printf,99,'\newcommand{\prctnonwagn}{'+strtrim(round(total(iidet_wac and iiqual_non)/total(iiqual)*100),2)+'\%\xspace}'
    printf,99,'\newcommand{\prctctwagn}{${\sim}$'+strtrim(round(total(iidet_wac and iiqual_non and llnon le (rl2nh(alog10(1.5e24),model='BORUS',/lum_out))[0])/total(iidet_wac and iiqual_non)*100),2)+'\%\xspace}'
    printf,99,'\newcommand{\nexages}{'+commas(nmages)+'\xspace}'
    printf,99,'\newcommand{\ninages}{'+commas(nmqages)+'\xspace}'
    printf,99,'\newcommand{\nexsoff}{'+commas(nmagesoff)+'\xspace}'
    printf,99,'\newcommand{\zexsig}{'+string(rnd(zsig,3),format="(d5.3)")+'\xspace}'
    printf,99,'\newcommand{\zinsig}{'+string(rnd(zqsig,3),format="(d5.3)")+'\xspace}'
    CLOSE,99
endif


if keyword_set(record) then begin
    print, '    PRINTING RECORDS TABLE'

    OPENW,1,'tables/table1.tex'
    printf,1,'\begin{deluxetable*}{lrrrrrrr}'
    printf,1,'\tablecaption{Number of sources with various telescope coverage.\label{sample}}'
    printf,1,'\tablehead{\colhead{} & \colhead{SDSS} & \colhead{XDQSOz} & \colhead{\wise} & \colhead{unWISE} & \colhead{UKIDSS} & \colhead{2MASS} & \colhead{\galex}}'
    printf,1,'\startdata'
    printf,1,'Initial & \nisdss & \nixdqso & \niwise & \niunwise & \niukidss & \nitwom & \nigalex  \\'
    printf,1,'Final   & \nsdss  & \nxdqso  & \nwise  & \nunwise  & \nukidss  & \ntwom  & \ngalex  \\'
    printf,1,'\hline'
    printf,1,' & \multicolumn{2}{r}{\chandra} & \multicolumn{2}{r}{\xmm} & \multicolumn{2}{r}{\nustar}  \\'
    printf,1,'\hline'
    printf,1,'Initial & \multicolumn{2}{r}{\nicha}   & \multicolumn{2}{r}{\nixmm}   & \multicolumn{2}{r}{\ninst}  \\'
    printf,1,'Final   & \multicolumn{2}{r}{\ninfcha} & \multicolumn{2}{r}{\ninfxmm} & \multicolumn{2}{r}{\ninfnst}  \\'
    printf,1,'\enddata'
    printf,1,'%\tablecomments{}'
    printf,1,'\end{deluxetable*}'
    CLOSE,1
endif



if keyword_set(fits) then begin
    print, '    PRINTING FIT OUTPUT TABLE'
    rows = 10
    ;; shorthand 
    ind = where(iiqual,ct)
    if (ct eq 0.) then stop
    inst = [['-','c'],['-','x'],['-','n']]

    ;; LX field
    lx_str = strarr(nsrc)
    lx_str[where(iiqual_non)] = '\textless \,'+string(loglxlim[where(iiqual_non)],format='(f5.2)')
    lx_str[where(iiqual_det)] = ' '+string(loglx[where(iiqual_det)],format='(f5.2)')+'$\pm$'+string(e_loglx[where(iiqual_det)],format='(f5.2)')+' '
    
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
    
    OPENW,3,'tables/table3.tex' 
    printf,3,'\begin{deluxetable*}{rrrrrrrrrr}'
    printf,3,'\tablecolumns{12}'
    printf,3,'\tablecaption{Source data and SED modeling output parameters.\label{sedfits}}'
    printf,3,'\tablehead{\colhead{RA} & \colhead{Dec} & \colhead{$z$} & \colhead{$E(\bv)_{\text{AGN}}$} & \colhead{$\chi^2 /$ DoF} & \colhead{$P_{\text{AGN}}$} & \colhead{$\log L_{\text{IR}}$} & \colhead{$\log L_{\text{X}}$} & \colhead{\wise} & \colhead{X-ray} \\' 
    printf,3,'\colhead{[deg]} & \colhead{[deg]} & \colhead{} & \colhead{} & \colhead{} & \colhead{[\%]} & \colhead{[erg\,s$^{-1}$]} & \colhead{[erg\,s$^{-1}$]} & \colhead{AGN} & \colhead{Field}}'
    printf,3,'\startdata'
    for i = 0,rows-1 do printf,3,line[i]
    printf,3,'\enddata'
    printf,3,'\tablecomments{Source data and modeling output parameters: right ascension, declination, redshift, color excess \ebv, chi-squared/degrees of freedom, percentage of resampled realizations with AGN contribution ($C_{\text{AGN}} > 0$), AGN $6\mum$ luminosity, X-ray $2{\text{--}}10\kev$ luminosity, boolean flag if matched to the \wise AGN Catalog (0=false, 1=true), and observations by X-ray field (\chandra, \xmm, \nustar). Uncertainties on \ebv, $\log L_{\text{IR}}$, and $\log L_{\text{X}}$ were estimated via the median absolute deviation of all source realizations containing AGN contribution.}'
    printf,3,'%\tablerefs{}'
    printf,3,'\end{deluxetable*}'
    CLOSE,3
endif


if keyword_set(mrt) then begin
    print, '    PRINTING FIT OUTPUT TABLE'
    ;rows = 10
    ;; shorthand 
    ind = where(iiqual,ct)
    if (ct eq 0.) then stop
    inst = [['-','c'],['-','x'],['-','n']]

    ;; LX field
    lt_str = strarr(nsrc)
    lt_str[where(iiqual_non)] = '<'
    lx_str = strarr(nsrc)
    lx_str[where(iiqual_non)] = string(lxlim[where(iiqual_non)]*1e-7,format='(e10.4)')
    lx_str[where(iiqual_det)] = string(lx[where(iiqual_det)]*1e-7,format='(e10.4)')
    lxerr_str = strarr(nsrc)
    lxerr_str[where(iiqual_det)] = string(e_lx[where(iiqual_det)]*1e-7,format='(e10.4)')
    
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



if keyword_set(quality_cuts) then begin
    print, '    PRINTING QUALITY CUTS TABLE'

    OPENW,4,'tables/table4.tex'
    printf,4,'\begin{deluxetable}{lrr}'
    printf,4,'\tablecolumns{3}'
    printf,4,'\tablecaption{Number of sources with various selections.\label{sample}}'
    printf,4,'\tablehead{\colhead{Selection} & \colhead{Final sample} & \colhead{Analysis subset}\\ '
    printf,4,'\colhead{} & \colhead{\#} & \colhead{\#}} '
    printf,4,'\startdata'
    printf,4,'SED Galaxy                  &  \ngal     &  \nqgal  \\'
    printf,4,'SED Galaxy+AGN              &  \nagn     &  \nqagn  \\'
    printf,4,'X-ray detected              &  \ndet     &  \nqdet  \\'
    printf,4,'X-ray non-det.              &  \nnon     &  \nqnon  \\'
    printf,4,'\cutinhead{\wise AGNs}'
    printf,4,'in \emph{WISE} AGN Catalog  &  \ninfwac  &  \nqinfwac  \\'
    printf,4,'X-ray detected              &  \ndetwac  &  \nqdetwac  \\'
    printf,4,'X-ray non-det.              &  \nnonwac  &  \nqnonwac  \\'
    printf,4,'\cutinhead{\chandra}'
    printf,4,'in X-ray field              &  \ninfcha  &  \nqinfcha  \\'
    printf,4,'X-ray detected              &  \ndetcha  &  \nqdetcha  \\'
    printf,4,'X-ray non-det.              &  \nnoncha  &  \nqnoncha  \\'
    printf,4,'\cutinhead{\xmm}'
    printf,4,'in X-ray field              &  \ninfxmm  &  \nqinfxmm  \\'
    printf,4,'X-ray detected              &  \ndetxmm  &  \nqdetxmm  \\'
    printf,4,'X-ray non-det.              &  \nnonxmm  &  \nqnonxmm  \\'
    printf,4,'\cutinhead{\nustar}'
    printf,4,'in X-ray field              &  \ninfnst  &  \nqinfnst  \\'
    printf,4,'X-ray detected              &  \ndetnst  &  \nqdetnst  \\'
    printf,4,'X-ray non-det.              &  \nnonnst  &  \nqnonnst  \\'
    printf,4,'\enddata'
    printf,4,'\tablecomments{Total initial sample: \nsrc sources, which exist within \chandra, \xmm, and/or \nustar observations.}'
    printf,4,'%\tablerefs{}'
    printf,4,'\end{deluxetable}'
    CLOSE,4
endif






END


