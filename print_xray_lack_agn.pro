PRO print_xray_lack_agn, COMMANDS = commands, $
                         RECORD = record, $
                         PROPERTY_CUTS = property_cuts, $
                         XRAY = xray, $
                         MRT_XRAY = mrt_xray, $
                         SED = sed, $
                         MRT_SED = mrt_sed, $
                         QUALITY_CUTS = quality_cuts, $
                         XSTACK = xstack
                         
                         
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

restore,'phot_spec_errs.sav'

file_mkdir,'tables'

if keyword_set(commands) then begin
    print, '    PRINTING NEWCOMMANDS LIST'

    openw,1,'tables/my_commands.tex' 
    printf,1,'%% Instruments'
    printf,1,'\newcommand{\chandra}{\emph{Chandra}\xspace}'
    printf,1,'\newcommand{\xmmn}{\emph{XMM-Newton}\xspace}'
    printf,1,'\newcommand{\xmm}{\emph{XMM}\xspace}'
    printf,1,'\newcommand{\nustar}{\emph{NuSTAR}\xspace}'
    printf,1,'\newcommand{\galex}{\emph{GALEX}\xspace}'
    printf,1,'\newcommand{\wise}{\emph{WISE}\xspace}'
    printf,1,'\newcommand{\xdqso}{\emph{XDQSOz}\xspace}'
    printf,1,'\newcommand{\chanmaster}{\textsc{chanmaster}\xspace}'
    printf,1,'\newcommand{\xmmmaster}{\textsc{xmmmaster}\xspace}'
    printf,1,'\newcommand{\numaster}{\textsc{numaster}\xspace}'
    printf,1,'\newcommand{\borus}{\textsc{borus}\xspace}'
    printf,1,'%% Parameters'
    printf,1,'\newcommand{\ebv}{$E(\bv)_{\text{AGN}}$\xspace}'
    printf,1,'\newcommand{\w}[1]{\emph{W}#1\xspace}'
    printf,1,'\newcommand{\nh}{$N_{\text{H}}$\xspace}'
    printf,1,'\newcommand{\lx}{$L_{\text{X}}\xspace$}'
    printf,1,'\newcommand(\lxlim}{$L_{\text{X,lim}}\xspace$}'
    printf,1,'\newcommand{\lir}{$L_{\text{MIR}}$\xspace}'
    printf,1,'\newcommand{\lxir}{$L_{\text{X}}(L_{\text{MIR}})$\xspace}'
    printf,1,'\newcommand{\rlum}{$R_{L_{\text{X}}}$\xspace}'
    printf,1,'\newcommand{\llim}{$L_{{\text{X}}_{\text{lim}}}$\xspace}'
    printf,1,'%% Units'
    printf,1,'\newcommand{\kev}{\,{\text{keV}}\xspace}'
    printf,1,'\newcommand{\mum}{\,\mu{\text{m}}\xspace}'
    printf,1,'\newcommand{\cmcm}{\,{\text{cm}}^{-2}\xspace}'
    printf,1,'%\newcommand{\minus}{\scalebox{0.75}[1.0]{$-$}}'
    printf,1,'%%'
    printf,1,'%%'
    printf,1,'%%'
    printf,1,'%% Quality cuts for analysis'    
    printf,1,'\newcommand{\ngal}{'+commas(ngal)+'\xspace}'
    printf,1,'\newcommand{\nagn}{'+commas(nagn)+'\xspace}'
    printf,1,'\newcommand{\ndet}{'+commas(long(total(iidet)))+'\xspace}'
    printf,1,'\newcommand{\nnon}{'+commas(long(total(~iidet)))+'\xspace}'
    printf,1,'\newcommand{\nqgal}{'+commas(0)+'\xspace}'
    printf,1,'\newcommand{\nqagn}{'+commas(long(total(iiqual)))+'\xspace}'
    printf,1,'\newcommand{\nqdet}{'+commas(long(total(iiqual_det)))+'\xspace}'
    printf,1,'\newcommand{\nqnon}{'+commas(long(total(iiqual_non)))+'\xspace}'
    printf,1,'%% WISE AGNs'
    printf,1,'\newcommand{\ninfwac}{'+commas(long(total(iiwac)))+'\xspace}'
    printf,1,'\newcommand{\ndetwac}{'+commas(long(total(iiwac and iidet)))+'\xspace}'
    printf,1,'\newcommand{\nnonwac}{'+commas(long(total(iiwac and ~iidet)))+'\xspace}'
    printf,1,'\newcommand{\nqinfwac}{'+commas(long(total(iiwac and iiqual)))+'\xspace}'
    printf,1,'\newcommand{\nqdetwac}{'+commas(long(total(iiwac and iiqual_det)))+'\xspace}'
    printf,1,'\newcommand{\nqnonwac}{'+commas(long(total(iiwac and iiqual_non)))+'\xspace}'
    printf,1,'%% Chandra'
    printf,1,'\newcommand{\ninfcha}{'+commas(long(total(iiinf_cha)))+'\xspace}'
    printf,1,'\newcommand{\ndetcha}{'+commas(long(total(iiinf_cha and iidet_cha)))+'\xspace}'
    printf,1,'\newcommand{\nnoncha}{'+commas(long(total(iiinf_cha and ~iidet_cha)))+'\xspace}'
    printf,1,'\newcommand{\nqinfcha}{'+commas(long(total(strmatch(xdet,'CHA') or strmatch(xnon,'CHA'))))+'\xspace}'
    printf,1,'\newcommand{\nqdetcha}{'+commas(long(total(strmatch(xdet,'CHA'))))+'\xspace}'
    printf,1,'\newcommand{\nqnoncha}{'+commas(long(total(strmatch(xnon,'CHA'))))+'\xspace}'
    printf,1,'%% XMM'
    printf,1,'\newcommand{\ninfxmm}{'+commas(long(total(iiinf_xmm)))+'\xspace}'
    printf,1,'\newcommand{\ndetxmm}{'+commas(long(total(iiinf_xmm and iidet_xmm)))+'\xspace}'
    printf,1,'\newcommand{\nnonxmm}{'+commas(long(total(iiinf_xmm and ~iidet_xmm)))+'\xspace}'
    printf,1,'\newcommand{\nqinfxmm}{'+commas(long(total(strmatch(xdet,'XMM') or strmatch(xnon,'XMM'))))+'\xspace}'
    printf,1,'\newcommand{\nqdetxmm}{'+commas(long(total(strmatch(xdet,'XMM'))))+'\xspace}'
    printf,1,'\newcommand{\nqnonxmm}{'+commas(long(total(strmatch(xnon,'XMM'))))+'\xspace}'
    printf,1,'%% NuSTAR'
    printf,1,'\newcommand{\ninfnst}{'+commas(long(total(iiinf_nst)))+'\xspace}'
    printf,1,'\newcommand{\ndetnst}{'+commas(long(total(iiinf_nst and iidet_nst)))+'\xspace}'
    printf,1,'\newcommand{\nnonnst}{'+commas(long(total(iiinf_nst and ~iidet_nst)))+'\xspace}'
    printf,1,'\newcommand{\nqinfnst}{'+commas(long(total(strmatch(xdet,'NST') or strmatch(xnon,'NST'))))+'\xspace}'
    printf,1,'\newcommand{\nqdetnst}{'+commas(long(total(strmatch(xdet,'NST'))))+'\xspace}'
    printf,1,'\newcommand{\nqnonnst}{'+commas(long(total(strmatch(xnon,'NST'))))+'\xspace}'
    printf,1,'% Extraneous numbers'
    printf,1,'\newcommand{\nsrc}{'+commas(nsrc)+'\xspace}'
    printf,1,'\newcommand{\nqwac}{'+commas(long(total(iiqual and iiwac)))+'\xspace}'
    printf,1,'\newcommand{\nqres}{'+commas(long(total(iiqual and ~iiwac)))+'\xspace}'
    printf,1,'\newcommand{\nqdetres}{'+commas(long(total(iiqual_det and ~iiwac)))+'\xspace}'
    printf,1,'\newcommand{\nqnonres}{'+commas(long(total(iiqual_non and ~iiwac)))+'\xspace}'
    printf,1,'\newcommand{\nnqual}{'+commas(long(total(iiqual_non)))+'\xspace}'
    printf,1,'\newcommand{\pukidss}{'+strtrim(round(total(total(bin[where(strmatch(band,'UK?')),*],1) gt 0)/nsrc*100),2)+'\%\xspace}'
    printf,1,'\newcommand{\ptwom}{'+strtrim(round(total(total(bin[where(strmatch(band,'TWOM?')),*],1) gt 0)/nsrc*100),2)+'\%\xspace}'
    printf,1,'\newcommand{\pgalex}{'+strtrim(round(total(total(bin[where(strmatch(band,'GALEX?')),*],1) gt 0)/nsrc*100),2)+'\%\xspace}'
    printf,1,'\newcommand{\pcorr}{${\sim}$'+strtrim(round(total(iicorr ne 0)/nsrc*100),2)+'\%\xspace}'
    printf,1,'\newcommand{\lcorr}{'+string(rnd(alog10(mean(lcorr[where(iicorr)])),3),format='(d6.3)')+'\pm'+string(rnd(stddev(lcorr[where(iicorr)])/(alog(10.)*mean(lcorr[where(iicorr)])),3),format='(d5.3)')+'\xspace}'
    printf,1,'\newcommand{\pctwagn}{${\sim}$'+strtrim(round(total(iiwac and iiqual_non and llnon le (rl2nh(alog10(1.5e24),model='BORUS',/rl_out))[0])/total(iiwac and iiqual_non)*100),2)+'\%\xspace}'
    printf,1,'\newcommands{\nages}{'+commas(nex+nin)+'\xspace}'
    printf,1,'\newcommand{\nexages}{'+commas(nex)+'\xspace}'
    printf,1,'\newcommand{\ninages}{'+commas(nin)+'\xspace}'
    printf,1,'\newcommand{\nexoff}{'+commas(nexoff)+'\xspace}'
    printf,1,'\newcommand{\zexsig}{'+string(rnd(zexsig,3),format="(d5.3)")+'\xspace}'
    printf,1,'\newcommand{\zinsig}{'+string(rnd(zinsig,3),format="(d5.3)")+'\xspace}'
    close,1
endif


if keyword_set(record) then begin
    print, '    PRINTING RECORDS TABLE'

    OPENW,1,'tables/table1.tex'
    printf,1,'\begin{deluxetable*}{lrrrrrrrrrr}'
    printf,1,'\tablecaption{Number of sources with various telescope coverage.\label{sample}}'
    printf,1,'\tablehead{\colhead{} & \colhead{SDSS} & \colhead{XDQSOz} & \colhead{\wise} & \colhead{UKIDSS} & \colhead{2MASS} & \colhead{\galex} & \colhead{\chandra} & \colhead{\xmm} & \colhead{\nustar}}'
    printf,1,'\startdata'
    printf,1,'Initial & \nisdss & \nixdqso & \niwise & \niunwise & \nitwom & \nigalex & \nicha   & \nixmm   & \ninst   \\'
    printf,1,'Final   & \nsdss  & \nxdqso  & \nwise  & \nunwise  & \ntwom  & \ngalex  & \ninfcha & \ninfxmm & \ninfnst \\'
    printf,1,'\enddata'
    printf,1,'\tablecomments{Our initial sample of \ninit sources with optical and IR coverage, reduced to \nfinal sources for our final sample. All \wise photometry was replaced with forced photometry from unWISE where available (\niunwise initial, \nunwise final).}'
    printf,1,'\end{deluxetable*}'
    CLOSE,1
endif


if keyword_set(property_cuts) then begin
    print, '    PRINTING PROPERTY CUTS TABLE'
    
    openw,1,'tables/table2.tex'
    printf,1,'\begin{deluxetable*}{lrrrr}'
    printf,1,'\tablecaption{Various property selection criteria and effects on sample numbers.\label{prop_cuts}}'
    printf,1,'\tablehead{\colhead{Property Cut} & \colhead{$N_{\text{tot}}$} & \colhead{Tot. Loss} & \colhead{$N_{\text{cum}}$} & \colhead{Cum. Loss}}'
    printf,1,'\startdata'
    printf,1,'initial sample          & \ninit  & 0\%         & \ninit   & 0\%         \\'
    printf,1,'valid redshift          & \nzgood & \nzgoodtotf & \nrzgood & \nzgoodcumf \\'
    printf,1,'clean photometry        & \nclean & \ncleantotf & \nrclean & \ncleancumf \\'
    printf,1,'seven photometric bands & \nbands & \nbandstotf & \nrbands & \nbandscumf \\'
    printf,1,'four \wise bands        & \nfourw & \nfourwtotf & \nrfourw & \nfourwcumf \\'
    printf,1,'not removed by mask     & \nnomsk & \nnomsktotf & \nrnomsk & \nnomskcumf \\'
    printf,1,'not a duplicate         & \nnodup & \nnoduptotf & \nrnodup & \nnodupcumf \\'
    printf,1,'\enddata'
    printf,1,'\tablecomments{Individual and cumulative effect of each property cut on our sample. For each property cut, the first two subsequent columns list the number of sources which pass the designated cut and the fractional loss from our initial sample (\ninit sources). The last two columns track the cumulative effects of each property cut on sample size and our cumulative fractional loss.}'
    printf,1,'\end{deluxetable*}'
    close,1
endif


if keyword_set(xray) then begin
    print, '    PRINTING X-RAY SOURCE TABLE'
    ;; number of rows to print for table
    rows = 10
    ;inst = [['-','c'],['-','x'],['-','n']]
    iq = where(iiqual)
        
    id_str = string(objid,format='(i19)')
    ra_str = string(ra,format='(d8.4)')
    dec_str = string(dec,format='(d8.4)')
    ;fld_str = inst[iiinf_cha,0]+'/'+inst[iiinf_xmm,1]+'/'+inst[iiinf_nst,2]
    cha_str = string(loglx_cha>loglxlim_cha,format='(f5.2)')+'$\pm$'+string(e_loglx_cha>e_loglxlim_cha,format='(f4.2)')
    cha_str[where(loglxlim_cha ne -9999.,/null)] = '$<$'+cha_str[where(loglxlim_cha ne -9999.,/null)]
    cha_str[where(loglx_cha ne -9999.,/null)] = '   '+cha_str[where(loglx_cha ne -9999.,/null)]
    cha_str[where(strmatch(cha_str,'*\**'),/null)] = '          \nodata'
    xmm_str = string(loglx_xmm>loglxlim_xmm,format='(f5.2)')+'$\pm$'+string(e_loglx_xmm>e_loglxlim_xmm,format='(f4.2)')
    xmm_str[where(loglxlim_xmm ne -9999.,/null)] = '$<$'+xmm_str[where(loglxlim_xmm ne -9999.,/null)]
    xmm_str[where(loglx_xmm ne -9999.,/null)] = '   '+xmm_str[where(loglx_xmm ne -9999.,/null)]
    xmm_str[where(strmatch(xmm_str,'*\**'),/null)] = '          \nodata'
    nst_str = string(loglx_nst>loglxlim_nst,format='(f5.2)')+'$\pm$'+string(e_loglx_nst>e_loglxlim_nst,format='(f4.2)')
    nst_str[where(loglxlim_nst ne -9999.,/null)] = '$<$'+nst_str[where(loglxlim_nst ne -9999.,/null)]
    nst_str[where(loglx_nst ne -9999.,/null)] = '   '+nst_str[where(loglx_nst ne -9999.,/null)]
    nst_str[where(strmatch(nst_str,'*\**'),/null)] = '          \nodata'
    ;q_str = strarr(nsrc)+'0'
    ;q_str[where(iiqual_det)] = ' 1'
    ;q_str[where(iiqual_non)] = '-1'
    x_str = strarr(nsrc)
    x_str[where(iiqual_det)] = '    '+xdet[where(iiqual_det)]
    x_str[where(iiqual_non)] = '    '+xnon[where(iiqual_non)]
    x_str[where(x_str eq '')] = '\nodata'

    ;; all output lines
    line = id_str +' & ' +ra_str+' & '+dec_str+' & '+cha_str+' & '+xmm_str+' & '+nst_str+' & '+x_str+' \\'
    ;isort = sort(dec[iq])
    ;line = line[iq[isort]]
    line = line[iq]

    openw,1,'tables/table3.tex' 
    printf,1,'\begin{deluxetable*}{rrrrrrrrr}'
    printf,1,'\tablecolumns{9}'
    printf,1,'\tablecaption{Source positions and X-ray luminosity estimates by observatory.\label{xlum}}'
    printf,1,'\tablehead{\colhead{SDSS ObjID} & \colhead{RA} & \colhead{Dec} & \colhead{$\log L_{\text{X}}(\chandra)$} & \colhead{$\log L_{\text{X}}(\xmm)$} & \colhead{$\log L_{\text{X}}(\nustar)$} & \colhead{} \\' 
    printf,1,'\colhead{} & \colhead{[deg]} & \colhead{[deg]} & \colhead{[erg\,s$^{-1}$]} & \colhead{[erg\,s$^{-1}$]} & \colhead{[erg\,s$^{-1}$]} & \colhead{}}'
    printf,1,'\startdata'
    for i = 0,rows-1 do printf,1,line[i]
    printf,1,'\enddata'
    printf,1,'\tablecomments{Upper limit X-ray luminosity estimates are labeled by ``$<$''. The field chosen for our analysis is shown in the final column.}'
    printf,1,'%\tablerefs{}'
    printf,1,'\end{deluxetable*}'
    close,1


    ;; MACHINE READABLE TABLE
    ;; number of rows to print for table
    rows = nsrc
    inst = [['-','c'],['-','x'],['-','n']]
    iq = where(iiqual)
    
    id_str = string(objid,format='(i19)')
    ra_str = string(ra,format='(d12.8)')
    dec_str = string(dec,format='(d14.10)')
    fld_str = inst[iiinf_cha,0]+'/'+inst[iiinf_xmm,1]+'/'+inst[iiinf_nst,2]
    chalim_str = strarr(nsrc)+' '
    chalim_str[where(lxlim_cha)] = '<'
    chalx_str = string((lx_cha>lxlim_cha)*1e-7,format='(e10.4)')
    chalxerr_str = string((e_lx_cha>e_lxlim_cha)*1e-7,format='(e10.4)')
    xmmlim_str = strarr(nsrc)+' '
    xmmlim_str[where(lxlim_xmm)] = '<'
    xmmlx_str = string((lx_xmm>lxlim_xmm)*1e-7,format='(e10.4)')
    xmmlxerr_str = string((e_lx_xmm>e_lxlim_xmm)*1e-7,format='(e10.4)')
    nstlim_str = strarr(nsrc)+' '
    nstlim_str[where(lxlim_nst)] = '<'
    nstlx_str = string((lx_nst>lxlim_nst)*1e-7,format='(e10.4)')
    nstlxerr_str = string((e_lx_nst>e_lxlim_nst)*1e-7,format='(e10.4)')
    x_str = strarr(nsrc)
    x_str[where(iiqual_det)] = xdet[where(iiqual_det)]
    x_str[where(iiqual_non)] = xnon[where(iiqual_non)]
    x_str[where(x_str eq '')] = '   '

    line = id_str+' '+ra_str+' '+dec_str+' '+fld_str+' '+chalim_str+' '+chalx_str+' '+chalxerr_str+' '+xmmlim_str+' '+xmmlx_str+' '+xmmlxerr_str+' '+nstlim_str+' '+nstlx_str+' '+nstlxerr_str+' '+x_str
    isort = sort(dec)
    line = line[isort]
    
    print, '    PRINTING MACHINE READABLE FIT OUTPUT TABLE'

    table_name='tables/table3.txt' 
    openw,1,table_name

    printf,1,'Title: A Large Population of Luminous Active Galactic Nuclei Lacking X-ray '
    printf,1,'       Detections: Evidence for Heavy Obscuration? '
    printf,1,'Authors: Carroll C.M., Hickox R.C., Masini A., Lanz L., Assef R.J., '
    printf,1,'    Stern D., Chen C.-T.J., Ananna T.T. '
    printf,1,'Table: Source positions and X-ray luminosity estimates by observatory'
    printf,1,'================================================================================'
    printf,1,'Byte-by-byte Description of file: datafile3.txt'
    printf,1,'--------------------------------------------------------------------------------'
    printf,1,'   Bytes Format Units   Label        Explanations'
    printf,1,'--------------------------------------------------------------------------------'
    printf,1,'   1- 19 I19    ---     objidsdss    SDSS ObjID'
    printf,1,'  21- 32 F12.8  deg     RAdeg        Right Ascension in degrees (J2000)'
    printf,1,'  34- 47 F14.10 deg     DEdeg        Declination in degrees (J2000)'
    printf,1,'  49- 49 A1     ---   l_lxcha        Upper limit flag for Chandra source'
    printf,1,'  51- 60 E10.4  W       lxcha        ?=0.0000e+00 X-ray luminosity for Chandra source'
    printf,1,'  62- 71 E10.4  W     e_lxcha        ?=0.0000e+00 X-ray luminosity error for Chandra source'
    printf,1,'  73- 73 A1     ---   l_lxxmm        Upper limit flag for XMM source'
    printf,1,'  75- 84 E10.4  W       lxxmm        ?=0.0000e+00 X-ray luminosity for XMM source'
    printf,1,'  86- 95 E10.4  W     e_lxxmm        ?=0.0000e+00 X-ray luminosity error for XMM source'
    printf,1,'  97- 97 A1     ---   l_lxnst        Upper limit flag for NuSTAR source'
    printf,1,'  99-108 E10.4  W       lxnst        ?=0.0000e+00 X-ray luminosity for NuSTAR source'
    printf,1,' 110-119 E10.4  W     e_lxnst        ?=0.0000e+00 X-ray luminosity error for NuSTAR source'
    printf,1,' 121-123 A3     ---     field        Final X-ray field used for analysis (1)'
    printf,1,'--------------------------------------------------------------------------------'
    printf,1,'Note (1): Column 13 (field): Final X-ray field used for analysis is blank if not used.'
    printf,1,'--------------------------------------------------------------------------------'
    for i = 0,rows-1 do printf,1,line[i]
    close,1

endif


if keyword_set(sed) then begin
    print, '    PRINTING SED OUTPUT TABLE'
    ;; number of rows to print in table
    rows = 10
    iq = where(iiqual)
    
    ;ra_str = string(ra,format='(d10.6)')
    ;dec_str = string(dec,format='(d10.6)')
    z_str = string(rnd(z,3),format='(f5.3)')+'$\pm$'+strtrim(string(rnd(zerr,3),format='(f5.3)'),2)
    ebv_str = string(rnd(ebv<50.0,2),format='(f5.2)')+'$\pm$'+strtrim(string(rnd(e_ebv,2),format='(f5.2)'),2)
    agn_str = strmid(string(param[2,*],format="(e8.2)"),0,4)+'$\times$10$^{'+strmid(string(param[2,*],format="(e8.2)"),5,4)+'}$'
    ;agn_str[where(param[2,*] eq 0)] = '0.00                   '
    ell_str = strmid(string(param[3,*],format="(e8.2)"),0,4)+'$\times$10$^{'+strmid(string(param[3,*],format="(e8.2)"),5,4)+'}$'
    ;ell_str[where(param[3,*] eq 0)] = '0.00                   '
    sfg_str = strmid(string(param[4,*],format="(e8.2)"),0,4)+'$\times$10$^{'+strmid(string(param[4,*],format="(e8.2)"),5,4)+'}$'
    ;sfg_str[where(param[4,*] eq 0)] = '0.00                   '
    irr_str = strmid(string(param[5,*],format="(e8.2)"),0,4)+'$\times$10$^{'+strmid(string(param[5,*],format="(e8.2)"),5,4)+'}$'
    ;irr_str[where(param[5,*] eq 0)] = '0.00                   '
    chi_str = string(rnd(chi,2),format='(f6.2)')+' / '+strtrim(fix(dof),2)
    perc_str = string(perc_agn,format='(i3)')
    lir_str = string(rnd(loglir,2),format='(f5.2)')+'$\pm$'+strtrim(string(rnd(e_loglir,2),format='(f5.2)'),2)
    lir_str[where(loglir eq -9999.,/null)] = '\nodata       '
    wagn_str = strtrim(fix(iiwac),2)
    
    ;; all output lines
    ;line = id_str+' & '+ra_str+' & '+dec_str+' & '+z_str+' & '+ebv_str+' & '+chi_str+' & '+perc_str+' & '+lir_str+' & '+lx_str+' & '+wagn_str+' & '+field_str+' \\ '
    line = z_str+' & '+ebv_str+' & '+agn_str+' & '+ell_str+' & '+sfg_str+' & '+irr_str+' & '+chi_str+' & '+perc_str+' & '+lir_str+' \\ '
    ;isort = sort(dec[iq])
    ;line = line[iq[isort]]
    line = line[iq]
    
    OPENW,1,'tables/table4.tex' 
    printf,1,'\begin{deluxetable*}{rrrrrrrrr}'
    printf,1,'\tablecolumns{9}'
    printf,1,'\tablecaption{SED modeling output parameters.\label{sed_output}}'
    printf,1,'\tablehead{\colhead{$z$} & \colhead{$E(\bv)_{\text{AGN}}$} & \colhead{$C_{\text{AGN}}$} & \colhead{$C_{\text{ELL}}$} & \colhead{$C_{\text{SFG}}$} & \colhead{$C_{\text{IRR}}$} & \colhead{$\chi^2 /$ DoF} & \colhead{$P_{\text{AGN}}$} & \colhead{$\log L_{\text{MIR}}$} \\' 
    printf,1,'\colhead{} & \colhead{} & \colhead{} & \colhead{} & \colhead{} & \colhead{} & \colhead{} & \colhead{[\%]} & \colhead{[erg\,s$^{-1}$]}}'
    printf,1,'\startdata'
    for i = 0,rows-1 do printf,1,line[i]
    printf,1,'\enddata'
    printf,1,'\tablecomments{Modeling output parameters redshift, color excess \ebv, template normalizations (AGN, elliptical, star-forming, irregular/starburst), chi-squared per degrees of freedom, percentage of resampled realizations with AGN contribution ($C_{\text{AGN}} > 0$), and AGN $6\mum$ luminosity. Uncertainties on \ebv and $\log L_{\text{IR}}$ were estimated via the median absolute deviation of all source realizations containing AGN contribution.}'
    printf,1,'%\tablerefs{}'
    printf,1,'\end{deluxetable*}'
    close,1


    ;; MACHINE READABLE TABLE
    rows=nsrc
    inst = [['-','c'],['-','x'],['-','n']]
    
    z_str = string(rnd(z,3),format='(f5.3)')
    zerr_str = string(rnd(zerr,3),format='(f5.3)')
    ebv_str = string(rnd(ebv,2),format='(f5.2)')
    ebverr_str = string(rnd(e_ebv,2),format='(f5.2)')
    agn_str = string(param[2,*],format='(e13.7)')
    ell_str = string(param[3,*],format='(e13.7)')
    sfg_str = string(param[4,*],format='(e13.7)')
    irr_str = string(param[5,*],format='(e13.7)')
    chi_str = string(rnd(chi,2),format='(f10.2)')
    dof_str = string(fix(dof),format='(i2)')
    perc_str = string(perc_agn,format='(i3)')
    lir_str = string(lir*1e-7,format='(e10.4)')
    lirerr_str = string(e_lir*1e-7,format='(e10.4)')
    wagn_str = strtrim(fix(iiwac),2)
    
    line = z_str+' '+zerr_str+' '+ebv_str+' '+ebverr_str+' '+agn_str+' '+ell_str+' '+sfg_str+' '+irr_str+' '+chi_str+' '+dof_str+' '+perc_str+' '+lir_str+' '+lirerr_str+' '+wagn_str
    isort = sort(dec)
    line = line[isort]

    print, '    PRINTING MACHINE READABLE FIT OUTPUT TABLE'

    table_name='tables/table4.txt' 
    openw,1,table_name

    printf,1,'Title: A Large Population of Luminous Active Galactic Nuclei Lacking X-ray '
    printf,1,'       Detections: Evidence for Heavy Obscuration? '
    printf,1,'Authors: Carroll C.M., Hickox R.C., Masini A., Lanz L., Assef R.J., '
    printf,1,'    Stern D., Chen C.-T.J., Ananna T.T. '
    printf,1,'Table: SED modeling output parameters.'
    printf,1,'================================================================================'
    printf,1,'Byte-by-byte Description of file: datafile4.txt'
    printf,1,'--------------------------------------------------------------------------------'
    printf,1,'   Bytes Format Units   Label       Explanations'
    printf,1,'--------------------------------------------------------------------------------'
    printf,1,'   1-  5 F5.3   ---      z          Redshift (1)'
    printf,1,'   7- 11 F5.3   ---    e_z          Redshift error'
    printf,1,'  13- 17 F5.2   ---      ebv        Color excess (1)'
    printf,1,'  19- 23 F5.2   ---    e_ebv        Color excess error (2)'
    printf,1,'  25- 37 E13.7  ---      cagn       ?=0.0000000e+00 AGN template normalization (1)'
    printf,1,'  39- 51 E13.7  ---      cell       ?=0.0000000e+00 Elliptical template normalization (1)'
    printf,1,'  53- 65 E13.7  ---      csfg       ?=0.0000000e+00 Star-forming template normalization (1)'
    printf,1,'  67- 79 E13.7  ---      cirr       ?=0.0000000e+00 Irregular/starburst template normalization (1)'
    printf,1,'  81- 90 F10.2  ---      chi2       Chi-squared (1)'
    printf,1,'  92- 93 I2     ---    d_chi2       Degrees of Freedom'
    printf,1,'  95- 97 I3     %        pagn       Percent of realizations with AGN (2, 3)'
    printf,1,'  99-108 E10.4  W        lir        ?=0.0000e+00 AGN 6um IR luminosity (1)'
    printf,1,' 110-119 E10.4  W      e_lir        ?=0.0000e+00 AGN 6um IR luminosity error (2)'
    printf,1,' 121-122 I1     ---      wiseagn    WISE AGN flag (0=false, 1=true)'
    printf,1,'--------------------------------------------------------------------------------'
    printf,1,'Note (1): Modeling output parameters'
    printf,1,'Note (2): Uncertainties on ebv and lir were estimated via the median absolute '
    printf,1,'    deviation of all source realizations containing AGN contribution.'
    printf,1,'Note (3): Column 11 (pagn): Percent of realizations out of 1000 with cagn > 0.'
    printf,1,'--------------------------------------------------------------------------------'
    for i = 0,rows-1 do printf,1,line[i]
    close,1
endif


if keyword_set(quality_cuts) then begin
    print, '    PRINTING QUALITY CUTS TABLE'

    openw,1,'tables/table5.tex'
    printf,4,'\begin{deluxetable}{llrr}'
    printf,4,'\tablecolumns{3}'
    printf,4,'\tablecaption{Number of sources with various selections.\label{qual_cut}}'
    printf,4,'\tablehead{\colhead{Selection} & & \colhead{$N_{\text{final}}$} & \colhead{$N_{\text{subset}}$}}'
    printf,4,'\startdata'
    printf,4,'\multirow{4}{*}{Final sample} & SED Galaxy     & \ngal    & \nqgal    \\'
    printf,4,'                              & SED Galaxy+AGN & \nagn    & \nqagn    \\'
    printf,4,'                              & X-ray detected & \ndet    & \nqdet    \\'
    printf,4,'                              & X-ray non-det. & \nnon    & \nqnon    \\'
    printf,4,'                              \cline{2-4}'
    printf,4,'\multirow{3}{*}{\wise AGN}    & in catalog     & \ninfwac & \nqinfwac \\'
    printf,4,'                              & X-ray detected & \ndetwac & \nqdetwac \\'
    printf,4,'                              & X-ray non-det. & \nnonwac & \nqnonwac \\'
    printf,4,'                              \cline{2-4}'
    printf,4,'\multirow{3}{*}{\chandra}     & in X-ray field & \ninfcha & \nqinfcha \\'
    printf,4,'                              & X-ray detected & \ndetcha & \nqdetcha \\'
    printf,4,'                              & X-ray non-det. & \nnoncha & \nqnoncha \\'
    printf,4,'                              \cline{2-4}'
    printf,4,'\multirow{3}{*}{\xmm}         & in X-ray field & \ninfxmm & \nqinfxmm \\'
    printf,4,'                              & X-ray detected & \ndetxmm & \nqdetxmm \\'
    printf,4,'                              & X-ray non-det. & \nnonxmm & \nqnonxmm \\'
    printf,4,'                              \cline{2-4}'
    printf,4,'\multirow{3}{*}{\nustar}      & in X-ray field & \ninfnst & \nqinfnst \\'
    printf,4,'                              & X-ray detected & \ndetnst & \nqdetnst \\'
    printf,4,'                              & X-ray non-det. & \nnonnst & \nqnonnst \\'
    printf,4,'\enddata'
    printf,4,'\tablecomments{Our final sample of \nsrc sources, reduced to \nqagn sources for our analysis subset.}'
    printf,4,'%\tablerefs{}'
    printf,4,'\end{deluxetable}'
    close,1
endif


if keyword_set(xstack) then begin
    print, '    PRINTING X-RAY STACK TABLE'

    restore,'stack_output/stack_output.sav'
    
    ind = [0,1,2,3,5,4]
    exp_str = string(rnd(expv[ind,1]/1e6,2),format='(f5.2)')
    src_str = reform(string(fix(srcv[ind,1:2]),format='(i4)'),6,2)
    bg_str = reform(string(fix(bgv[ind,1:2]),format='(i4)'),6,2)
    net_str = reform(string(fix(net_srcv[ind,1:2]),format='(i4)'),6,2)
    flx_str = reform(string(rnd(strmid(strtrim(fluxv[ind,1:2],2),0,6),3),format='(f4.2)'),6,2)
    flxerr_str = reform(string(rnd(strmid(strtrim(flux_errv[ind,1:2],2),0,6),3),format='(f4.2)'),6,2)
    ;; move decimal place for errors with smaller powers
    flx_pow = float(strmid(fluxv[ind,1:2],13,15))
    err_pow = float(strmid(flux_errv[ind,1:2],13,15))
    iless = where(flx_pow/err_pow lt 1.)
    flxerr_str[iless] = string(rnd(flxerr_str[iless]/10.,3),format='(f4.2)')
    flx_str = '('+flx_str+'$\pm$'+flxerr_str+')$\times10^{'+strmid(fluxv[ind,1:2],13,15)+'}$'
    lx_str = reform(string(rnd(alog10(lxv[ind,1:2]),2),format='(f5.2)'),6,2)
    sn_str = fluxv[ind,1:2]/flux_errv[ind,1:2]
    lxerr_str = string(rnd((lxv[ind,1:2]/sn_str)/(alog(10.)*lxv[ind,1:2]),2),format='(f4.2)')
    lx_str = lx_str+'$\pm$'+lxerr_str
    ;lim_str = strarr(6,2)+'   '
    ;lim_str[-1] = '$<$'
    ;lx_str = lim_str+lx_str

    openw,1,'tables/table6.tex'
    printf,1,'\begin{deluxetable*}{lrrrrrrr}'
    printf,1,'\tablecolumns{8}'
    printf,1,'\tablecaption{X-ray stacking data and results.\label{stack_results}}'
    printf,1,'\tablehead{'
    printf,1,'\colhead{} & \colhead{$t_{\text{exp}}$} & \colhead{Energy} & \colhead{$N_{\text{src}}$} & \colhead{$N_{\text{bg}}$} & \colhead{$N_{\text{net}}$} & \colhead{$F_{\text{X}}$} & \colhead{$\log L_{\text{X}}$} \\'
    printf,1,'\colhead{} & \colhead{[Ms]} & \colhead{[keV]} & \colhead{} & \colhead{} & \colhead{} & \colhead{[erg\,s$^{-1}$\,cm$^{-2}$]} & \colhead{[erg\,s$^{-1}$]}'
    printf,1,'}'
    printf,1,'\startdata'
    printf,1,'\multirow{2}{*}{\wise AGN (X-ray det.)}     & \multirow{2}{*}{'+exp_str[0]+'} & 0.5--2 & '+src_str[0,0]+' & '+bg_str[0,0]+' & '+net_str[0,0]+' & '+flx_str[0,0]+' & '+lx_str[0,0]+' \\'
    printf,1,'                                            &                        & 2--7   & '+src_str[0,1]+' & '+bg_str[0,1]+' & '+net_str[0,1]+' & '+flx_str[0,1]+' & '+lx_str[0,1]+' \\'
    printf,1,'                                                                     \cline{2-8}'
    printf,1,'\multirow{2}{*}{\wise AGN (X-ray non-det.)} & \multirow{2}{*}{'+exp_str[1]+'} & 0.5--2 & '+src_str[1,0]+' & '+bg_str[1,0]+' & '+net_str[1,0]+' & '+flx_str[1,0]+' & '+lx_str[1,0]+' \\'
    printf,1,'                                            &                        & 2--7   & '+src_str[1,1]+' & '+bg_str[1,1]+' & '+net_str[1,1]+' & '+flx_str[1,1]+' & '+lx_str[1,1]+' \\'
    printf,1,'                                                                     \cline{2-8}'
    printf,1,'\multirow{2}{*}{Secondary (X-ray det.)}     & \multirow{2}{*}{'+exp_str[2]+'} & 0.5--2 & '+src_str[2,0]+' & '+bg_str[2,0]+' & '+net_str[2,0]+' & '+flx_str[2,0]+' & '+lx_str[2,0]+' \\'
    printf,1,'                                            &                        &   2--7 & '+src_str[2,1]+' & '+bg_str[2,1]+' & '+net_str[2,1]+' & '+flx_str[2,1]+' & '+lx_str[2,1]+' \\'
    printf,1,'                                                                     \cline{2-8}'
    printf,1,'\multirow{2}{*}{Secondary (X-ray non-det.)} & \multirow{2}{*}{'+exp_str[3]+'} & 0.5--2 & '+src_str[3,0]+' & '+bg_str[3,0]+' & '+net_str[3,0]+' & '+flx_str[3,0]+' & '+lx_str[3,0]+' \\'
    printf,1,'                                            &                        &   2--7 & '+src_str[3,1]+' & '+bg_str[3,1]+' & '+net_str[3,1]+' & '+flx_str[3,1]+' & '+lx_str[3,1]+' \\'
    printf,1,'                                                                     \cline{2-8}'
    printf,1,'\multirow{2}{*}{Removed AGN}                & \multirow{2}{*}{'+exp_str[4]+'} & 0.5--2 & '+src_str[4,0]+' & '+bg_str[4,0]+' & '+net_str[4,0]+' & '+flx_str[4,0]+' & '+lx_str[4,0]+' \\'
    printf,1,'                                            &                        & 2--7   & '+src_str[4,1]+' & '+bg_str[4,1]+' & '+net_str[4,1]+' & '+flx_str[4,1]+' & '+lx_str[4,1]+' \\'
    printf,1,'                                                                     \cline{2-8}'
    printf,1,'\multirow{2}{*}{SED galaxy}                 & \multirow{2}{*}{'+exp_str[5]+'} & 0.5--2 & '+src_str[5,0]+' & '+bg_str[5,0]+' & '+net_str[5,0]+' & '+flx_str[5,0]+' & '+lx_str[5,0]+' \\'
    printf,1,'                                            &                        &   2--7 & '+src_str[5,1]+' & '+bg_str[5,1]+' & '+net_str[5,1]+' & '+flx_str[5,1]+' & '+lx_str[5,1]+' \\'
    printf,1,'\enddata'
    printf,1,'\tablecomments{X-ray stacking results for energy ranges of \mbox{0.5--2} and \mbox{2--7 keV}. For each energy range, we present exposure time, photon counts, and flux and luminosity estimates. The flux and luminosity values listed in the \mbox{2--7 keV} energy range have been scaled to \mbox{2--10 keV} for direct comparison to the body of this work.}'
    printf,1,'%\tablerefs{}'
    printf,1,'\end{deluxetable*}'
    close,1
endif



END


