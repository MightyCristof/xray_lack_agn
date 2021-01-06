PRO print_samp_totals


common _init
common _inf_cha
common _inf_xmm
common _inf_nst


nsrcs = n_elements(obs)
table_name='samp_totals.tex' 
OPENW,1,table_name

printf,1,'%% Sample Numbers (initial)'
printf,1,'\newcommand{\ninit}{'+commas(n_elements(obs))+'\xspace}'
printf,1,'\newcommand{\nisdss}{'+commas(long(total(iisdss)))+'\xspace}'
printf,1,'\newcommand{\nizsupp}{'+commas(long(total(iizsupp)))+'\xspace}'
printf,1,'\newcommand{\nixdqso}{'+commas(long(total(iixdqso and ~iisdss)))+'\xspace}'
printf,1,'\newcommand{\nioptical}{'+commas(long(total(iisdss or iixdqso)))+'\xspace}'
printf,1,'\newcommand{\niwise}{'+commas(long(total(iiwise)))+'\xspace}'
printf,1,'\newcommand{\niunwise}{'+commas(long(total(iiunwise)))+'\xspace}'
printf,1,'\newcommand{\nimir}{'+commas(long(total(iiwise or iiunwise)))+'\xspace}'
printf,1,'\newcommand{\niukidss}{'+commas(long(total(iiukidss)))+'\xspace}'
printf,1,'\newcommand{\nitwom}{'+commas(long(total(iitwom)))+'\xspace}'
printf,1,'\newcommand{\ninir}{'+commas(long(total(iiukidss or iitwom)))+'\xspace}'
printf,1,'\newcommand{\nigalex}{'+commas(long(total(iigalex)))+'\xspace}'
printf,1,'\newcommand{\nixray}{'+commas(long(total(iiinf_cha or iiinf_xmm or iiinf_nst)))+'\xspace}'
printf,1,'\newcommand{\nicha}{'+commas(long(total(iiinf_cha)))+'\xspace}'
printf,1,'\newcommand{\nixmm}{'+commas(long(total(iiinf_xmm)))+'\xspace}'
printf,1,'\newcommand{\ninst}{'+commas(long(total(iiinf_nst)))+'\xspace}'
printf,1,'\newcommand{\nichaxmm}{'+commas(long(total(iiinf_cha and iiinf_xmm)))+'\xspace}'
printf,1,'\newcommand{\nichanst}{'+commas(long(total(iiinf_cha and iiinf_nst)))+'\xspace}'
printf,1,'\newcommand{\nixmmnst}{'+commas(long(total(iiinf_xmm and iiinf_nst)))+'\xspace}'

ii = ['IIZGOOD','IICLEAN','IIBANDS','IIFOURW','IINOMSK','IINODUP']
printf,1,'%% Property cuts'
for i = 0,n_elements(ii)-1 do begin
    com = 'n'+strlowcase(strsplit(ii[i],'II',/regex,/extract))
    rcom ='nr'+strlowcase(strsplit(ii[i],'II',/regex,/extract))
    re = execute('tot = total('+ii[i]+')')
    totf = (1.-tot/nsrcs)*100.
    iicum = strjoin(ii[0:i],' and ')
    re = execute('cum = total('+iicum+')')
    cumf = (1.-cum/nsrcs)*100.
    printf,1, '\newcommand{\'+com+'}{'+commas(long(tot))+'\xspace}'
    printf,1, '\newcommand{\'+com+'totf}{'+string(rnd(totf,2),format="(d4.1)")+'\%\xspace}'
    printf,1, '\newcommand{\'+rcom+'}{'+commas(long(cum))+'\xspace}'
    printf,1, '\newcommand{\'+com+'cumf}{'+string(rnd(cumf,2),format="(d4.1)")+'\%\xspace}'
endfor

iiaccept = obs.iiaccept
printf,1,'%% Sample Numbers (final)'
printf,1,'\newcommand{\nfinal}{'+commas(long(total(iiaccept)))+'\xspace}'
printf,1,'\newcommand{\nsdss}{'+commas(long(total(iiaccept and iisdss)))+'\xspace}'
printf,1,'\newcommand{\nzsupp}{'+commas(long(total(iiaccept and iizsupp)))+'\xspace}'
printf,1,'\newcommand{\nxdqso}{'+commas(long(total(iiaccept and iixdqso and ~iisdss)))+'\xspace}'
printf,1,'\newcommand{\noptical}{'+commas(long(total(iiaccept and (iisdss or iixdqso))))+'\xspace}'
printf,1,'\newcommand{\nwise}{'+commas(long(total(iiaccept and iiwise)))+'\xspace}'
printf,1,'\newcommand{\nunwise}{'+commas(long(total(iiaccept and iiunwise)))+'\xspace}'
printf,1,'\newcommand{\nmir}{'+commas(long(total(iiaccept and (iiwise or iiunwise))))+'\xspace}'
printf,1,'\newcommand{\nukidss}{'+commas(long(total(iiaccept and iiukidss)))+'\xspace}'
printf,1,'\newcommand{\ntwom}{'+commas(long(total(iiaccept and iitwom)))+'\xspace}'
printf,1,'\newcommand{\nnir}{'+commas(long(total(iiaccept and (iiukidss or iitwom))))+'\xspace}'
printf,1,'\newcommand{\ngalex}{'+commas(long(total(iiaccept and iigalex)))+'\xspace}'
;printf,1,'\newcommand{\ninfcha}{'+commas(long(total(iiaccept and iiinf_cha)))+'\xspace}'
;printf,1,'\newcommand{\ninfxmm}{'+commas(long(total(iiaccept and iiinf_xmm)))+'\xspace}'
;printf,1,'\newcommand{\ninfnst}{'+commas(long(total(iiaccept and iiinf_nst)))+'\xspace}'
;printf,1,'\newcommand{\ninfchaxmm}{'+commas(long(total(iiaccept and iiinf_cha and iiinf_xmm)))+'\xspace}'
;printf,1,'\newcommand{\ninfchanst}{'+commas(long(total(iiaccept and iiinf_cha and iiinf_nst)))+'\xspace}'
;printf,1,'\newcommand{\ninfxmmnst}{'+commas(long(total(iiaccept and iiinf_xmm and iiinf_nst)))+'\xspace}'

CLOSE,1


END



















;iiaccept = obs.iiaccept
;
;printf,1,'%% Sample Numbers (final)'
;printf,1,'\newcommand{\nfnit}{'+commas(long(total(iiaccept)))+'\xspace}'
;printf,1,'\newcommand{\nfsdss}{'+commas(long(total(iiaccept and iisdss)))+'\xspace}'
;printf,1,'\newcommand{\nfzsupp}{'+commas(long(total(iiaccept and iizsupp)))+'\xspace}'
;printf,1,'\newcommand{\nfxdqso}{'+commas(long(total(iiaccept and iixdqso and ~iisdss)))+'\xspace}'
;printf,1,'\newcommand{\nfwise}{'+commas(long(total(iiaccept and iiwise)))+'\xspace}'
;printf,1,'\newcommand{\nfunwise}{'+commas(long(total(iiaccept and iiunwise)))+'\xspace}'
;printf,1,'\newcommand{\nfukidss}{'+commas(long(total(iiaccept and iiukidss)))+'\xspace}'
;printf,1,'\newcommand{\nftwom}{'+commas(long(total(iiaccept and iitwom)))+'\xspace}'
;printf,1,'\newcommand{\nfgalex}{'+commas(long(total(iiaccept and iigalex)))+'\xspace}'
;printf,1,'\newcommand{\nfcha}{'+commas(long(total(iiaccept and iiinf_cha)))+'\xspace}'
;printf,1,'\newcommand{\nfxmm}{'+commas(long(total(iiaccept and iiinf_xmm)))+'\xspace}'
;printf,1,'\newcommand{\nfnst}{'+commas(long(total(iiaccept and iiinf_nst)))+'\xspace}'
;printf,1,'\newcommand{\nfchaxmm}{'+commas(long(total(iiaccept and iiinf_cha and iiinf_xmm)))+'\xspace}'
;printf,1,'\newcommand{\nfchanst}{'+commas(long(total(iiaccept and iiinf_cha and iiinf_nst)))+'\xspace}'
;printf,1,'\newcommand{\nfxmmnst}{'+commas(long(total(iiaccept and iiinf_xmm and iiinf_nst)))+'\xspace}'
;
;CLOSE,1
