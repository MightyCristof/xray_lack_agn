PRO print_samp_totals


common _tot
common _inf_cha
common _inf_xmm
common _inf_nst

table_name='samp_totals.tex' 
OPENW,1,table_name

printf,1,'%% Sample Numbers (initial)'
printf,1,'\newcommand{\ninit}{'+commas(n_elements(obs))+'\xspace}'
printf,1,'\newcommand{\nisdss}{'+commas(long(total(iisdss)))+'\xspace}'
printf,1,'\newcommand{\nizsupp}{'+commas(long(total(iizsupp)))+'\xspace}'
printf,1,'\newcommand{\nixdqso}{'+commas(long(total(iixdqso and ~iisdss)))+'\xspace}'
printf,1,'\newcommand{\niwise}{'+commas(long(total(iiwise)))+'\xspace}'
printf,1,'\newcommand{\niunwise}{'+commas(long(total(iiunwise)))+'\xspace}'
printf,1,'\newcommand{\niukidss}{'+commas(long(total(iiukidss)))+'\xspace}'
printf,1,'\newcommand{\nitwom}{'+commas(long(total(iitwom)))+'\xspace}'
printf,1,'\newcommand{\nigalex}{'+commas(long(total(iigalex)))+'\xspace}'
printf,1,'\newcommand{\nicha}{'+commas(long(total(iiinf_cha)))+'\xspace}'
printf,1,'\newcommand{\nixmm}{'+commas(long(total(iiinf_xmm)))+'\xspace}'
printf,1,'\newcommand{\ninst}{'+commas(long(total(iiinf_nst)))+'\xspace}'
printf,1,'\newcommand{\nichaxmm}{'+commas(long(total(iiinf_cha and iiinf_xmm)))+'\xspace}'
printf,1,'\newcommand{\nichanst}{'+commas(long(total(iiinf_cha and iiinf_nst)))+'\xspace}'
printf,1,'\newcommand{\nixmmnst}{'+commas(long(total(iiinf_xmm and iiinf_nst)))+'\xspace}'

iiaccept = obs.iiaccept

printf,1,'%% Sample Numbers (final)'
printf,1,'\newcommand{\nfnit}{'+commas(long(total(iiaccept)))+'\xspace}'
printf,1,'\newcommand{\nfsdss}{'+commas(long(total(iiaccept and iisdss)))+'\xspace}'
printf,1,'\newcommand{\nfzsupp}{'+commas(long(total(iiaccept and iizsupp)))+'\xspace}'
printf,1,'\newcommand{\nfxdqso}{'+commas(long(total(iiaccept and iixdqso and ~iisdss)))+'\xspace}'
printf,1,'\newcommand{\nfwise}{'+commas(long(total(iiaccept and iiwise)))+'\xspace}'
printf,1,'\newcommand{\nfunwise}{'+commas(long(total(iiaccept and iiunwise)))+'\xspace}'
printf,1,'\newcommand{\nfukidss}{'+commas(long(total(iiaccept and iiukidss)))+'\xspace}'
printf,1,'\newcommand{\nftwom}{'+commas(long(total(iiaccept and iitwom)))+'\xspace}'
printf,1,'\newcommand{\nfgalex}{'+commas(long(total(iiaccept and iigalex)))+'\xspace}'
printf,1,'\newcommand{\nfcha}{'+commas(long(total(iiaccept and iiinf_cha)))+'\xspace}'
printf,1,'\newcommand{\nfxmm}{'+commas(long(total(iiaccept and iiinf_xmm)))+'\xspace}'
printf,1,'\newcommand{\nfnst}{'+commas(long(total(iiaccept and iiinf_nst)))+'\xspace}'
printf,1,'\newcommand{\nfchaxmm}{'+commas(long(total(iiaccept and iiinf_cha and iiinf_xmm)))+'\xspace}'
printf,1,'\newcommand{\nfchanst}{'+commas(long(total(iiaccept and iiinf_cha and iiinf_nst)))+'\xspace}'
printf,1,'\newcommand{\nfxmmnst}{'+commas(long(total(iiaccept and iiinf_xmm and iiinf_nst)))+'\xspace}'

CLOSE,1


END





