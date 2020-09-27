PRO print_property_cuts


common _tot


ii = ['IIZGOOD','IICLEAN','IIBANDS','IIFOURW','IINOMSK','IINODUP']

table_name='propert_cuts.tex' 
OPENW,1,table_name
for i = 0,n_elements(ii)-1 do begin
    com = 'n'+strlowcase(strsplit(ii[i],'II',/regex,/extract))
    rcom ='nr'+strlowcase(strsplit(ii[i],'II',/regex,/extract))
    re = execute('tot = total('+ii[i]+')')
    totf = (1.-tot/nsrcs)*100.
    iicum = strjoin(ii[0:i],' and ')
    re = execute('cum = total('+iicum+')')
    cumf = (1.-cum/nsrcs)*100.
    printf,1, '\newcommand{\'+com+'}{'+commas(long(tot))+'\xspace}'
    printf,1, '\newcommand{\'+com+'tot}{'+string(rnd(totf,2),format="(d4.1)")+'\%\xspace}'
    printf,1, '\newcommand{\'+rcom+'}{'+commas(long(cum))+'\xspace}'
    if (i eq n_elements(ii)-1) then printf,1, '\newcommand{\'+com+'cum}{$>$0.01\%\xspace}' else $
                                    printf,1, '\newcommand{\'+com+'cum}{'+string(rnd(cumf,2),format="(d4.1)")+'\%\xspace}'
endfor
CLOSE,1


END











