PRO xray_lacking_agn_inf


;; load fit output, template components, sample matched NuSTAR detections, NuSTAR "in_field" match
common _fits
common _xray_nst
common _xray_xmm
common _xray_cha

;; all COMMON blocks called that match sample size
block = '_'+['fits','xray_nst','xray_xmm','xray_cha']

;; source is in an X-ray field
iiinf = []
for i = 1,n_elements(block)-1 do begin
    tags = scope_varname(common=block[i])
    iiinf = [iiinf,tags[where(strmatch(tags,'iiinf*'),/null)]]
endfor
iiinf = strjoin(iiinf,' or ')
re = execute('iiinf = '+iiinf)
iinf = where(iiinf,nsrc)

;; number of sources
nobj = n_elements(ra)

;; handle only sources in an X-ray field
for i = 0,n_elements(block)-1 do begin
    ;; remove '_FIELD' suffix
    block_vars = scope_varname(common=block[i])
    vars = (strsplit(block_vars,'_field',/regex,/extract)).ToArray()
    ;; replace FITS output with '_INF' suffix
    if (block[i] eq '_fits') then vars = vars+'_inf'
    for v = 0,n_elements(block_vars)-1 do begin
        re = execute('sz = size('+block_vars[v]+')')
        case sz[0] of
            2: re = execute(vars[v]+' = '+block_vars[v]+'[*,iinf]')
            1: begin
                if (sz[1] eq nobj) then re = execute(vars[v]+' = '+block_vars[v]+'[iinf]') else $
                                        re = execute(vars[v]+' = '+block_vars[v])
               end
            else: print, 'CASE EXCEPTION: VARIABLE SIZE MISMATCH'
        endcase
    endfor
    ;; add NSRC to _INF_FITS
    if (block[i] eq '_FITS') then vars = ['NSRC',vars]
    sav_str = strjoin(vars,',')
    sav_file = 'infield'+strsplit(block[i],'_xray',/regex,/extract)+'.sav'
    re = execute('save,'+sav_str+',/compress,file=sav_file')        
endfor



END




