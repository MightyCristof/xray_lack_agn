PRO combine_nustar_catalogs


field = file_search(/test_dir)
field = field[where(~strmatch(field,'ECDFS'),/null)]

;; restore NuSTAR fields and normalize data
for i = 0,n_elements(field)-1 do begin
    pushd, field[i]
    re = execute(field[i]+' = mrdfits("nu*.fits",1)')
    re = execute(field[i]+' = standardize_nustar_catalogs('+field[i]+',field[i])')
    popd
endfor

;; number of sources
str = strjoin('n_elements('+field+')',"+")
re = execute('len = '+str)

;; combined source structure
nst= {RA:0., DEC:0., $
      SBF:0., HBF:0., FBF:0., $
      E_SBF:0., E_HBF:0., E_FBF:0., $
      SEXP:0., HEXP:0., FEXP:0., $
      NST_FIELD:'' $
      }
nst = replicate(nst,len)
nutags = tag_names(nst)

;; individual field tags to match nustar tags
xtags = [['CXOXMM_RA','CXOXMM_DEC','SB_FLUX','HB_FLUX','BB_FLUX','SB_FLUX_ERROR','HB_FLUX_ERROR','BB_FLUX_ERROR','SB_EXPOSURE','HB_EXPOSURE','BB_EXPOSURE'], $
         ;['CTRPRT_RA','CTRPRT_DEC','SB_FLUX','HB_FLUX','FB_FLUX','SB_FLUX_ERROR','HB_FLUX_ERROR','FB_FLUX_ERROR','SB_EXPOSURE','HB_EXPOSURE','FB_EXPOSURE'], $
         ['SOFTXRAY_CTRPART_RA','SOFTXRAY_CTRPART_DEC','SB_FLUX','HB_FLUX','FB_FLUX','SB_FLUX_ERROR','HB_FLUX_ERROR','FB_FLUX_ERROR','SB_EXPOSURE','HB_EXPOSURE','FB_EXPOSURE'], $
         ['SOFTXRAY_CTRPART_RA','SOFTXRAY_CTRPART_DEC','SB_FLUX','HB_FLUX','FB_FLUX','SB_FLUX_ERROR','HB_FLUX_ERROR','FB_FLUX_ERROR','SB_EXPOSURE','HB_EXPOSURE','FB_EXPOSURE'], $
         ['XMMRA','XMMDE','_3_8F','_8_24F','_3_24F','E_POS_3_8F','E_POS_8_24F','E_POS_3_24F','_3_8EXP','_8_24EXP','_3_24EXP'] $
         ]
nxtags = n_elements(xtags[*,0])

;; match tags and combine fields
ct = 0
for f = 0,n_elements(field)-1 do begin
    re = execute('ctfield = ct+n_elements('+field[f]+')')
    re = execute('tags = tag_names('+field[f]+')')
    for t = 0,nxtags-1 do begin
        itag = where(strmatch(tags,xtags[t,f]))
        re = execute('nst[ct:ctfield-1].(t) = '+field[f]+'.(itag)')
    endfor
    nst[ct:ctfield-1].(nxtags) = field[f]
    ct = ctfield
endfor

mwrfits,nst,'nustar_combined_catalogs.fits'


END



