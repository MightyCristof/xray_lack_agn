PRO surv_analysis, ASURV = asurv


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


file_mkdir,'surv_analysis'
pushd,'surv_analysis'
spawn,'ln -s ~/Research/statistics/asurv-master/asurv .'
spawn,'ln -s /Users/ccarroll/IDLWorkspace/libraries/cmc/asurv_idl/all.com .'
spawn,'ln -s /Users/ccarroll/IDLWorkspace/libraries/cmc/asurv_idl/wagn.com .'
spawn,'ln -s /Users/ccarroll/IDLWorkspace/libraries/cmc/asurv_idl/ragn.com .'
spawn,'ln -s /Users/ccarroll/IDLWorkspace/libraries/cmc/asurv_idl/all2.com .'
spawn,'ln -s /Users/ccarroll/IDLWorkspace/libraries/cmc/asurv_idl/wagn2.com .'
spawn,'ln -s /Users/ccarroll/IDLWorkspace/libraries/cmc/asurv_idl/ragn2.com .'

;; prep det/non-det, luminosity ratio, obscured/unobscured/borderline
det = make_array(nsrc,/string,value='99')
det[where(iiqual_non)] = '-1'
det[where(iiqual_det)] = ' 0'
rl = dblarr(nsrc)-9999.
rl[where(iiqual_non)] = llnon[where(iiqual_non)]
rl[where(iiqual_det)] = lldet[where(iiqual_det)]
obsc = make_array(nsrc,/string,value='99')
obsc[where(ebv+e_ebv lt 0.15)] = ' 0'
obsc[where(ebv-e_ebv gt 0.15)] = ' 1'
obsc[where((ebv lt 0.15 and ebv+e_ebv gt 0.15) or (ebv gt 0.15 and ebv-e_ebv lt 0.15))] = ' 2'
obsc = obsc[where(iiqual)]

;; CUT TO RUN SURVIVAL ANALYSIS ON ONLY ANALYSIS SET
det = det[where(iiqual)]
rl = rl[where(iiqual)]
iwac = iidet_wac[where(iiqual)]
save,det,rl,iwac,file='surv_input.sav'

;; univariate, KM test
forprint,det,rl,format='2x,a2,2x,d7.4',textout='all.dat',/nocomment
forprint,det[where(iwac)],rl[where(iwac)],format='2x,a2,2x,d7.4',textout='wagn.dat',/nocomment
forprint,det[where(~iwac)],rl[where(~iwac)],format='2x,a2,2x,d7.4',textout='ragn.dat',/nocomment
;forprint,det[where(det eq 0)],rl[where(det eq 0)],format='2x,a2,2x,d7.4',textout='det1.dat',/nocomment
;forprint,det[where(det eq -1)],rl[where(det eq -1)],format='2x,a2,2x,d7.4',textout='non1.dat',/nocomment

;; univariate, two-sample test [0 unobscured, 1 obscured, 2 unsure]
forprint,obsc,det,rl,format='2x,a2,2x,a2,2x,d7.4',textout='all2.dat',/nocomment
forprint,obsc[where(iwac)],det[where(iwac)],rl[where(iwac)],format='2x,a2,2x,a2,2x,d7.4',textout='wagn2.dat',/nocomment
forprint,obsc[where(~iwac)],det[where(~iwac)],rl[where(~iwac)],format='2x,a2,2x,a2,2x,d7.4',textout='ragn2.dat',/nocomment
;forprint,obsc[where(det eq 0)],det[where(det eq 0)],rl[where(det eq 0)],format='2x,a2,2x,a2,2x,d7.4',textout='det2.dat',/nocomment
;forprint,obsc[where(det eq -1)],det[where(det eq -1)],rl[where(det eq -1)],format='2x,a2,2x,a2,2x,d7.4',textout='non2.dat',/nocomment


;; RUN ASURV
if keyword_set(asurv) then spawn,'./asurv'
    
;; RUN asurv_analysis.pro
surv_hist,/twosamp,/plt
restore,'surv_nh.sav'
;; RUN stck_nondet_flux

;; indicies for the four groups: (0) NON-DET WISE AGN, (1) DET WISE AGN, (2) NON-DET REMAINING SRCS, (3) DET REMAINING AGN SRCS
ind = ['ii = where(iidet_wac and xnon eq "CHA",nobj)','ii = where(iidet_wac and xdet eq "CHA",nobj)','ii = where(~iidet_wac and xnon eq "CHA",nobj)','ii = where(~iidet_wac and xdet eq "CHA",nobj)']

;; FLUX STACKING
stack_set = ['FX_WNON','FX_WDET','FX_RNON','FX_RDET']
for i = 0,3,2 do begin
    re = execute(ind[i])
    re = execute(stack_set[i]+' = stack_nondet_flux(loglir[ii],dl2[ii],wbinc,wd,det[where(iwac)],rl[where(iwac)])')
    re = execute(ind[i+1])
    re = execute(stack_set[i+1]+' = stack_det_flux(logfx[ii])')
endfor
sav_str = strjoin(stack_set,',')
re = execute('save,'+sav_str+',file="stack_fx.sav"')

;; OUTPUT SOURCES FOR X-RAY STACKING ANALYSIS
source_set = ['WAGN_NON','WAGN_DET','RAGN_NON','RAGN_DET']
vars = ['OBJID','RA','E_RA','DEC','E_DEC','Z','ZERR','LIR','E_LIR','EXP_CHA','FX','E_FX','TEXP_CHA','FXLIM','E_FXLIM']
source = {sdss_id:0ll,ra:0d,e_ra:0d,dec:0d,e_dec:0d,z:0d,zerr:0d,lir6:0d,e_lir6:0d,exp:0d,fx210:0d,e_fx210:0d,exp_lim:0d,fxlim:0d,e_fxlim:0d}
tags = tag_names(source)
for i = 0,n_elements(source_set)-1 do begin
    re = execute(ind[i])
    re = execute(source_set[i]+' = replicate(source,nobj)')
    for t = 0,n_elements(tags)-1 do re = execute(source_set[i]+'.'+tags[t]+' = '+vars[t]+'[ii]')
    if (i eq 0) then new = 1 else new = 0
    ;re = execute('mwrfits,'+source_set[i]+',"stack_data.fits",create=new')
endfor
stack_str = '{'+strjoin(source_set+':'+source_set,',')+'}'
re = execute('stacks = '+stack_str)
save,stacks,file='stack_data.sav'

popd


END



