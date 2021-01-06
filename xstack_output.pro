PRO xstack_output


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
common _surv
common _rsurv

;; indicies for the four groups: (0) NON-DET WISE AGN, (1) DET WISE AGN, (2) NON-DET REMAINING SRCS, (3) DET REMAINING AGN SRCS
;restore,'rsurv_input.sav'
ind = ['ii = where(iiwac and xnon eq "CHA",nobj)','ii = where(iiwac and xdet eq "CHA",nobj)','ii = where(~iiwac and xnon eq "CHA",nobj)','ii = where(~iiwac and xdet eq "CHA",nobj)']


;; FLUX STACKING (but actually this time!)
stack_set = ['FX_WNON','FX_WDET','FX_RNON','FX_RDET']
;; X-ray non-detected
iwn = where(iiwac and xnon eq "CHA",nobj)
if (nobj gt 0) then fx_wnon = stack_nondet_flux(loglir[iwn],dl2[iwn],reverse(-wbinc),reverse(wd),xdw,rlw)
irn = where(~iiwac and xnon eq "CHA",nobj)
if (nobj gt 0) then fx_rnon = stack_nondet_flux(loglir[irn],dl2[irn],reverse(-rbinc),reverse(rd),xdr,rlr)
;; X-ray detected
iwd = where(iiwac and xdet eq 'CHA',nobj)
if (nobj gt 0) then fx_wdet = stack_det_flux(logfx[iwd])
ird = where(~iiwac and xdet eq 'CHA',nobj)
if (nobj gt 0) then fx_rdet = stack_det_flux(logfx[ird])
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


END