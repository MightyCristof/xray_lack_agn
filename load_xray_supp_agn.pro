PRO load_xray_supp_agn


data_dir = '/Users/ccarroll/Research/sed_models/qso_sed/xray_supp_agn/run_20190909_1335'
print, '*********************************************************************************'
print, 'PATH: '+data_dir
print, '*********************************************************************************'

pushd,data_dir
load_vars,'fits.sav','_fits'            
load_vars,'infield_cha.sav','_inf_cha'   
load_vars,'infield_xmm.sav','_inf_xmm'   
load_vars,'infield_nst.sav','_inf_nst'   
load_vars,'detections_cha.sav','_det_cha'  
load_vars,'detections_xmm.sav','_det_xmm'  
load_vars,'detections_nst.sav','_det_nst'
load_vars,'detections_wac.sav','_det_wac'
load_vars,'softx_conversions.sav','_softx'   
load_vars,'catalog_flux_limits.sav','_fluxlim' 
load_vars,'../comp_*.sav','_comp'              
load_vars,'src_luminosity.sav','_agnlum'     
load_vars,'cleaned_cha.sav','_clean_cha'   
load_vars,'cleaned_xmm.sav','_clean_xmm'   
load_vars,'cleaned_nst.sav','_clean_nst'   
load_vars,'quality_src.sav','_quality'  
load_vars,'lum_ratio.sav','_lum_ratio'
load_vars,'nh_dist.sav','_nhdist'   
popd


END

;       common _fits        
;       common _inf_cha     
;       common _inf_xmm     
;       common _inf_nst     
;       common _det_cha     
;       common _det_xmm     
;       common _det_nst     
;       common _det_wac     
;       common _softx       
;       common _fluxlim     
;       common _comp        
;       common _agnlum      
;       common _clean_cha   
;       common _clean_xmm   
;       common _clean_nst   
;       common _quality     
;       common _lum_ratio   
;       common _nhdist      




   



