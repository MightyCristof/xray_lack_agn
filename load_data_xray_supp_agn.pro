PRO load_data_xray_supp_agn


data_dir = '/Users/ccarroll/Research/sed_models/qso_sed/xray_supp_agn/run_20190715_1654'
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
load_vars,'xray_lum.sav','_agn_lum'     
load_vars,'cleaned_cha.sav','_clean_cha'   
load_vars,'cleaned_xmm.sav','_clean_xmm'   
load_vars,'cleaned_nst.sav','_clean_nst'   
load_vars,'quality_src.sav','_quality'     
popd


END

