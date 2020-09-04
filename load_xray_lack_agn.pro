PRO load_xray_lack_agn, CURRENT = current


if keyword_set(current) then cd, current=data_dir else $
                             data_dir = '/Users/ccarroll/Research/projects/xray_lack_agn/workspace/run_20200507_final_fiore'

path = 'PATH: '+data_dir
str = strjoin(make_array(strlen(path),value='*'))
print, str
print, path
print, str

pushd,data_dir
load_vars,'fits.sav','_fits'
load_vars,'resamp.sav','_resamp'            
load_vars,'infield_cha.sav','_inf_cha'   
load_vars,'infield_xmm.sav','_inf_xmm'   
load_vars,'infield_nst.sav','_inf_nst'   
load_vars,'detections_cha.sav','_det_cha'  
load_vars,'detections_xmm.sav','_det_xmm'  
load_vars,'detections_nst.sav','_det_nst'
load_vars,'detections_wac.sav','_det_wac'
load_vars,'xband_conversions.sav','_xconv'   
load_vars,'catalog_flux_limits.sav','_fxlim' 
load_vars,'src_luminosities.sav','_agnlum'     
load_vars,'cleaned_cha.sav','_clean_cha'   
load_vars,'cleaned_xmm.sav','_clean_xmm'   
load_vars,'cleaned_nst.sav','_clean_nst'   
load_vars,'quality_src.sav','_quality'  
load_vars,'combined_lum.sav','_combined'
load_vars,'nh_dist.sav','_nhdist'   
popd


END






   



