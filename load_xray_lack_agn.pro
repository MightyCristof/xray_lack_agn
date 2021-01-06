PRO load_xray_lack_agn, subdir, $
                        PHOT = phot


if (n_elements(subdir) eq 0) then cd,current=path else $
                                  path = file_search(subdir,/fully)
path += '/'

dir = 'LOAD: '+path
str = strjoin(make_array(strlen(dir),value='*'))
print, ''
print, str
print, dir
print, str
print, ''

if keyword_set(phot) then load_vars,'../data_prep/sdssXwise_initial.sav','_phot'
load_vars,'fits.sav','_fits'
load_vars,'resamp.sav','_resamp' 
load_vars,'../data_prep/comp*.sav','_comp'
cd,path
load_vars,path+'infield_cha.sav','_inf_cha'   
load_vars,path+'infield_xmm.sav','_inf_xmm'   
load_vars,path+'infield_nst.sav','_inf_nst'   
load_vars,path+'detections_cha.sav','_det_cha'  
load_vars,path+'detections_xmm.sav','_det_xmm'  
load_vars,path+'detections_nst.sav','_det_nst'
load_vars,path+'detections_wac.sav','_wac'
load_vars,path+'xband_conversions.sav','_xconv'   
load_vars,path+'catalog_flux_limits.sav','_fxlim' 
load_vars,path+'src_luminosities.sav','_agnlum'     
load_vars,path+'quality_src.sav','_quality'  
load_vars,path+'combined_lum.sav','_combined'
load_vars,path+'nh_dist.sav','_nhdist'   
load_vars,path+'surv_anal.sav','_surv'
load_vars,path+'stack_fx.sav','_xstack'

END






   



