PRO compute_nh_distribution


common _quality
common _lratio


;; NH distribution from all detections/non-detections
nhxdet = ll2nh(llxdet,'2-10')
nhxnon = ll2nh(llxnon,'2-10')






sav_vars = ['NHXDET','NHXNON']
sav_inds = []

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="nh_dist.sav"')


END




