PRO compute_nh_distribution


common _wac
common _quality
common _combined


;; histogram bin size
nhbinsz = 0.2

;; POWER-LAW MODEL
;; NH distribution from all det_powerections/non_power-det_powerections
nhdet_power = rl2nh(lldet,model='power')
e_nhdet_power = nhdet_power
e_nhdet_power[where(iiqual_det)] = resamp_rlnh_err(ret='NH_MAD',model='POWER')
nhnon_power = rl2nh(llnon,model='power')


yhist_det_power = histogram(nhdet_power[where(iiqual_det)],locations=xhist_det_power,bin=nhbinsz)
yhist_non_power = histogram(nhnon_power[where(iiqual_non)],locations=xhist_non_power,bin=nhbinsz)

xhist_det_power = [xhist_det_power[0]-nhbinsz,xhist_det_power,xhist_det_power[-1]+nhbinsz]
xhist_non_power = [xhist_non_power[0]-nhbinsz,xhist_non_power,xhist_non_power[-1]+nhbinsz]
yhist_det_power = [0.,yhist_det_power,0.]
yhist_non_power = [0.,yhist_non_power,0.]

sav_vars = ['NHBINSZ','NHDET_POWER','E_NHDET_POWER','NHNON_POWER', $
                      'YHIST_DET_POWER','XHIST_DET_POWER','YHIST_NON_POWER','XHIST_NON_POWER']
sav_inds = []

;; construct distributions by type WISE AGN/Remaining
;; indices of WISE AGN detections/non-detections
iiwdet_power = iiqual_det and iiwac
iiwnon_power = iiqual_non and iiwac
;; indices of remaining detections/non-detections
iirdet_power = iiqual_det and ~iiwac
iirnon_power = iiqual_non and ~iiwac

yhist_wdet_power = histogram(nhdet_power[where(iiwdet_power)],locations=xhist_wdet_power,bin=nhbinsz)
yhist_wnon_power = histogram(nhnon_power[where(iiwnon_power)],locations=xhist_wnon_power,bin=nhbinsz)
yhist_rdet_power = histogram(nhdet_power[where(iirdet_power)],locations=xhist_rdet_power,bin=nhbinsz)
yhist_rnon_power = histogram(nhnon_power[where(iirnon_power)],locations=xhist_rnon_power,bin=nhbinsz)

xhist_wdet_power = [xhist_wdet_power[0]-nhbinsz,xhist_wdet_power,xhist_wdet_power[-1]+nhbinsz]
xhist_wnon_power = [xhist_wnon_power[0]-nhbinsz,xhist_wnon_power,xhist_wnon_power[-1]+nhbinsz]
xhist_rdet_power = [xhist_rdet_power[0]-nhbinsz,xhist_rdet_power,xhist_rdet_power[-1]+nhbinsz]
xhist_rnon_power = [xhist_rnon_power[0]-nhbinsz,xhist_rnon_power,xhist_rnon_power[-1]+nhbinsz]
yhist_wdet_power = [0.,yhist_wdet_power,0.]
yhist_wnon_power = [0.,yhist_wnon_power,0.]
yhist_rdet_power = [0.,yhist_rdet_power,0.]
yhist_rnon_power = [0.,yhist_rnon_power,0.]

sav_vars = [sav_vars,'YHIST_WDET_POWER','XHIST_WDET_POWER','YHIST_WNON_POWER','XHIST_WNON_POWER', $
                     'YHIST_RDET_POWER','XHIST_RDET_POWER','YHIST_RNON_POWER','XHIST_RNON_POWER']
sav_inds = [sav_inds,'IIWDET_POWER','IIWNON_POWER','IIRDET_POWER','IIRNON_POWER']


;; BORUS MODEL
;; NH distribution from all det_powerections/non_power-det_powerections
nhdet_borus = rl2nh(lldet,model='BORUS')
e_nhdet_borus = nhdet_borus
e_nhdet_borus[where(iiqual_det)] = resamp_rlnh_err(ret='NH_MAD',model='BORUS')
nhnon_borus = rl2nh(llnon,model='BORUS')

yhist_det_borus = histogram(nhdet_borus[where(iiqual_det)],locations=xhist_det_borus,bin=nhbinsz)
yhist_non_borus = histogram(nhnon_borus[where(iiqual_non)],locations=xhist_non_borus,bin=nhbinsz)

xhist_det_borus = [xhist_det_borus[0]-nhbinsz,xhist_det_borus,xhist_det_borus[-1]+nhbinsz]
xhist_non_borus = [xhist_non_borus[0]-nhbinsz,xhist_non_borus,xhist_non_borus[-1]+nhbinsz]
yhist_det_borus = [0.,yhist_det_borus,0.]
yhist_non_borus = [0.,yhist_non_borus,0.]

sav_vars = [sav_vars,'NHDET_BORUS','E_NHDET_BORUS','NHNON_BORUS', $
                     'YHIST_DET_BORUS','XHIST_DET_BORUS','YHIST_NON_BORUS','XHIST_NON_BORUS']
sav_inds = [sav_inds]

;; construct distributions by type WAC/Remaining
;; indices of WISE AGN det_borusections/non_borus-det_borusections
iiwdet_borus = iiqual_det and iiwac
iiwnon_borus = iiqual_non and iiwac
;; indices of remaining det_borusections/non_borus-det_borusections
iirdet_borus = iiqual_det and ~iiwac
iirnon_borus = iiqual_non and ~iiwac

yhist_wdet_borus = histogram(nhdet_borus[where(iiwdet_borus)],locations=xhist_wdet_borus,bin=nhbinsz)
yhist_wnon_borus = histogram(nhnon_borus[where(iiwnon_borus)],locations=xhist_wnon_borus,bin=nhbinsz)
yhist_rdet_borus = histogram(nhdet_borus[where(iirdet_borus)],locations=xhist_rdet_borus,bin=nhbinsz)
yhist_rnon_borus = histogram(nhnon_borus[where(iirnon_borus)],locations=xhist_rnon_borus,bin=nhbinsz)

xhist_wdet_borus = [xhist_wdet_borus[0]-nhbinsz,xhist_wdet_borus,xhist_wdet_borus[-1]+nhbinsz]
xhist_wnon_borus = [xhist_wnon_borus[0]-nhbinsz,xhist_wnon_borus,xhist_wnon_borus[-1]+nhbinsz]
xhist_rdet_borus = [xhist_rdet_borus[0]-nhbinsz,xhist_rdet_borus,xhist_rdet_borus[-1]+nhbinsz]
xhist_rnon_borus = [xhist_rnon_borus[0]-nhbinsz,xhist_rnon_borus,xhist_rnon_borus[-1]+nhbinsz]
yhist_wdet_borus = [0.,yhist_wdet_borus,0.]
yhist_wnon_borus = [0.,yhist_wnon_borus,0.]
yhist_rdet_borus = [0.,yhist_rdet_borus,0.]
yhist_rnon_borus = [0.,yhist_rnon_borus,0.]

sav_vars = [sav_vars,'YHIST_WDET_BORUS','XHIST_WDET_BORUS','YHIST_WNON_BORUS','XHIST_WNON_BORUS', $
                     'YHIST_RDET_BORUS','XHIST_RDET_BORUS','YHIST_RNON_BORUS','XHIST_RNON_BORUS']
sav_inds = [sav_inds,'IIWDET_BORUS','IIWNON_BORUS','IIRDET_BORUS','IIRNON_BORUS']


sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="nh_dist.sav"')


END




