FUNCTION create_nustar_master_obsid_list


dir = '/Users/ccarroll/Research/surveys/NuSTAR/'
files = dir+['COSMOS/nst-obsid-list.fits', $
             ;'ECDFS/', $
             'SSC/nst-ssc-table1.fits', $
             ;'SSC2/', $
             'UDS/nst-uds-table1.fits']

obsid = []
for i = 0,n_elements(files)-1 do begin
    r = mrdfits(files[i],1)
    ii = where(r.obsid gt 0,/null)
    obsid = [obsid,r[ii].obsid]
endfor
uid = obsid[uniq(obsid,sort(obsid))]

return, uid


END