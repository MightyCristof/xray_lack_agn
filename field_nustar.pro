PRO field_nustar


common _fits


;; Combined NuSTAR Fields
pth1 = '/Users/ccarroll/Research/surveys/NuSTAR/nustar_fields.fits'

;; OBSERVATIONS
;; *******************************
;; taken from FIND_SRC_IN_NSTR.PRO
;; *******************************
    nsrc = n_elements(ra) ;; How many srcs
    ;;;; ========= X-ray Database Information ==============
    ;; Read in the NuSTAR observation information (HEASARC);
    pth1 = '/Users/ccarroll/Research/surveys/NuSTAR/'
    nm = file_search(pth1+'*master*.fits')
    arch = mrdfits(nm,1)
    obsid = arch.obsid
    xra = arch.ra
    xdec = arch.dec
    rot_angle = arch.roll_angle
    instrument = arch.observation_mode
    ontime = arch.ontime_a
    ;; We only require ACIS imaging pointings so...(remove the ones that aren't ACIS-I; rather than sub-select; 
    ;; for NuSTAR, select just the science subset
    scien=where(instrument eq 'SCIENCE',complement=ibad)  ;  SCIENCE instead of 'ACIS-I'
    remove,ibad,obsid,xra,xdec,rot_angle,instrument,ontime  ; , grating,seq_num,
    nnstr = n_elements(scien)  ; go from 4222 to 2065
    ;; ==========================================================
    ;; Now we need to find which sources are within the FOV of a
    ;; NuSTAR observation. 
    ;; Initialize the output arrays
    src_infield = lonarr(nsrc)
    ;src_nnstr_obs = intarr(nsrc)
    ;src_nstr_obsid = strarr(nsrc)
    src_nstr_expt  = strarr(nsrc)
    src_nstr_dist = strarr(nsrc)
    ;; initial rough cut on NuSTAR field centers using spherematch, 
    ;; not necessary to do one object at a time
    nufov = 13.
    hypot = sqrt(2.*nufov^2)
    spherematch,ra,dec,xra,xdec,hypot/60.,is,ix,sepnu,maxmatch=0
    isu = is[uniq(is,sort(is))]
    ixu = ix[uniq(ix,sort(ix))]
    ;; Step through each test source
    for n=0L,n_elements(ra[isu])-1 do begin 
       ;; For each test source, check if it's in the FOV of the X-ray observation
       for m=0L,n_elements(xra[ixu])-1 do begin 
          ;; Quickly work out the distance between source and NuSTAR pointing
          ;; if greater than 0.75 degrees then this doesn't concern us (need to double check the max dist. bt pointing and edge of NuSTAR field)
          ;; GCIRC is idlastro for rigourus great circle distance, returns in arcsec
          GCIRC, 2, xra[ixu[m]],xdec[ixu[m]],ra[isu[n]],dec[isu[n]],dist_test
          if dist_test le 2000. then begin  ; acis-i is 16.9'x16.9'; nustar is 13'x13', assuming similar pointing center 2000/2700 ~ 13/16.9
             ;; Calculate Ra,Dec of the FOV corners of Chandra ACIS-I for
             ;; a given pointing ra,dec and roll angle
             ;; this will need to be changed to NuSTAR fov (also is fpma/fpmb same fov?)
             nustar_fov,xra[ixu[m]],xdec[ixu[m]],rot_angle[ixu[m]],box_enc_x,box_enc_y   
             ;; Is Ra,Dec for src inside the box enclosed by NuSTAR's FOV?
             ;; this part should not require any changes
             dummy=IsPointInPolygon(box_enc_x,box_enc_y,ra[isu[n]],dec[isu[n]])
             if dummy eq -1 then begin ;; Src in Chandra FOV
                ;; Keep track of Src with X-ray observation
                ;; strip removes leading and trailing blank spaces
                src_infield[isu[n]] = 1
                src_nstr_expt[isu[n]] = src_nstr_expt[isu[n]]+strip(ontime[ixu[m]])+','
                src_nstr_dist[isu[n]] = src_nstr_dist[isu[n]]+strip(dist_test)+','
             endif
          endif
       endfor
       if n mod 1000 eq 0 then print, 'At Object ', n, ' of ', n_elements(ra[isu])
    endfor
;; *******************************
;; 
;; *******************************

;; prep NuSTAR fields (multiple pointings)
iiinf_nst_field = byte(src_infield)
texp_nst_field = dblarr(nsrc)
sdst_nst_field = dblarr(nsrc)
inst = where(iiinf_nst_field,nnst)
for i = 0,nnst-1 do begin
    temp_time = double(strsplit(src_nstr_expt[inst[i]],',',/extract))
    temp_dist = double(strsplit(src_nstr_dist[inst[i]],',',/extract))
    ;; take the observation where source is closest to pointing center
    sdst_nst_field[inst[i]] = min(temp_dist,imin)
    texp_nst_field[inst[i]] = temp_time[imin]
    ;; take the observation where source has longest exposure time
    ;expt_nst[inus[i]] = max(temp_time,imax)
    ;sdst_nst[inus[i]] = temp_dist[imax]
endfor

save,iiinf_nst_field,texp_nst_field,sdst_nst_field, $
     /compress,filename='xfield_nst.sav' 


END





