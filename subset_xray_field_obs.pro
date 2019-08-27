;; takes in OBS files, this is before fitting
PRO subset_xray_field_obs, in_files

;; output SAV file string
sav_str = ((strsplit(in_files,'/',/extract,/regex)).toArray())[*,-1]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; CHANDRA MASTER
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Chandra Master Path
mast_cha_path = '/Users/ccarroll/Resecha_arch/surveys/Chandra/*master*.fits'
;; 3XMM-DR8 Serendip Catalog Per-Observation Source Table
cat_cha_path = '/Users/ccarroll/Resecha_arch/surveys/Chandra/observation-source-2.fits'
;; Chandra Master cha_archive
cha_arch = mrdfits(mast_cha_path,1)
;; Master Catalog is updated more frequently than CSC2! 
;; avoid spurious non-detections!
cat_cha = mrdfits(cat_cha_path,1)
;; use only OBSID that are in cat_chaalots
mast_cha_id = cha_arch.obsid
cat_cha_id = cat_cha[where(cat_cha.instrument eq 'ACIS',/null)].obsid
cat_cha_id = cat_cha_id[uniq(cat_cha_id,sort(cat_cha_id))]
match,mast_cha_id,cat_cha_id,imast_cha,icat_cha
iicha_arch = bytarr(n_elements(cha_arch))
iicha_arch[imast_cha] = 1
cha_arch = cha_arch[where(iicha_arch,/null)]
;; use only cha_archived sources
cha_arch = cha_arch[where(cha_arch.status eq 'cha_archIVED' or cha_arch.status eq 'OBSERVED',/null)]
cha_arch = cha_arch[where(cha_arch.detector eq 'ACIS-I',/null)]
;; ACIS-I FOV is 16'x16'
;; https://heasarc.gsfc.nasa.gov/docs/chandra/chandra.html
fov_cha = 16.*60./2.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; XMM-NEWTON MASTER
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; XMM Master Path
mast_xmm_path = '/Users/ccarroll/Resexmm_arch/surveys/XMM/*master*.fits'
;; 3XMM-DR8 Serendip xmm_catalog
xmm_cat_xmm_path = '/Users/ccarroll/Resexmm_arch/surveys/XMM/3XMM_DR8xmm_cat_v1.0.fits'
;; XMM Master xmm_archive
xmm_arch = mrdfits(mast_xmm_path,1)
;; Master xmm_catalog is updated more frequently than 3XMM-DR8! 
;; avoid spurious non-detections!
xmm_cat = mrdfits(xmm_cat_xmm_path,1)
;; use only OBSID that are in xmm_catalots
mast_xmm_id = xmm_arch.obsid
xmm_cat_xmm_id = xmm_cat.obs_id
xmm_cat_xmm_id = xmm_cat_xmm_id[uniq(xmm_cat_xmm_id,sort(xmm_cat_xmm_id))]
match,mast_xmm_id,xmm_cat_xmm_id,imast_xmm,icat_xmm
iixmm_arch = bytarr(n_elements(xmm_arch))
iixmm_arch[imast_xmm] = 1
xmm_arch = xmm_arch[where(iixmm_arch,/null)]
;; use only xmm_archived sources (possibly use )
xmm_arch = xmm_arch[where(xmm_arch.status eq 'xmm_archIVED' or xmm_arch.status eq 'OBSERVED',/null)]  ;; observed sources
xmm_arch = xmm_arch[where(xmm_arch.pn_time gt 0.,/null)]                                      ;; ensure PN observation
xmm_arch = xmm_arch[where(xmm_arch.duration gt 0.,/null)]                                     ;; sanity check
iimode = strmatch(xmm_arch.pn_mode,'*FLG*',/fold) or $                          ;; ensure Large-Window or Full-Frame mode
         strmatch(xmm_arch.pn_mode,'*FF*',/fold) or $
         strmatch(xmm_arch.pn_mode,'*EFF*',/fold)
xmm_arch = xmm_arch[where(iimode,/null)]
;; XMM PN MOS FOV is ~27.5'x27.5'; use FOV inscribed circle--being conservative
;; https://heasarc.gsfc.nasa.gov/docs/xmm/xmm.html
fov_xmm = 27.5*60./2.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; NuSTAR MASTER
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Combined NuSTAR Fields Path
mast_path = '/Users/ccarroll/Research/surveys/NuSTAR/*master*.fits'
;; NuSTAR Catalogs
cat_path = '/Users/ccarroll/Research/surveys/NuSTAR/combined_nustar_fields.fits'
;; Read in the NuSTAR observation information (HEASARC);
arch = mrdfits(mast_path,1)
;; Master Catalog is updated more frequently than 3XMM-DR8! 
;; avoid spurious non-detections!
cat = mrdfits(cat_path,1)
;; NuSTAR FOV is 13'x13'
;; https://heasarc.gsfc.nasa.gov/docs/nustar/nustar.html
fov_nst = 13.*60./2.
spherematch,arch.ra,arch.dec,cat.ra,cat.dec,fov_nst/3600.,imast,icat,sep,maxmatch=0
iiarch = bytarr(n_elements(arch))
iiarch[imast] = 1
arch = arch[where(iiarch,/null)]
;; for NuSTAR, select just the science subset
arch = arch[where(arch.observation_mode eq 'SCIENCE',/null)]
;; pull data from ARCH
obsid = arch.obsid
ra_nst = arch.ra
dec_nst = arch.dec
rot_angle = arch.roll_angle
ontime = arch.ontime_a

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; LOOP OVER EACH FILE
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
for i = 0,n_elements(in_files)-1 do begin
    print, 'WORKING FIELD: '+in_files[i]
    
    restore,in_files[i]
    nsrc = n_elements(obs)

    iiinf_cha = bytarr(nsrc)
    iiinf_xmm = bytarr(nsrc)
    iiinf_nst = bytarr(nsrc)

    ;; CHANDRA 
    spherematch,ra,dec,cha_arch.ra,cha_arch.dec,fov_cha/3600.,isamp_cha,ifield,sep_cntr,maxmatch=0
    iiinf_cha[isamp_cha] = 1

    ;; XMM
    spherematch,ra,dec,xmm_arch.ra,xmm_arch.dec,fov_xmm/3600.,isamp_xmm,ifield,sep_cntr,maxmatch=0
    iiinf_xmm[isamp_xmm] = 1

    ;; NUSTAR
    spherematch,ra,dec,ra_nst,dec_nst,fov_nst/3600.,is,ix,sepnu,maxmatch=0
    isu = is[uniq(is,sort(is))]
    ixu = ix[uniq(ix,sort(ix))]
    ;; for each sample source
    for n=0L,n_elements(ra[isu])-1 do begin 
       ;; for each x-ray observation
       for m=0L,n_elements(ra_nst[ixu])-1 do begin 
          GCIRC, 2, ra_nst[ixu[m]],dec_nst[ixu[m]],ra[isu[n]],dec[isu[n]],dist_test
          if (dist_test le 2000.) then begin
             nustar_fov,ra_nst[ixu[m]],dec_nst[ixu[m]],rot_angle[ixu[m]],box_enc_x,box_enc_y   
             dummy=IsPointInPolygon(box_enc_x,box_enc_y,ra[isu[n]],dec[isu[n]])
             if (dummy eq -1) then iiinf_nst[isu[n]] = 1
          endif
       endfor
    endfor

    ;; COMBINE ALL FIELDS NOW
    iiinf = iiinf_cha or iiinf_xmm or iiinf_nst
    iinf = where(iiinf,ct)
    if (ct eq 0) then continue
    obs = obs[iinf]
    save,obs,band,/compress,file=sav_str[i]    
endfor


END








