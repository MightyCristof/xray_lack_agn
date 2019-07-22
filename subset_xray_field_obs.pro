PRO subset_xray_field_obs, in_files


;; output SAV file string
sav_str = ((strsplit(in_files,'/',/extract,/regex)).toArray())[*,-1]


;; X-RAY FIELD MASTER ARCHIVES

;; Chandra Master
cha_pth = '/Users/ccarroll/Research/surveys/Chandra/*master*.fits'
;; OBSERVATIONS
;; use only archived sources
cha_arch = mrdfits(cha_pth,1)
cha_arch = cha_arch[where(cha_arch.status eq 'ARCHIVED' or cha_arch.status eq 'OBSERVED',/null)]
cha_arch = cha_arch[where(cha_arch.detector eq 'ACIS-I',/null)]
;; ACIS-I FOV is 16'x16'; use inscribed circle--being conservative
;; https://heasarc.gsfc.nasa.gov/docs/chandra/chandra.html
fov_cha = 16.*60./2.


;; XMM Master
xmm_pth = '/Users/ccarroll/Research/surveys/XMM/*master*.fits'
;; OBSERVATIONS
;; use only archived sources (possibly use )
xmm_arch = mrdfits(xmm_pth,1)
xmm_arch = xmm_arch[where(xmm_arch.status eq 'ARCHIVED' or xmm_arch.status eq 'OBSERVED',/null)]  ;; observed sources
xmm_arch = xmm_arch[where(xmm_arch.pn_time gt 0.,/null)]                                      ;; ensure PN observation
xmm_arch = xmm_arch[where(xmm_arch.duration gt 0.,/null)]                                     ;; sanity check
iimode = strmatch(xmm_arch.pn_mode,'*FLG*',/fold) or $                          ;; ensure Large-Window or Full-Frame mode
         strmatch(xmm_arch.pn_mode,'*FF*',/fold) or $
         strmatch(xmm_arch.pn_mode,'*EFF*',/fold)
xmm_arch = xmm_arch[where(iimode,/null)]
;; XMM MOS FOV is ~33'x33'; use inscribed circle--being conservative
;; https://heasarc.gsfc.nasa.gov/docs/xmm/xmm.html
fov_xmm = 27.5*60./2.


;; NuSTAR Master
nst_pth = '/Users/ccarroll/Research/surveys/NuSTAR/*master*.fits'
;; OBSERVATIONS
;; use only ACIS-I SCIENCE
nst_arch = mrdfits(nst_pth,1)
nst_arch = nst_arch[where(nst_arch.observation_mode eq 'SCIENCE',/null)]
xra = nst_arch.ra
xdec = nst_arch.dec
rot_angle = nst_arch.roll_angle
ontime = nst_arch.ontime_a
;; NuSTAR FOV; FOV inscribed circle--being conservative
;; https://heasarc.gsfc.nasa.gov/docs/nustar/nustar.html
fov_nst = 13.*60./2.
hypot = sqrt(2*fov_nst^2.)
fov_nst = hypot

for i = 0,n_elements(in_files)-1 do begin
    print, 'WORKING FIELD: '+in_files[i]
    
    restore,in_files[i]
    nsrc = n_elements(obs)
    
    ;; Chandra
    iiinf_cha = bytarr(nsrc)
    texp_cha = dblarr(nsrc)
    sdst_cha = dblarr(nsrc)
    spherematch,obs.ra,obs.dec,cha_arch.ra,cha_arch.dec,fov_cha/3600.,isamp,icha,sep_cha,maxmatch=0
    ;; tag main sample sources as "in field"
    if (icha[0] ne -1) then begin
        iiinf_cha[isamp] = 1
        texp_cha[isamp] = cha_arch[icha].exposure
        sdst_cha[isamp] = sep_cha
    endif
    
    
    ;; XMM 
    iiinf_xmm = bytarr(nsrc)
    texp_xmm = dblarr(nsrc)
    sdst_xmm = dblarr(nsrc)
    spherematch,obs.ra,obs.dec,xmm_arch.ra,xmm_arch.dec,fov_xmm/3600.,isamp,ixmm,sep_xmm,maxmatch=0
    ;; tag main sample sources as "in field"
    if (ixmm[0] ne -1) then begin
        iiinf_xmm[isamp] = 1
        texp_xmm[isamp] = xmm_arch[ixmm].duration
        sdst_xmm[isamp] = sep_xmm
    endif
    
    
    ;; NuSTAR
    iiinf_nst = bytarr(nsrc)
    texp_nst = dblarr(nsrc)
    sdst_nst = dblarr(nsrc)
    
    src_infield = lonarr(nsrc)
    src_nstr_expt  = strarr(nsrc)
    src_nstr_dist = strarr(nsrc)
    spherematch,obs.ra,obs.dec,xra,xdec,fov_nst/3600.,is,ix,sepnu,maxmatch=0
    if (is[0] ne -1) then begin
        isu = is[uniq(is,sort(is))]
        ixu = ix[uniq(ix,sort(ix))]
        for n=0L,n_elements(obs[isu].ra)-1 do begin 
           ;; For each test source, check if it's in the FOV of the X-ray observation
           for m=0L,n_elements(xra[ixu])-1 do begin 
              ;; Quickly work out the distance between source and NuSTAR pointing
              ;; if greater than 0.75 degrees then this doesn't concern us (need to double check the max dist. bt pointing and edge of NuSTAR field)
              ;; GCIRC is idlastro for rigourus great circle distance, returns in arcsec
              GCIRC, 2, xra[ixu[m]],xdec[ixu[m]],obs[isu[n]].ra,obs[isu[n]].dec,dist_test
              if dist_test le 2000. then begin  ; acis-i is 16.9'x16.9'; nustar is 13'x13', assuming similar pointing center 2000/2700 ~ 13/16.9
                 ;; Calculate Ra,Dec of the FOV corners of Chandra ACIS-I for
                 ;; a given pointing ra,dec and roll angle
                 ;; this will need to be changed to NuSTAR fov (also is fpma/fpmb same fov?)
                 nustar_fov,xra[ixu[m]],xdec[ixu[m]],rot_angle[ixu[m]],box_enc_x,box_enc_y   
                 ;; Is Ra,Dec for src inside the box enclosed by NuSTAR's FOV?
                 ;; this part should not require any changes
                 dummy=IsPointInPolygon(box_enc_x,box_enc_y,obs[isu[n]].ra,obs[isu[n]].dec)
                 if dummy eq -1 then begin ;; Src in Chandra FOV
                    ;; Keep track of Src with X-ray observation
                    ;; strip removes leading and trailing blank spaces
                    src_infield[isu[n]] = 1
                    src_nstr_expt[isu[n]] = src_nstr_expt[isu[n]]+strip(ontime[ixu[m]])+','
                    src_nstr_dist[isu[n]] = src_nstr_dist[isu[n]]+strip(dist_test)+','
                 endif
              endif
           endfor
           ;if n mod 1000 eq 0 then print, 'At Object ', n, ' of ', n_elements(obs[isu].ra)
        endfor
        isrc = where(src_infield,src_ct)
        if (src_ct gt 0) then iiinf_nst[where(src_infield)] = 1
    endif
    
    
    ;; COMBINE ALL FIELDS NOW
    iiinf = iiinf_cha or iiinf_xmm or iiinf_nst
    iinf = where(iiinf,ct)
    if (ct eq 0) then continue
    obs = obs[iinf]
    save,obs,band,/compress,file=sav_str[i]    
endfor


END








