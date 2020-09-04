FUNCTION xray_lack_sky_area, phi, theta, rad, IND=ind


if keyword_set(ind) then ind = ind else ind = lindgen(n_elements(phi))

nside = 256;512;256
npix = nside2npix(nside)
sky_area = 41252.96
pix_area = sky_area/npix

ang2vec,theta[ind],phi[ind],vec,/astro
map = dblarr(npix,3)

for i = 0,n_elements(vec[*,0])-1 do begin

    query_disc,nside,vec[i,*],rad/3600.,listpix,nlist,/deg
    if (nlist gt 0) then begin
        map[listpix,0] = pix_area
        pix2vec_ring,nside,listpix,vec0
        vec2ang,vec0,theta_out,phi_out,/astro
        map[listpix,1]=phi_out
        map[listpix,2]=theta_out
    endif
    
endfor

return, map


END




;IDL> print, total(map)
;       263.06611
;
;IDL> print, 41252.96/total(map)
;       156.81595