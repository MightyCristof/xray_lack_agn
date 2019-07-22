FUNCTION interp_6micr_lum, nufnu_in, $
                           z_in, $
                           LOG = log


nufnu = nufnu_in
z = z_in

dl2 = dlum(z)^2         ;; luminosity distance in cm^2

lum_out = 4.*!const.pi*dl2*nufnu
if keyword_set(log) then lum_out = alog10(lum_out) > 0.

return, lum_out


END




