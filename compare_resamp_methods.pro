dir = ['run_20200130_ebv50_flat/','run_20200129_ebv50/']

pushd,dir[0]

load_vars,'fits.sav','_fits'
fvars = scope_varname(common='_fits')
fvars_flat = fvars+'_FLAT'
restore,'fits.sav'
for i = 0,n_elements(fvars)-1 do re = execute(fvars_flat[i]+' = '+fvars[i])

load_vars,'quality_src.sav','_qual'
qvars = scope_varname(common='_qual')
qvars_flat = qvars+'_FLAT'
restore,'quality_src.sav'
for i = 0,n_elements(qvars)-1 do re = execute(qvars_flat[i]+' = '+qvars[i])

load_vars,'src_luminosity.sav','_lum'
lvars = scope_varname(common='_lum')
lvars_flat = lvars+'_FLAT'
restore,'src_luminosity.sav'
for i = 0,n_elements(lvars)-1 do re = execute(lvars_flat[i]+' = '+lvars[i])

popd

pushd,dir[1]
restore,'fits.sav'
restore,'quality_src.sav'
restore,'src_luminosity.sav'
popd


END











