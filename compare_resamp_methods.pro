dir = ['run_20200129_ebv50/','run_20200130_ebv50_flat/']

pushd,dir[1]
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
popd

pushd,dir[0]
restore,'fits.sav'
restore,'quality_src.sav'
popd


