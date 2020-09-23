FUNCTION stack_det_flux,  logfx


yfx_log = histogram(logfx,locations=xfx_log,bin=scott(logfx))
yfx_lin = histogram(10.^logfx,locations=xfx_lin,bin=scott(10.^logfx))

fx_stack = {fxlog:mean(logfx),ylog:yfx_log,xlog:xfx_log, $
            fxlin:mean(10.^logfx),ylin:yfx_lin,xlin:xfx_lin}

return, fx_stack


END




