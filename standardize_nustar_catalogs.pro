FUNCTION standardize_nustar_catalogs, in_rr, $
                                      in_field


field = strupcase(in_field)
case field of 
    'COSMOS': begin
        ;; use NuSTAR positions for sources without soft X-ray counterpart; Section 5.4 of source paper
        rr = in_rr
        snum = rr.source_number
        match,snum,[135,141,245],irr,inoc
        rr[irr].cxoxmm_ra = rr[irr].ra
        rr[irr].cxoxmm_dec = rr[irr].dec
        ;; flag upper limits with flux_error == -99
        rr[where(rr.sb_flux_limit eq '<')].sb_flux_error = -99.
        rr[where(rr.hb_flux_limit eq '<')].hb_flux_error = -99.
        rr[where(rr.bb_flux_limit eq '<')].bb_flux_error = -99.
    end
    'ECDFS': begin
        stop
        ;; remove duplicate NuSTAR sources with dual soft X-ray counterparts; Section 2.3.1 of source paper
        rr = in_rr[sort(in_rr.source_number)]
        snum = rr.source_number
        w = width(snum,all=all)
        idup = where(~all)
        offset = rr.offset
        ;; keep source with smallest counterpart offset
        for i = 0,n_elements(idup)-1 do if (offset[idup[i]] lt offset[idup[i]+1]) then idup[i]++
        iu = exclude(rr,idup)
        rr = rr[iu]
        ;; use NuSTAR positions for sources without soft X-ray counterpart; Section 4 of source paper
        match,rr.source_number,[37],irr,inoc
        rr[irr].ctrprt_ra = rr[irr].ra
        rr[irr].ctrprt_dec = rr[irr].dec
        ;; net exposure time A+B
        rr.sb_exposure *= 0.5
        rr.hb_exposure *= 0.5
        rr.fb_exposure *= 0.5
    end
    'SSC': begin
        rr = in_rr
        ;; use optical counterpart positions for sources without soft X-ray counterpart; Section 3.2 of source paper
        irr = where(rr.softxray_ctrpart_ra eq 0. and rr.opt_ra gt 0.)
        rr[irr].softxray_ctrpart_ra = rr[irr].opt_ra
        rr[irr].softxray_ctrpart_dec = rr[irr].opt_dec
        ;; net exposure time A+B
        rr.sb_exposure *= 0.5
        rr.hb_exposure *= 0.5
        rr.fb_exposure *= 0.5
    end
    'SSC2': begin
        rr = in_rr
        ;; A.6 of source paper: 4 sources are left blank as in these cases the A+B data prohibit reliable photometric constraints
        rr = rr[where(rr.sb_counts and rr.hb_counts and rr.fb_counts)]
        ;; use optical counterpart positions for sources without soft X-ray counterpart; Section 3.2 of source paper
        irr = where(rr.softxray_ctrpart_ra eq 0. and rr.opt_ra gt 0.)
        rr[irr].softxray_ctrpart_ra = rr[irr].opt_ra
        rr[irr].softxray_ctrpart_dec = rr[irr].opt_dec
        ;; net exposure time A+B
        rr.sb_exposure *= 0.5
        rr.hb_exposure *= 0.5
        rr.fb_exposure *= 0.5
    end
    'UDS': begin
        rr = in_rr
        ;; counterpart RA/Dec
        irr = where(rr.xmmra eq -99.)
        rr[irr].xmmra = rr[irr].radeg
        rr[irr].xmmde = rr[irr].dedeg
        ;; net exposure time A+B, exposure time in ks
        rr._3_8exp  *= 500.
        rr._8_24exp *= 500.
        rr._3_24exp *= 500.        
    end
endcase

return, rr


END



