FUNCTION mc_nondet_dist, logl6um, $
                         wbinc, $
                         wdc, $
                         det, $
                         rl, $
                         ;pick, $
                         SEED = seed, $
                         PLT = plt


pick = 2
;; number of sources == number of random pulls
nobj = n_elements(logl6um)

;; read differential KM
;readcol,'wagn1_10x05.out',wbinc,wdc,format='d,d',/silent,skipline=587

;; CHOICE
;; bin start vs. bin center?
case pick of
    0:  begin
        ;; BINS CENTERED; NO END BINS
        wbin = wbinc
        wd = wdc
        end
        ;; ends work--abrupt cutoff at -3
        ;; more detections than total in certain bins
    1:  begin
        ;; BINS CENTERED; ADD END BINS
        wbin = [wbinc[0]-0.5,wbinc,wbinc[-1]+0.5]
        wd = [0.,wdc,0.]
        ;; ends don't work--goes passed -3
        ;; still more detections than total in certain bins
        end
    2:  begin
        ;; BIN START
        wbin = wbinc-0.25
        wd = wdc
        ;; ends don't work--goes passed -3
        ;; fixes detections in excess of total in certain bins
        end
    3:  begin
        ;; BIN START; ADD END BINS
        wbin = [wbinc[0]-0.75,wbinc-0.25,wbinc[-1]+0.25]
        wd = [0.,wdc,0.]
        end
    else:   begin
            stop
            end
endcase

;; match input data to same bins
;readcol,'wagn1.dat',det,rl,format='d,d'
yrl = histogram(rl,locations=xrl,binsize=0.5,min=wbin[0],max=wbin[-1])
yrld = histogram(rl[where(det eq 0)],locations=xrld,binsize=0.5,min=wbin[0],max=wbin[-1])
yrln = histogram(rl[where(det eq -1)],locations=xrln,binsize=0.5,min=wbin[0],max=wbin[-1])
;; subtract detections from underlying distribution
if (n_elements(wd) ne n_elements(yrld)) then stop
ynon = wd-yrld > 0.
xnon = wbin
if (pick eq 2 or pick eq 3) then begin
    xrld += 0.25
    xrln += 0.25
    xnon += 0.25
endif
;; create finer grid
dx = [xnon[0]:xnon[-1]:diff(minmax(xnon))/999.]
;; interpolate grid to create PDF
pdf = spline(xnon,ynon,dx) > 0.
pdf /= total(pdf)
;; sum PDF to get CDF
cdf = total(pdf,/cumulative)
;; unique elements of CDF
iu = uniq(cdf)
udx = dx[iu]
cdf = cdf[iu]
;; set seed
if keyword_set(seed) then seed = seed
;; draw random values from 0-1 from CDF
draw_cdf = randomu(seed,nobj)
;; project CDF back to PDF
;; remember newY=interpol(Y,X,newX)
samp_pdf = interpol(udx,cdf,draw_cdf)
;; histogram sampled PDF
ypdf = histogram(samp_pdf,locations=xpdf,bin=0.5,min=xnon[0],max=xnon[-1])

if keyword_set(plt) then begin
    p = plot(wbinc,wdc,/stairstep,xra=[-3.5,2.])
    p = plot(xrld,yrld,/stairstep,/ov,col='dodger blue')
    p = plot(xnon,ynon,/stairstep,/ov,col='orange')
    p = plot(xpdf,nm(ypdf)*max(wdc),/stairstep,/ov,col='red')
    ;; plot fine grid
    ;samp_plt = interpol(udx,cdf,randomu(seed,100000))
    ;yplt = histogram(samp_plt,locations=xplt,bin=width(dx),min=dx[0],max=dx[-1])
    ;p = plot(xplt,nm(yplt)*max(wdc),/stairstep,/ov,col='red',transparency=75)

    ;p = plot(xpdf,ypdf,/stairstep)
    ;plothist,samp_pdf,bin=0.5,xhist,yhist
    ;p = plot(xhist,yhist,/stairstep,/ov,col='red')
endif

return, samp_pdf

END








