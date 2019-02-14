FUNCTION clean_source_nustar, in_data


data = in_data

;; remove Galactic Center
ig = where(data.field ne 'GalCen',/null)
data = data[ig]
;; exposure > 0 (Fexp = Sexp + Hexp)
ig = where(data.fexp gt 0.,/null)
data = data[ig]
;; detections and not upper limits
ig = where(data.fbf gt 0. and data.e_fbf gt 0.,/null)
data = data[ig]
;; S/N doesn't matter when only interested in detections
;ig = where(data.fbf/data.e_fbf gt 2.,/null)
;data = data[ig]

;; sort by exposure time
isort = sort(data.fexp)
data = data[isort]

return, data


END