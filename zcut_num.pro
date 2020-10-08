

pec = photoerrorclass
iizo = ztype ne 'ZP'
iizp = ztype eq 'ZP'
iipec3 = pec ge -1 and pec le 3
iipec1 = pec eq 1
iiz08 = z le 0.8
iiz06 = z le 0.6

print,''
print,'                 -------      -------      -------'
print,'                   ALL          DET        NON-DET'
print,'                 -------      -------      -------'
print,'QUAL:      ', total(iiqual), total(iiqual_det), total(iiqual_non)
print,'QUAL/=ZP:  ', total(iiqual and iizo), total(iiqual_det and iizo), total(iiqual_non and iizo)
print,'QUAL+ZP:   ', total(iiqual and iizp), total(iiqual_det and iizp), total(iiqual_non and iizp)
print,''
print,'-1<PEC<3:  ', total(iiqual and (iizo or (iizp and iipec3))), total(iiqual_det and (iizo or (iizp and iipec3))), total(iiqual_non and (iizo or(iizp and iipec3)))
print,'QUAL+Z<0.8:', total(iiqual and (iizo or (iizp and iipec3)) and iiz08), total(iiqual_det and (iizo or (iizp and iipec3)) and iiz08), total(iiqual_non and (iizo or (iizp and iipec3)) and iiz08)
print,'QUAL+Z<0.6:', total(iiqual and (iizo or (iizp and iipec3)) and iiz06), total(iiqual_det and (iizo or (iizp and iipec3)) and iiz06), total(iiqual_non and (iizo or (iizp and iipec3)) and iiz06)
print,''
print,'PEC=1:     ',total(iiqual and (iizo or (iizp and iipec1))), total(iiqual_det and (iizo or (iizp and iipec1))), total(iiqual_non and (iizo or(iizp and iipec1)))
print,'QUAL+Z<0.8:', total(iiqual and (iizo or (iizp and iipec1)) and iiz08), total(iiqual_det and (iizo or (iizp and iipec1)) and iiz08), total(iiqual_non and (iizo or (iizp and iipec1)) and iiz08)
print,'QUAL+Z<0.6:', total(iiqual and (iizo or (iizp and iipec1)) and iiz06), total(iiqual_det and (iizo or (iizp and iipec1)) and iiz06), total(iiqual_non and (iizo or (iizp and iipec1)) and iiz06)
print,''

END














