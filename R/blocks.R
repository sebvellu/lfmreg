blocks <- function(wndw, lgth, slct = 1) { 
	tval <- seq_len(lgth)
	if (slct == 1) {
		# Non-overlapping windows
		blck <- lgth %/% wndw
		fnsh <- seq_len(blck) * wndw
		strt <- fnsh - wndw + 1
		# Extend final block estimate to remaining observations
		indx <- pmin((tval - 1) %/% wndw + 1, blck)
	} else if (slct == 2) {
		# Overlapping windows (non-truncated)
		blck <- lgth - wndw + 1
		strt <- seq_len(blck)
		fnsh <- strt + wndw - 1
		indx <- pmin(pmax(tval - (wndw - 1) %/% 2, 1), blck)
	} else {
		# Overlapping windows (truncated)
		blck <- lgth
		strt <- pmax(tval - (wndw - 1) %/% 2, 1)
		fnsh <- pmin(tval + wndw %/% 2, blck)
		indx <- tval
	}
	return(list(blck = blck, wndw = wndw, strt = strt, fnsh = fnsh, indx = indx))
}
