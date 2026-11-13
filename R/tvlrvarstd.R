tvlrvarstd <- function(
	yvls, zvls, vvls, wndw, krnl = "ba", band = "band", slct = 1
) {
	yvls <- as.matrix(yvls)
	#
	rols <- cbind(qr.resid(qr(zvls), yvls)[-1, ], vvls)
	#rols <- yvls - zvls %*% qr_safe_solve(zvls, yvls)
	#rols <- cbind(rols[-1, ], vvls)
	blck <- blocks(wndw, nrow(rols), slct)
    lrvr <- lrvar::tvlrvar(rols, blck$strt, blck$fnsh, krnl, band)
	long <- lrvr$longvar
	#
	xidx <- ncol(yvls) + seq_len(ncol(vvls))
	yidx <- seq_len(ncol(yvls))
	#
	lrso <- array(NA_real_, c(ncol(vvls), ncol(yvls), blck$blck))
	clrv <- array(NA_real_, c(ncol(yvls), ncol(yvls), blck$blck))
	#
	for (indx in seq_len(blck$blck)) {
		lrso[, , indx] <- helperkit::safesymsol(
			mata = lrvr$longvar[xidx, xidx, indx], 
			matb = lrvr$longvar[xidx, yidx, indx]
		)
		temp <- lrvr$longvar[yidx, xidx, indx] %*% lrso[, , indx]
		clrv[, , indx] <- lrvr$longvar[yidx, yidx, indx] - temp
	}
    #
    return(list(
        rols = rols,
        lrvr = lrvr,
        lrso = lrso,
        clrv = clrv,
		blck = blck
    ))
}
