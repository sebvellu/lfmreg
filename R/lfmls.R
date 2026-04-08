#' @export

lfmls <- function(
	yvls, zvls, xvls, sgmh, dzdx, krnl = "ba", band = "and",
    rmat = NULL, rvec = NULL, wght = NULL, 
	dtrn = NULL, step = 2, slct = 1
) {
	yvls <- as.matrix(yvls)
	zvls <- as.matrix(zvls)
	xvls <- as.matrix(xvls)
	#
	lgth <- nrow(yvls)
	ynum <- ncol(yvls)
	znum <- ncol(zvls)
	xnum <- ncol(xvls)
	#
    rslt <- lfmlsfit(
		yvls = yvls, 
		zvls = zvls,
		xvls = xvls,
		sgmh = sgmh,
		dzdx = dzdx,
        krnl = krnl,
        band = band,
		rmat = NULL, 
		rvec = NULL,
		wght = NULL,
		dtrn = dtrn,
		step = step,
		slct = slct
	)
	vvls <- rslt$vvls
	long <- rslt$long
	clrv <- rslt$clrv
    ymod <- rslt$ymod
	yadj <- rslt$yadj
    yols <- rslt$yols
	wvls <- rslt$wvls
	qrwv <- rslt$qrwv
	cffs <- rslt$cffs
	fits <- rslt$fits
	rsds <- rslt$rsds
	ufit <- rslt$ufit
	ursd <- rslt$ursd
	bidx <- rslt$bidx
	#
	cvls <- wvls
	wnum <- ncol(wvls)
	#
	cent <- crossprod(helperkit::sqrtinvprod(qrwv, FALSE)) 
	#
	bred <- diag(1, ynum) %x% cent
	meat <- 0
	for (indx in seq_len(nrow(wvls))) {
		#bidx <- tvblockidx(indx, lgth, mval)
		temp <- helperkit::array3tomat(clrv, bidx$lidx[indx])
		meat <- meat + (temp %x% tcrossprod(wvls[indx, ]))
	}
	cvrs <- bred %*% (meat %*% bred)
	csds <- matrix(sqrt(diag(cvrs)), wnum, ynum)
	csds <- t(csds)
	#
	# Restricted Case:
	#
	if (is.null(rmat)) {
		rmat <- diag(1, wnum * ynum)
	}
	if (is.null(rvec)) {
		rvec <- numeric(wnum * ynum)
	}
	if (is.null(wght)) {
		wght <- diag(1, ynum)
	}
	#
	rslt <- helperkit::rfmls(ymod, yadj, yols, wvls, rmat, rvec, wght)
	ahlf <- rslt$ahlf
	cfff <- rslt$cfff
	cffr <- rslt$cffr
	fitr <- rslt$fitr
	rsdr <- rslt$rsdr
	#
	rfit <- tcrossprod(zvls, cffr)
	rrsd <- yvls - rfit
	#
    ahiv <- helperkit::sqrtinvprod(qr(ahlf), TRUE)
	#
	mult <- wght %x% diag(1, wnum)
	meat <- mult %*% (meat %*% mult)
	meat <- crossprod(rmat, tcrossprod(meat, rmat))
	#
	cvrf <- ahiv %*% (meat %*% ahiv)
	cvrr <- tcrossprod(rmat %*% cvrf, rmat)
	#
	csdr <- matrix(sqrt(diag(cvrr)), wnum, ynum)
	csdr <- t(csdr)
	#
	# Not yet implemented
	#
	return(list(
		cffs = cffs, #cffu
		fits = fits, #fitu
		rsds = rsds, #rsdu
		cvrs = cvrs, #cvru
		csds = csds, #csdu
		clrv = clrv,
		cffr = cffr,
		fitr = fitr,
		rsdr = rsdr,
		cvrr = cvrr, 
		csdr = csdr,
		cfff = cfff,
		cvrf = cvrf,
		ufit = ufit,
		ursd = ursd,
		rfit = rfit,
		rrsd = rrsd
	))
}
