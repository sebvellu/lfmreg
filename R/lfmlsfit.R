#' @export

lfmlsfit <- function(
	yvls, zvls, xvls, wndw, dzdx, krnl = "ba", band = "and",
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
	vvls <- helperkit::diffdtrnd(xvls, dtrn, step)
	lrvr <- tvlrvarstd(yvls, zvls, vvls, wndw, krnl, band, slct)
    blck <- lrvr$blck
    lrso <- lrvr$lrso
    clrv <- lrvr$clrv
    lrvr <- lrvr$lrvr
	#
	xidx <- ynum + 1:xnum
	yidx <- 1:ynum
	#
	bico <- array(NA_real_, c(xnum, ynum, blck$blck))
	for (indx in seq_len(blck$blck)) {
		temp <- lrvr$deltvar[xidx, xidx, indx] %*% lrso[, , indx]
		temp <- lrvr$deltvar[xidx, yidx, indx] - temp
		bico[, , indx] <- temp
	}
	bico <- rowMeans(bico, dims = 2)
	bico <- rowSums(dzdx, dims = 2) %*% bico
	#
	ymod <- yvls[-1, , drop = FALSE]
	yadj <- matrix(NA_real_, lgth - 1, ynum)
	for (indx in seq_len(lgth - 1)) {
		#bidx <- tvblockidx(indx, lgth, mval) # not indx - 1, since t = 2,...,T not t=1,...,T, indx = t + 1!!
		temp <- vvls[indx, , drop = FALSE] %*% lrso[, , blck$indx[indx]]
		yadj[indx, ] <- ymod[indx, , drop = FALSE] - temp
	}
	wvls <- zvls[-1, , drop = FALSE]
	#
    qrwv <- qr(wvls)
	rqrw <- qr.R(qrwv)
	#whiv <- qr_safe_solve(t(rqrw), diag(1, ncol(rqrw)))
	whiv <- helperkit::safe_forwardsolve(t(rqrw), diag(1, ncol(rqrw)))
	cent <- crossprod(whiv)
    #
    # lsfit
    yols <- yadj - (wvls %*% cent) %*% bico
    #cffs <- qr.coef(qrwv, yadj) - cent %*% bico # ==
    cffs <- qr.coef(qrwv, yols)
    rsds <- yadj - wvls %*% cffs
    fits <- ymod - rsds
    cffs <- t(cffs)
	#
    if (is.null(rmat) || is.null(rvec) || is.null(wght)) {
        ahlf <- NULL
        cfff <- c(t(cffs))
        cffr <- cffs
        fitr <- fits
        rsdr <- rsds
    } else {
	    rslt <- helperkit::rfmls(ymod, yadj, yols, xvls, rmat, rvec, wght)
	    ahlf <- rslt$ahlf
        cfff <- rslt$cfff
	    cffr <- rslt$cffr
	    fitr <- rslt$fitr
	    rsdr <- rslt$rsdr
    }
    #
    ufit <- tcrossprod(zvls, cffs)
    ursd <- yvls - ufit
    rfit <- tcrossprod(zvls, cffr)
    rrsd <- yvls - rfit
	#
    return(list(
        iobs = seq_len(lgth - 1) + 1,
        long = lrvr$longvar,
        clrv = clrv,
        lrso = lrso,
		blck = blck,
        #
        vvls = vvls,
        #
        ymod = ymod,
        yadj = yadj,
        yols = yols,
        wvls = wvls,
        qrwv = qrwv,
        #
        cffs = cffs,
        fits = fits,
        rsds = rsds,
        #
        ahlf = ahlf,
        cfff = cfff,
        cffr = cffr,
        fitr = fitr,
        rsdr = rsdr,
        #
        ufit = ufit,
        ursd = ursd,
        rfit = rfit,
        rrsd = rrsd
    ))
}
