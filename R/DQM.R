##******************************************************************************
# Detrended quantile matching algorithm for direct climate downscaling
# Alex Cannon (acannon@uvic.ca)
##******************************************************************************

DQM <-
# Detrended quantile matching with delta-method extrapolation
function(o.h, g.h, g.f, ratio=TRUE, detrend=TRUE, n.max=NULL)
{
    if(ratio){
        o.h[o.h < sqrt(.Machine$double.eps)] <-
            runif(sum(o.h < sqrt(.Machine$double.eps)), min=0,
                  max=sqrt(.Machine$double.eps))
        g.h[g.h < sqrt(.Machine$double.eps)] <-
            runif(sum(g.h < sqrt(.Machine$double.eps)), min=0,
                  max=sqrt(.Machine$double.eps))
        g.f[g.f < sqrt(.Machine$double.eps)] <-
            runif(sum(g.f < sqrt(.Machine$double.eps)), min=0,
                  max=sqrt(.Machine$double.eps))
    }
    o.h.mn <- mean(o.h)
    g.h.mn <- mean(g.h)
    if(ratio){
        g.f <- g.f/g.h.mn*o.h.mn
    } else{
        g.f <- g.f-g.h.mn+o.h.mn
    }
    if(detrend){
        g.f.mn <- lm.fit(cbind(1, seq_along(g.f)), g.f)$fitted
    } else{
        g.f.mn <- o.h.mn
    }
    if(is.null(n.max)) n.max <- max(length(o.h), length(g.h))
    tau <- c(0, (1:n.max)/(n.max+1), 1)
    if(ratio & any(o.h < sqrt(.Machine$double.eps))){
        x <- quantile(g.h/g.h.mn, tau)
        y <- quantile(o.h/o.h.mn, tau)
        yout <- approx(x, y, xout=g.f/g.f.mn, rule=2:1)$y
        extrap <- is.na(yout)
        yout[extrap] <- max(o.h/o.h.mn)*((g.f/g.f.mn)[extrap]/max(g.h/g.h.mn))
        yout <- yout*g.f.mn
        yout.h <- approx(x, y, xout=g.h/g.h.mn, rule=1)$y*o.h.mn
    } else if(ratio & !any(o.h < sqrt(.Machine$double.eps))){
        x <- quantile(g.h/g.h.mn, tau)
        y <- quantile(o.h/o.h.mn, tau)
        yout <- approx(x, y, xout=g.f/g.f.mn, rule=1)$y
        extrap.lower <- is.na(yout) & ((g.f/g.f.mn) < min(g.h/g.h.mn))
        extrap.upper <- is.na(yout) & ((g.f/g.f.mn) > max(g.h/g.h.mn))
        yout[extrap.lower] <- min(o.h/o.h.mn)*((g.f/g.f.mn)[extrap.lower]/
                              min(g.h/g.h.mn))
        yout[extrap.upper] <- max(o.h/o.h.mn)*((g.f/g.f.mn)[extrap.upper]/
                              max(g.h/g.h.mn))
        yout <- yout*g.f.mn
        yout.h <- approx(x, y, xout=g.h/g.h.mn, rule=1)$y*o.h.mn
    } else{
        x <- quantile(g.h-g.h.mn, tau)
        y <- quantile(o.h-o.h.mn, tau)
        yout <- approx(x, y, xout=g.f-g.f.mn, rule=1)$y
        extrap.lower <- is.na(yout) & ((g.f-g.f.mn) < min(g.h-g.h.mn))
        extrap.upper <- is.na(yout) & ((g.f-g.f.mn) > max(g.h-g.h.mn))
        yout[extrap.lower] <- min(o.h-o.h.mn) + ((g.f-g.f.mn)[extrap.lower]-
                              min(g.h-g.h.mn))
        yout[extrap.upper] <- max(o.h-o.h.mn) + ((g.f-g.f.mn)[extrap.upper]-
                              max(g.h-g.h.mn))
        yout <- yout+g.f.mn
        yout.h <- approx(x, y, xout=g.h-g.h.mn, rule=1)$y+o.h.mn
    }
    if(ratio){
        yout[yout < sqrt(.Machine$double.eps)] <- 0
        yout.h[yout.h < sqrt(.Machine$double.eps)] <- 0
    }
    list(g.h.bc=yout.h, g.f.bc=yout)
}

mnDQM <-
# Monthly detrended quantile matching with delta-method extrapolation
# h = historical period (bias correction calibration)
# f = future period (application of bias correction)
# p = prior to historical period (optional)
# ratio = ratio variable?
# detrend = detrend prior to bias correction?
# missing.flag = replace NAs with missing.flag
# n.max = maximum number of tau-quantiles to estimate
function(obs.h, gcm.h, gcm.f, months.obs.h, months.gcm.h, months.gcm.f,
         gcm.p=NULL, months.gcm.p=NULL, ratio=TRUE, detrend=TRUE,
         missing.flag=-99.9, n.max=NULL) {

    if (length(gcm.p) == 0) {
        gcm.p <- NULL
    }
    if (length(gcm.h) == 0 || length(gcm.f) == 0) {
        stop(paste("mnDQM received either a historical or future period of length",
             "zero, which is an error. Check your configuration options",
             "'calibration.start' and 'calibration.end' (",
             getOption('calibration.start'), getOption('calibration.end'),
             ") and ensure that the GCM time covers the range and extends",
             "beyond it into the future")
         )
    }

    obs.h.months <- split(obs.h, months.obs.h)
    gcm.h.months <- split(gcm.h, months.gcm.h)
    gcm.f.months <- split(gcm.f, months.gcm.f)
    indices.f.months <- split(seq_along(gcm.f), months.gcm.f)
    indices.h.months <- split(seq_along(gcm.h), months.gcm.h)
    gcm.f.correct <- gcm.f*0
    gcm.h.correct <- gcm.h*0

    if(!is.null(gcm.p)){
        indices.p.months <- split(seq_along(gcm.p), months.gcm.p)
        gcm.p.months <- split(gcm.p, months.gcm.p)
        gcm.p.correct <- gcm.p*0
    }
    months <- sort(unique(c(months.gcm.f, months.gcm.p)))
    for(mn in months){
        bc <- DQM(na.omit(obs.h.months[[mn]]), na.omit(gcm.h.months[[mn]]),
                  na.omit(gcm.f.months[[mn]]), ratio, detrend, n.max)
        indices.h.valid <- indices.h.months[[mn]]
        indices.h.valid <- indices.h.valid[!is.na(gcm.h.months[[mn]])]
        indices.f.valid <- indices.f.months[[mn]]
        indices.f.valid <- indices.f.valid[!is.na(gcm.f.months[[mn]])]
        gcm.h.correct[indices.h.valid] <- bc$g.h.bc
        gcm.f.correct[indices.f.valid] <- bc$g.f.bc
        if(!is.null(gcm.p)){
            bc <- DQM(na.omit(obs.h.months[[mn]]), na.omit(gcm.h.months[[mn]]),
                      na.omit(gcm.p.months[[mn]]), ratio, detrend, n.max)
            indices.p.valid <- indices.p.months[[mn]]
            indices.p.valid <- indices.p.valid[!is.na(gcm.p.months[[mn]])]
            gcm.p.correct[indices.p.valid] <- bc$g.f.bc
        }
    }
    gcm.h.correct[is.na(gcm.h.correct)] <- missing.flag
    gcm.f.correct[is.na(gcm.f.correct)] <- missing.flag
    if(!is.null(gcm.p)){
        gcm.p.correct[is.na(gcm.p.correct)] <- missing.flag
        out <- list(g.h.bc=gcm.h.correct, g.f.bc=gcm.f.correct,
                    g.p.bc=gcm.p.correct)
    } else{
        out <- list(g.h.bc=gcm.h.correct, g.f.bc=gcm.f.correct)
    }
    out
}

##******************************************************************************

