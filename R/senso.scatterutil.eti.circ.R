"senso.scatterutil.eti.circ" <-
function (x, y, label, clabel, origin = c(0, 0), boxes = FALSE,col="black") {
    if (is.null(label))    return(invisible())
    if (any(is.na(label)))    return(invisible())
    if (any(label == ""))    return(invisible())
    xref <- x - origin[1]
    yref <- y - origin[2]
    for (i in 1:(length(x))) {
        cha <- as.character(label[i])
        cha <- paste(" ", cha, " ", sep = "")
        cex0 <- par("cex") * clabel
        xh <- strwidth(cha, cex = cex0)
        yh <- strheight(cha, cex = cex0) * 5/6
        if ((xref[i] > yref[i]) & (xref[i] > -yref[i])) {
            x1 <- x[i] + xh/2
            y1 <- y[i]
        }
        else if ((xref[i] > yref[i]) & (xref[i] <= (-yref[i]))) {
            x1 <- x[i]
            y1 <- y[i] - yh
        }
        else if ((xref[i] <= yref[i]) & (xref[i] <= (-yref[i]))) {
            x1 <- x[i] - xh/2
            y1 <- y[i]
        }
        else if ((xref[i] <= yref[i]) & (xref[i] > (-yref[i]))) {
            x1 <- x[i]
            y1 <- y[i] + yh
        }
if (boxes) {
            rect(x1 - xh/2, y1 - yh, x1 + xh/2, y1 + yh, col = "white", border = 1)
        }
        text(x1, y1, cha, cex = cex0,col=col)
    }
}
