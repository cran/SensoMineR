"senso.scatterutil.eti" <- function (x, y, label, clabel, boxes = FALSE, coul = rep(1, length(x)), font = 2, pos = 4, offset = 0.2 ) {
    if (length(label) == 0)    return(invisible())
    if (is.null(label))    return(invisible())
    if (any(label == ""))    return(invisible())
    for (i in 1:(length(x))) {
        cha <- as.character(label[i])
        cha <- paste(" ", cha, " ", sep = "")
        cex0 <- par("cex") * clabel
        x1 <- x[i]
        y1 <- y[i]
        if (boxes) {
            xh <- strwidth(cha, cex = cex0)
            yh <- strheight(cha, cex = cex0) * 5/3
            rect(x1 - xh/2, y1 - yh/2, x1 + xh/2, y1 + yh/2, col = "white", border = coul[i])
        }
        text(x1, y1, cha, cex = cex0, col = coul[i], pos = pos, offset = offset ,font = font )
    }
}
