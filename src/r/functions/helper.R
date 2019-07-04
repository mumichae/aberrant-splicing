
plotNBySample <- function(resSignif, allSamples=unique(resSignif$sampleID),
                    main="Aberrant Events by sample", zeroOff=0.7){
    plotdt <- resSignif[,.(sampleID, .N),by="hgnc_symbol"][,.N,by="sampleID"][order(N)]
    setkey(plotdt, "sampleID")
    plotdt <- plotdt[J(unique(allSamples, sampleID))]
    plotdt$N <- as.double(plotdt$N)
    plotdt[is.na(N), N:=zeroOff]
    plotdt <- plotdt[order(N)]

    pos <- plotdt[,barplot(N, col=grepl("GTEX", sampleID), log="y",
            xlab="Sample Rank", ylab="Number of events by gene", main=main)]
    grid()
    abline(h=plotdt[,median(N)], col="gray70")

    axis(side=1, at=pos[,1], labels=FALSE)
    xtickIdx <- ceiling(seq(1, nrow(plotdt), length.out=20))
    axis(side=1, at=pos[xtickIdx,1], tick = FALSE, labels = xtickIdx)
}
