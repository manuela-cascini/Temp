fst_plots_mc <- function (pFst, pS) 
{
    outputloc <- paste0(RandRbase, species, "/output_", analysis, 
        "/")
    if (!dir.exists(outputloc)) {
        dir.create(outputloc)
    }
    tempdir <- paste0(outputloc, "temp/")
    if (!dir.exists(tempdir)) {
        dir.create(tempdir)
    }
    fst_dir <- paste0(outputloc, "fst")
    if (!dir.exists(fst_dir)) {
        dir.create(fst_dir)
    }
    library(reshape2)
    library(vegan)
    library(ade4)
    library(ggplot2)
    library(cowplot)
    tiff(paste0(fst_dir, "/", species, " fst plot.tiff"), units = "in", 
        width = 10, height = 5, res = 300)
    par(mfrow = c(1, 2))
    diag(pS$S) <- NA
    diag(pFst$Fst) <- NA
    Fst_sig <- cbind(melt(pS$S), unlist(as.list(pFst$Fst)))
    colnames(Fst_sig)[3] <- "Geo_dist"
    colnames(Fst_sig)[4] <- "Fst"
    Fst_sig$Geo_dist2 <- Fst_sig$Geo_dist/1000
    plot(Fst_sig$Geo_dist2, Fst_sig$Fst, xlab = "distance (km)", 
        ylab = "Fst", cex = 1, font = 4, cex.main = 1)
    title(main = paste0(species, " pairwise fst plots"), adj = 0.001, 
        font.main = 4)
    man <- mantel(xdis = pS$S, ydis = pFst$Fst, permutations = 999, 
        na.rm = TRUE)
    legend("bottomright", bty = "n", cex = 0.75, text.col = "blue", 
        legend = paste("Mantel statistic r is ", format(man$statistic, 
            digits = 4), " P =", format(man$signif)))
    Fst_sig <- cbind(melt(pS$S), unlist(as.list(pFst$Fst)))
    colnames(Fst_sig)[3] <- "Geo_dist"
    colnames(Fst_sig)[4] <- "Fst"
    Fst_sig$lin_fst <- Fst_sig$Fst/(1 - Fst_sig$Fst)
    Fst_sig$Geo_dist2 <- Fst_sig$Geo_dist/1000
    Fst_sig$log10Geo_dist <- log10(Fst_sig$Geo_dist2)
    write.table(Fst_sig, paste0(fst_dir, "/", species, " fst_geo data.csv"))
    plot(Fst_sig$log10Geo_dist, Fst_sig$lin_fst, xlab = "log10(distance)", 
        ylab = "Linearised_Fst", font = 4, cex.main = 1)
    man2 <- mantel(xdis = pS$S, ydis = (pFst$Fst)/(1 - pFst$Fst), 
        permutations = 999, na.rm = TRUE)
    legend("bottomright", bty = "n", cex = 0.75, text.col = "blue", 
        legend = paste("Mantel statistic r is ", format(man2$statistic, 
            digits = 4), " P =", format(man2$signif)))
    dev.off()
    cat(paste0("pairwise fst and linearised fst drawn!", "\n"))
    par(mfrow = c(2, 1), oma = c(0, 0, 1, 0))
    geo_d <- pS$S
    geo_d[upper.tri(geo_d)] <- NA
    rownames(geo_d) <- colnames(pS$S)
    dimnames <- list(var1 = colnames(pS$S), var2 = colnames(pS$S))
    mat <- matrix(geo_d, ncol = length(colnames(geo_d)), nrow = length(colnames(geo_d)), 
        dimnames = dimnames)
    df <- as.data.frame(as.table(mat))
    p1 <- ggplot(df, aes(var1, var2)) + geom_tile(aes(fill = Freq), 
        colour = "white") + scale_fill_gradient(low = "white", 
        high = "steelblue") + theme(axis.text.x = element_text(angle = 90, 
        hjust = 1)) + xlab("") + ylab("") + theme(legend.position = "none", 
        axis.text = element_text(size = 5), plot.title = element_text(face = "bold.italic")) + 
        ggtitle(paste0(species, " heatmaps", "\n", "heatmap of pairwise geographic distance"))
    genetic_d <- pFst$Fst
    genetic_d[upper.tri(genetic_d)] <- NA
    rownames(genetic_d) <- colnames(pFst$Fst)
    dimnames2 <- list(var1 = colnames(pFst$Fst), var2 = colnames(pFst$Fst))
    mat2 <- matrix(genetic_d, ncol = length(colnames(geo_d)), 
        nrow = length(colnames(geo_d)), dimnames = dimnames)
    df2 <- as.data.frame(as.table(mat2))
    df3 <- df2[complete.cases(df2$Freq), ]
    p2 <- ggplot(df3, aes(var1, var2)) + geom_tile(aes(fill = Freq), 
        colour = "white", na.rm = TRUE) + scale_fill_gradient(low = "white", 
        high = "red") + geom_text(aes(label = round(Freq, 3)), 
        size = 3, df3) + theme(axis.text.x = element_text(angle = 90, 
        hjust = 1)) + xlab("") + ylab("") + theme(legend.position = "none", 
        axis.text = element_text(size = 5), plot.title = element_text(face = "bold.italic")) + 
        ggtitle("heatmap of pairwise Fst")
    plot_grid(p1, p2, ncol = 1)
    ggsave(paste0(fst_dir, "/", species, " fst heatmaps.tiff"), 
        width = 8, height = 9, dpi = 300, units = "in", device = "tiff")
    cat(paste0("geographic and pairwise fst heatmap drawn!", 
        "\n"))
    
    genetic_d <- as.matrix(as.dist(pFst$Fst))
    write.csv(genetic_d, "Rmaiden_pairwisefst.csv") 

    #genetic_d <- pFst$Fst
    #write.table(genetic_d, paste0(fst_dir, "/", species, " fst matrix.csv"), 
     #   sep = ",")
    geo_d <- pS$S
    geo_d[upper.tri(geo_d)] <- geo_d[lower.tri(geo_d)]
    new <- matrix(NA, nrow = dim(geo_d)[1], ncol = dim(geo_d)[2])
    new[upper.tri(new)] <- geo_d[upper.tri(geo_d)]
    new[lower.tri(new)] <- genetic_d[lower.tri(genetic_d)]
    colnames(new) <- colnames(geo_d)
    rownames(new) <- rownames(geo_d)
    write.table(new, paste0(fst_dir, "/", species, " heatmap matrix.csv"), 
        sep = ",")
}
