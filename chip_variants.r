#!/usr/bin/env Rscript
library(data.table)
args = commandArgs(trailingOnly=TRUE)
# the size of simulated chip (Mega base)
if (is.na(args[2])) {
    wg_chip <- 1
} else {
    wg_chip <- as.numeric(args[2])
}
print(wg_chip)

chip_maf <- data.table(maf=suppressWarnings(fread("../pier_biobank.frq")$V5))
seq_pos <- fread(paste0(args[1], ".pos"))


chip_maf <- chip_maf[maf!=0]
#rg <- c(0, 1e-4, seq(1e-3, 9.5e-2, length.out=20), seq(2e-2, 0.5, length.out=30))
#rg <- seq(0, 1, length.out=60)
rg <- c(0, 0.007, 0.01, 0.04, 0.08,  0.1, 0.2, 0.3, 0.4, 0.5)
chip_maf$range <- cut(chip_maf$maf, rg)
chip_maf$ctg <- as.integer(chip_maf$range)

names(seq_pos) <- c("pos", "an", "ac")
seq_pos$af <- seq_pos$ac / seq_pos$an
print(seq_pos)
seq_pos$maf <- seq_pos$af
seq_pos[af>0.5]$maf<- 1 - seq_pos[af>0.5]$af

seq_pos$range <- cut(seq_pos$maf, rg)
seq_pos$ctg <- as.integer(seq_pos$range)

frq <- chip_maf[, .N, ctg]
frq$frq <- frq$N / nrow(chip_maf)
setkey(frq, "ctg")


chip_size <- wg_chip / 3088
chip_size <- chip_size * max(seq_pos$pos)
# chip_size <- chip_size * nrow(seq_pos)
result <- list()
for (i in seq(nrow(frq))) {
    nsample <- ceiling(frq$frq[i] * chip_size)
    nc <- frq$ctg[i]
    spc <- seq_pos[ctg == nc]
    result[[nc]] <- spc[sample(seq(nrow(spc)), nsample, replace=F)]
    cat(nc, "\t", nsample, frq$frq[i], "\n")
}

result <- rbindlist(result)
outtb <- data.table(chr=1, pos_from=result$pos, pos_to=result$pos)
setkey(outtb, "pos_from")
fwrite(file=paste0(args[1], ".chip.variants"), outtb, col.names=F, sep="\t")

bin <- c(0, fread("../winni.bin", head=F)$V1)
seq_pos$catg <- as.integer(cut(seq_pos$af, bin))
print(seq_pos)
fwrite(data.table(seq_pos$catg - 1), col.name=F, file="freq")
