load('/path/to/ELIGOS_PRcurve_window_50bp.Rdata')
pr_ELIGOS <- pr
load('/path/to/m6Anet_PRcurve_window_50bp.Rdata')
pr_m6anet <- pr

ELIGOS_recall <- pr_ELIGOS$curve[,1]
ELIGOS_precision <- pr_ELIGOS$curve[,2]

# identify which row of the pr_ELIGOS$curve has the adjusted p-value closest to 0.0001 
difference <- 1
value <- NA
index <- NA
for (i in 1:length(pr_ELIGOS$curve[,3])) {
  dif <- abs(pr_ELIGOS$curve[,3][i] + 0.0001)
  if (dif < difference) {
    difference <- dif
    value <- pr_ELIGOS$curve[,3][i]
    index <- i
  }
}

pr_ELIGOS$curve[index,]

# identify which row of the pr_ELIGOS$curve has the adjusted p-value closest to 0.05 
difference <- 1
value <- NA
index <- NA
for (i in 1:length(pr_ELIGOS$curve[,3])) {
  dif <- abs(pr_ELIGOS$curve[,3][i] + 0.05)
  if (dif < difference) {
    difference <- dif
    value <- pr_ELIGOS$curve[,3][i]
    index <- i
  }
}

pr_ELIGOS$curve[index,]

m6anet_recall <- pr_m6anet$curve[,1]
m6anet_precision <- pr_m6anet$curve[,2]

# identify which row of the pr_m6anet$curve has the probability of modification closest to 0.75
difference <- 1
value <- NA
index <- NA
for (i in 1:length(pr_m6anet$curve[,3])) {
  dif <- abs(pr_m6anet$curve[,3][i] - 0.75)
  if (dif < difference) {
    difference <- dif
    value <- pr_m6anet$curve[,3][i]
    index <- i
  }
}

pr_m6anet$curve[index,]

# identify which row of the pr_m6anet$curve has the probability of modification closest to 0.9
difference <- 1
value <- NA
index <- NA
for (i in 1:length(pr_m6anet$curve[,3])) {
  dif <- abs(pr_m6anet$curve[,3][i] - 0.9)
  if (dif < difference) {
    difference <- dif
    value <- pr_m6anet$curve[,3][i]
    index <- i
  }
}

pr_m6anet$curve[index,]

pdf('/path/to/PR_ELIGOS_m6Anet.pdf')
plot(ELIGOS_recall, ELIGOS_precision,type = 'l', col='red',xlab='Recall', ylab = 'Precision', main='ELIGOS and m6Anet recall vs precision')
points(x = c(0.0378195648,0.07775195), y= c(0.9421052632,0.83023125), pch = 15, cex = 2,  col=c('green', 'orange'))
lines(m6anet_recall, m6anet_precision, type='l', col='grey')
points(x = c(0.09539404,0.1595711), y= c(0.88616290,0.7730297), pch = 15, cex = 2, col=c('black', 'violet'))
legend(x='topright', legend = c('ELIGOS', 'm6Anet', 'adj.Pvalue ELIGOS < 0.0001', 'adj.Pvalue ELIGOS < 0.05', 
                                'Probability of modification m6Anet > 0.9', 'Probability of modification m6Anet > 0.75' ),
       pch = c(NA,NA,15,15,15,15),lty = c(1,1,NA,NA,NA,NA),cex=1,col = c('red','grey','green', 'orange','black', 'violet'))
dev.off()

