# ADMIXTURE results are visualized

# copy - pasted values for CV error depending on number of clusters
K <- c(0.90734,0.85226,0.80625,0.78304,0.76245,0.75030,0.73184,0.71927,0.71094,0.70559)

png(filename = "ADMIXTURE_CV_error.png", width = 500, height = 500, res=100)
plot(K,xlab="Number of Clusters", ylab="Cross-Validation error")#, main="Without MAF filter")
dev.off()
