options(warn=-1)

invisible(if (file.exists("jcv_cluster.jpg")) file.remove("jcv_cluster.jpg"))
invisible(if (file.exists("jcv_cluster.mx")) file.remove("jcv_cluster.mx"))
invisible(if (file.exists("jcv_cluster.sif")) file.remove("jcv_cluster.sif"))
invisible(if (file.exists("jcv_cluster.noa")) file.remove("jcv_cluster.noa"))
invisible(if (file.exists("clusters.txt")) file.remove("clusters.txt"))
invisible(if (file.exists("stats.txt")) file.remove("stats.txt"))

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
        stop("Wrong parameters. Run as: Rscript JaccardClusters.R <species-ortholog list> <estimated number of clusters> <list of outlier species>", call.=FALSE)
}

jcv <- function(list1,list2)
{
        length(intersect(list1,list2))/(length(list1)+length(list2)-length(intersect(list1,list2)))
}

message("Reading data...")
data <- read.delim(args[1], header=F, sep="\t")
k <- as.numeric(args[2])
outliers <- as.matrix(read.delim(args[3], header=F, sep="\t"))
datamx <- as.matrix(data)
cols = c("V1","V2")
datamx2 = datamx[,cols]
#species = sort(unique(datamx2[,1]))
species = as.matrix(sort(unique(datamx2[,1])))
species_no_outliers <- species[!species %in% outliers]

jcvs=list()
message("Calculating JCV values...")
for (i in 1:length(species)) {
        for (j in 1:length(species)) {
                x=datamx2[datamx2[,1]==species[i],2]
                x=x[x!='']
                x=unique(x)
                y=datamx2[datamx2[,1]==species[j],2]
                y=y[y!='']
                y=unique(y)
                jcv_ij = jcv(x,y)
                jcvs=c(jcvs,jcv_ij)
                if (i<j) {
                        spi = species[i]
                        spj = species[j]
                        cat(spi,"pp",spj,'\n',file="jcv_clusters.sif",sep="\t",append=TRUE)
                        cat(spi,"(pp)",spj,'=',jcv_ij,'\n',file="jcv_clusters.noa",sep=" ",append=TRUE)
                }
        }
}
message("The noa and sif files are complete...")

jcvs2 = matrix(jcvs,length(species),nrow=length(species))
class(jcvs2) <- "numeric"
dim(jcvs2) <- c(length(species),length(species))
colnames(jcvs2) = species
rownames(jcvs2) = species

# matrix without outliers
jcvs2b = data.frame(jcvs2)
jcvs2c = jcvs2b[!rownames(jcvs2) %in% outliers,!colnames(jcvs2) %in% outliers]

jcvs2_no = as.matrix(jcvs2c,length(species_no_outliers),nrow=length(species_no_outliers))
class(jcvs2_no) <- "numeric"
dim(jcvs2_no) <- c(length(species_no_outliers),length(species_no_outliers))
colnames(jcvs2_no) = species_no_outliers
rownames(jcvs2_no) = species_no_outliers

message("Creating heatmap...")

myBreaks <- c(seq(0,1,by=0.01))
cexx = 1-length(species)/200
ceyy = cexx

jpeg(filename = "jcv_cluster.jpg", height = 5000, width = 5000, units = "px", res = 600)
heatmap(jcvs2, symkey =F, symbreaks=F, scale="none", dendrogram = F, Rowv=F, Colv=F,col = gray.colors(100), breaks = myBreaks, na.color="white", margin = c(12,16), 
cexRow=cexx,cexCol=ceyy, key=T, trace="none", lmat=rbind( c(4, 3), c(2,1), c(0,0) ), lhei=c(1.8,6.5, 1 ))
invisible(dev.off())

write.table(jcvs2, file="jcv_cluster.mx", sep="\t", row.names=species, col.names=species)

message("Calculating baramins and statistics...")

cluster_results = kmeans(jcvs2_no,k)
write.table(cluster_results$cluster, file="clusters.txt", col.names=F, quote=F, sep="\t")

header = "baramin\tspecies\tmean\tstdev\tmin\tmax\tp-value"
write(header, file="stats.txt", sep="\t", append=T)

cluster_sizes = cluster_results$size
for (n_cluster in 1:k) {
        csize = cluster_sizes[n_cluster]
        if (csize >= 3) {
                m1 = as.matrix(jcvs2_no[cluster_results$cluster == n_cluster,cluster_results$cluster == n_cluster])
                x = m1[upper.tri(m1)]
                ll = dim(m1)[1]

                m2 = as.matrix(cbind(jcvs2_no[cluster_results$cluster != n_cluster,cluster_results$cluster == n_cluster],t(jcvs2_no[cluster_results$cluster == n_cluster,cluster_results$cluster != n_cluster])))
                m2b = m2[!duplicated(colnames(m2))]

                t = t.test(x,m2b)
                pval = t$p.value
                min = min(x)
                max = max(x)

                mean2 = sprintf("%.3f", mean(x))
                sd2 = sprintf("%.3f", sd(x))
                min2 = sprintf("%.3f", min)
                max2 = sprintf("%.3f", max)
                pval2 = sprintf("%.3f", pval)

                stats = paste(n_cluster, ll, mean2, sd2, min2, max2, pval, sep="\t")
                stats2 = gsub("\n","\t",stats)
                write(stats, file="stats.txt", sep="\t", append = T)
        }
}

message("Complete!")
