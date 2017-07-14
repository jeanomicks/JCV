options(warn=-1)

invisible(if (file.exists("jcv_clusters.jpg")) file.remove("jcv_clusters.jpg"))
invisible(if (file.exists("jcv_clusters.mx")) file.remove("jcv_clusters.mx"))
invisible(if (file.exists("jcv_clusters.sif")) file.remove("jcv_clusters.sif"))
invisible(if (file.exists("jcv_clusters.noa")) file.remove("jcv_clusters.noa"))
invisible(if (file.exists("clusters.txt")) file.remove("clusters.txt"))
invisible(if (file.exists("clusters_pgq.txt")) file.remove("clusters.txt"))
invisible(if (file.exists("clusters_mxcut.txt")) file.remove("clusters.txt"))
invisible(if (file.exists("stats.txt")) file.remove("stats.txt"))
invisible(if (file.exists("stats_pgq.txt")) file.remove("stats.txt"))
invisible(if (file.exists("stats_mxcut.txt")) file.remove("stats.txt"))

args = commandArgs(trailingOnly=TRUE)
if ((args[1]=="k-means")&&(!(length(args)==4))) {
        stop("Wrong parameters. Run as: Rscript JaccardClustersMulti.R k-means <species-ortholog list> <list of outlier species> <estimated no. clusters>", call.=FALSE)
} else {
	  k = args[4]
	  if (!(k == round(k))) {
		stop("Cluster number not an integer!")
	  }
}
if ((args[1]=="mxcut")&&(!(length(args)==4))) {
        stop("Wrong parameters. Run as: Rscript JaccardClustersMulti.R mxcut <species-ortholog list> <list of outlier species> <JCV cutoff>", call.=FALSE)
} else {
	  if ((args[4]<0)||(args[4]>1)) {
		stop("Cutoff not between 0 and 1")
	  }
}
if ((args[1]=="pqg")&&(!(length(args)==3))) {
        stop("Wrong parameters. Run as: Rscript JaccardClustersMulti.R pgq <species-ortholog list> <list of outlier species>", call.=FALSE)
}
if (!(args[1] %in% c("k-means","pgq","mxcut"))) {
	  stop("First parameter must be either k-means, pgq, or mxcut !")
}

jcv <- function(list1,list2)
{
        length(intersect(list1,list2))/(length(list1)+length(list2)-length(intersect(list1,list2)))
}

jcv_multi <- function(mx, list)
{
	sum(as.matrix(table(datamx2[datamx2[,1] %in% list,2]),2)==length(list))/length(unique(mx[mx[,1] %in% list,2]))
}

get_p_val <- function(mx_no, allspecs_no, specs)
{
	m1 = as.matrix(mx_no[specs, specs])
	x = m1[upper.tri(m1)]
	non_specs = allspecs_no[which(!allspecs_no %in% specs)]
	m2 = as.matrix(mx_no[specs, non_specs])
	t = t.test(x,m2)
	pval = t$p.value
	pval
}

vec2str <- function(vec)
{
	paste(vec, collapse=" ")
}

str2vec <- function(str)
{
	as.list(strsplit(str," "))[[1]]
}

mxcuts <- function(mx, cut, s)
{
	list <- names(which(mx[s,] >= cut))
	list <- sort(unique(list))
	while (length(s) < length(list)) {
		s = list
		lists = c()
		for (i in 1:length(list)) {
			elem = list[i]
			listx <- names(which(mx[elem,] >= cut))
			lists = c(lists, listx)
		}
		list <- sort(unique(lists))
	}
	list
}

cluster_members <- function(datamx, jcvs, spec)
{
	rr=as.matrix(rev(sort(jcvs[spec,])))
	rrr=rownames(rr)
	jcvmp=1
	djcvm=1
	bmembers=c()
	for (i in 1:length(rrr)) {
		l=rrr[1:i]
		jcvm=jcv_multi(datamx,l)
		dd = 0
		if (i>1) {
			dd = abs(jcvm - jcvmp)
			if (dd > djcvm) {
				bmembers = l[1:i-1]
			}
		}
		jcvmp=jcvm
		djcvm = dd
	}
	sort(bmembers)
}

message("Reading data...")

progname <- args[1]
data <- read.delim(args[2], header=F, sep="\t")
outliers <- as.matrix(read.delim(args[3], header=F, sep="\t"))
param <- args[4]

datamx <- as.matrix(data)
cols = c("V1","V2")
datamx2 = datamx[,cols]
datamx2 = unique(datamx2)
species = as.matrix(sort(unique(datamx2[,1])))
species_no_outliers <- species[!species %in% outliers]
datamx2_no = matrix(datamx2[datamx2[,1] %in% species_no_outliers],ncol=2)

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

jpeg(filename = "jcv_clusters.jpg", height = 5000, width = 5000, units = "px", res = 600)
heatmap(jcvs2, symkey =F, symbreaks=F, scale="none", dendrogram = F, Rowv=F, Colv=F,col = gray.colors(100), breaks = myBreaks, na.color="white", margin = c(12,16), 
cexRow=cexx,cexCol=ceyy, key=T, trace="none", lmat=rbind( c(4, 3), c(2,1), c(0,0) ), lhei=c(1.8,6.5, 1 ))
invisible(dev.off())

write.table(jcvs2, file="jcv_clusters.mx", sep="\t", row.names=species, col.names=species)

message("Calculating baramins and statistics...")

# PGQ
if (args[1] == "pgq") {
	cluster_list = c()
	for (s in species_no_outliers) {
		sx = cluster_members(datamx2_no, jcvs2_no, s)
		pgq = jcv_multi(datamx2_no,sx)
		sxpgq = c(as.character(sx),as.character(pgq))
		sxpgqs = vec2str(sxpgq)
		cluster_list = c(cluster_list,sxpgqs)
	}
	cluster_list <- unique(cluster_list)

	# write states file
	header = "baramin\tspecies\tmean\tstdev\tmin\tmax\tPGQ\tp-value"
	write(header, file="stats_pgq.txt", sep="\t", append=T)
	pgqc = 0
	for (param in cluster_list) {
		pgqc = pgqc + 1
		sxpgqv = str2vec(param)
		specs = sxpgqv[1:length(sxpgqv)-1]
		if (length(specs) >= 3) {
			pgq = sxpgqv[length(sxpgqv)]
			submx = jcvs2[specs, specs]

			submx_uppertr = submx[upper.tri(submx)]
			mean2 = sprintf("%.3f", mean(submx_uppertr))
           		sd2 = sprintf("%.3f", sd(submx_uppertr))
           		min2 = sprintf("%.3f", min(submx_uppertr))
           		max2 = sprintf("%.3f", max(submx_uppertr))
			p = get_p_val(jcvs2_no, species_no_outliers, specs)

	           	stats = paste(pgqc, length(specs), mean2, sd2, min2, max2, pgq, p, sep="\t")
	           	stats2 = gsub("\n","\t",stats)
	           	write(stats, file="stats_pgq.txt", sep="\t", append = T)

			for (spec in specs) {
				isp = paste(spec,pgqc,sep="\t")
				isp2 = gsub("\n","\t",isp)
				write(isp2, "clusters_pgq.txt", sep="\t", append=T)
			}
		}
	}
}

# mxcuts
if (args[1] == "mxcut") {
	message("Running mxcuts...")
	cut = args[4]
	allnames = species_no_outliers
	s = allnames[1]

	header = "baramin\tspecies\tmean\tstdev\tmin\tmax\tp-value"
	write(header, file="stats_mxcut.txt", sep="\t", append=T)

	nclust = 0
	while(length(allnames) > 0){
		xlist = mxcuts(jcvs2_no,cut,s)
		lxlist = length(xlist)
		if (lxlist >= 3) {
			nclust = nclust + 1
			submx = jcvs2[xlist, xlist]
			submx_uppertr = submx[upper.tri(submx)]
			mean2 = sprintf("%.3f", mean(submx_uppertr))
           		sd2 = sprintf("%.3f", sd(submx_uppertr))
           		min2 = sprintf("%.3f", min(submx_uppertr))
           		max2 = sprintf("%.3f", max(submx_uppertr))
			p = get_p_val(jcvs2_no, species_no_outliers, xlist)
			stats = paste(nclust, lxlist, mean2, sd2, min2, max2, p, sep="\t")
			stats2 = gsub("\n","\t",stats)
			write(stats2, file="stats_mxcut.txt", sep="\t", append = T)
			for (i in xlist) {
				isp = paste(i,nclust,sep="\t")
				isp2 = gsub("\n","\t",isp)
				write(isp2, "clusters_mxcut.txt", sep="\t", append=T)
			}
		}
		allnames = allnames[!allnames %in% xlist]
		s=allnames[1]
	}
}

# k-means
if (args[1] == "k-means") {
	k = args[4]

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
}

message("Complete!")
