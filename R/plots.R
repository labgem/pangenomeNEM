#!/usr/bin/env R
options(show.error.locations = TRUE)

library("ggplot2")
library("reshape2")
library("data.table")
library("ggrepel")
library("grDevices")
library("RColorBrewer")

ggcolors <- function(n = 6){
    h = c(0, 360) + 15
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    grDevices::hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

args = commandArgs(trailingOnly=TRUE)

color_chart = c(pangenome="black", accessory="#808080", shell = "#F0E442", persistant="#CC79A7", strict_core ="#E0162B", cloud = "#80CCFF")

#cloud = Bleu
#shell = Jaune foncé
#persistant = Rouge clair
#strict_core = Rouge foncé
#pangenome = noir
#accesoire = gris

#EVOLUTION
print(paste0(args[1],"/","evolution_stats_nem_3.txt"))
data <- read.table(paste0(args[1],"/","evolution_stats_nem_3.txt"))
print(nrow(data))
if(nrow(data)>=2){
	colnames(data) <- c("c","persistant","shell","cloud","pangenome")
	data <- melt(data, id = "c")

	colnames(data) <- c("c","cluster","value")


	max_c <- max(data$c,na.rm=T)
	print(max_c)
	final_state = data[data$c == max_c,]
	print(final_state)

	final_pangenome = final_state[final_state$cluster == "pangenome", "value"]
	final_shell = final_state[final_state$cluster == "shell", "value"]
	final_cloud = final_state[final_state$cluster == "cloud", "value"]
	final_persistant = final_state[final_state$cluster == "persistant", "value"]

	final<- as.vector(c(final_pangenome,final_shell,final_cloud,final_persistant))
	names(final) <- c("pangenome","shell","cloud","persistant")
	print(final)
	print(names(sort(final, decreasing = TRUE)))

	#gamma and kappa are calculated according to the Tettelin et al. 2008 approach
	median_by_comb <- setDT(data)[,list(med=as.numeric(median(value))),by=c("c","cluster")]
	print(median_by_comb)
	colnames(median_by_comb) <- c("comb","cluster","med")
	regression <- nls(med~kapa*(comb^gama),median_by_comb[median_by_comb$cluster == "pangenome",],start=list(kapa=1000,gama=1))
	print(regression)
	coefficient <- coef(regression)
	kappa  <- coefficient["kapa"]
	gamma <- coefficient["gama"]


	p <- ggplot(data = data, aes_string(x="c",y="value", colour = "cluster"))+
		ggtitle(bquote(list("Rarefaction curve. Heaps-law parameters based on Tettelin et al. 2008 approach", kappa==.(kappa), gamma==.(gamma))))+
		geom_smooth(data = median_by_comb[median_by_comb$cluster %in% c("pangenome","shell","cloud") ,], aes_string(x="comb",y="med",colour = "cluster"), method="nls",formula=y~kapa*(x^gama),method.args =list(start=c(kapa=1000,gama=1)), linetype="twodash",size = 1.5,se=FALSE, show.legend = FALSE)+       
		stat_summary(fun.ymin = function(z) { quantile(z,0.25) },  fun.ymax = function(z) { quantile(z,0.75) }, geom="ribbon", alpha=0.1,size=0.1, linetype="dashed", show.legend = FALSE)+
		stat_summary(fun.y=median, geom="line",size=0.5)+
		stat_summary(fun.y=median, geom="point",shape=4,size=1, show.legend = FALSE)+
		stat_summary(fun.ymax=max,fun.ymin=min,geom="errorbar",linetype="dotted",size=0.1,width=0.2)+
		scale_x_continuous(breaks = as.numeric(unique(data$c)))+
		scale_y_continuous(limits=c(0,max(data$value,na.rm=T)), breaks = seq(0,max(data$value,na.rm=T),1000))+
		scale_colour_manual(name = "nem classification", values = color_chart, breaks=names(sort(final, decreasing = TRUE)))+
		geom_label_repel(data = final_state, aes_string(x="c", y="value", colour = "cluster", label = "value"), show.legend = FALSE,
				  fontface = 'bold', fill = 'white',
				  box.padding = unit(0.35, "lines"),
				  point.padding = unit(0.5, "lines"),
				  segment.color = 'grey50',
				  nudge_x = 45) +
		xlab("# of organism")+
		ylab("# of MICFAM")+
		ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

	ggsave(paste0(args[1],"/figures/","evolution.pdf"), device = "pdf", width = (par("din")[1]*2) ,p)

	data <- read.table(paste0(args[1],"/","evolution_stats_exact.txt")) 
	colnames(data) <- c("c","strict_core","accessory","pangenome")
	data <- melt(data, id = "c")
	colnames(data) <- c("c","cluster","value")


	max_c <- max(data$c,na.rm=T)
	print(max_c)
	final_state = data[data$c == max_c,]

	final_pangenome = final_state[final_state$cluster == "pangenome", "value"]
	final_accessory = final_state[final_state$cluster == "accessory", "value"]
	final_core_strict = final_state[final_state$cluster == "strict_core", "value"]

	final<- c(final_pangenome,final_accessory,final_core_strict)
	names(final) <- c("pangenome","accessory","strict_core")

	median_by_comb <- setDT(data)[,list(med=as.numeric(median(value))),by=c("c","cluster")]
	colnames(median_by_comb) <- c("comb","cluster","med")

	p <- ggplot(data = data, aes_string(x="c",y="value", colour = "cluster"))+
		ggtitle(bquote(list("Rarefaction curve. Heaps-law parameters based on Tettelin et al. 2008 approach", kappa==.(kappa), gamma==.(gamma))))+
		geom_smooth(data = median_by_comb[median_by_comb$cluster %in% c("pangenome","accessory") ,], aes_string(x="comb",y="med",colour = "cluster"), method="nls",formula=y~kapa*(x^gama),method.args =list(start=c(kapa=1000,gama=1)),linetype="twodash",size = 1.5,se=FALSE, show.legend = FALSE)+
		stat_summary(fun.ymin = function(z) { quantile(z,0.25) },  fun.ymax = function(z) { quantile(z,0.75) }, geom="ribbon", alpha=0.1,size=0.1, linetype="dashed", show.legend = FALSE)+
		stat_summary(fun.y=median, geom="line",size=0.5)+
		stat_summary(fun.y=median, geom="point",shape=4,size=1, show.legend = FALSE)+
		stat_summary(fun.ymax=max,fun.ymin=min,geom="errorbar",linetype="dotted",size=0.1,width=0.2)+
		scale_x_continuous(breaks = as.numeric(unique(data$c)))+
		scale_y_continuous(limits=c(0,max(data$value,na.rm=T)), breaks = seq(0,max(data$value,na.rm=T),1000))+
		scale_color_manual(name = "exact classification", values = color_chart, breaks=names(sort(final, decreasing = TRUE)))+	
		geom_label_repel(data = final_state, aes_string(x="c", y="value", colour = "cluster", label = "value"), show.legend = FALSE,
				  fontface = 'bold', fill = 'white',
				  box.padding = unit(0.35, "lines"),
				  point.padding = unit(0.5, "lines"),
				  segment.color = 'grey50',
				  nudge_x = 45) +
		xlab("# of organism") +
		ylab("# of MICFAM")+
		ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


	ggsave(paste0(args[1],"/figures/","evolution_exact.pdf"), device = "pdf", width = (par("din")[1]*2) ,p)

}else{
	max_c <- data[1,1]
}
#CLUSTER

binary_matrix <- read.table(paste0(args[1],"NEM_results/nborg",max_c,"_k3_i0/file.dat"), header=FALSE)

head(binary_matrix)
binary_matrix <- ifelse(binary_matrix != 0, TRUE, FALSE)


occurences <- rowSums(binary_matrix)

head(occurences)

classification_vector <- unlist(strsplit(readLines(paste0(args[1],"NEM_results/nborg",max_c,"_k3_i0/file.cf")), " "))


head(classification_vector)


means <- data.frame(cluster = c("1","2","3"), mean = rep(NA,3))

means[means$cluster == "1","mean"] <- mean(occurences[classification_vector == "1"])
means[means$cluster == "2","mean"] <- mean(occurences[classification_vector == "2"])
means[means$cluster == "3","mean"] <- mean(occurences[classification_vector == "3"])

means <- means[order(means$mean),]

classification_vector[classification_vector == means[1,"cluster"]] <- "cloud"
classification_vector[classification_vector == means[2,"cluster"]] <- "shell"
classification_vector[classification_vector == means[3,"cluster"]] <- "persistant"

c = data.frame(nb_org = occurences, cluster = classification_vector)

plot <- ggplot(data = c) + 
	geom_bar(aes_string(x = "factor(nb_org)", fill = "cluster")) +
	coord_flip() +
	scale_fill_manual(name = "nem classification", values = color_chart, breaks=c("persistant","shell","cloud")) +
	xlab("# of organims")+
	ylab("# of MICFAM")

ggsave(paste0(args[1],"/figures/","clusters.pdf"), device = "pdf", height= (par("din")[2]*1.5),plot)

#cluster with GO
onto <- paste0(args[1],"ontology.txt")
if (file.exists(onto)){
	go <- read.csv(onto,row.names=1)
	print(head(go))
	go_mf <- go$"Molecular Function"

	c_go = data.frame(nb_org = occurences, go_mf = go_mf)

	plot <- ggplot(data = c_go) + 
		geom_bar(aes_string(x = "factor(nb_org)", fill = "go_mf")) +
		coord_flip() +
		xlab("# of organims")+
		ylab("# of MICFAM")

	ggsave(paste0(args[1],"/figures/","ontology.pdf"), device = "pdf", height= (par("din")[2]*1.5),plot)
}

#cluster with COG
cog <- paste0(args[1],"COG.txt")
if (file.exists(cog)){

	funcat <- read.csv(cog,row.names=1)
	print(head(funcat))

	funcat[funcat$funcat == "Unknown","funcat"] <- "S"
	print(head(funcat))
	c_funcat = data.frame(nb_org = occurences, funcat = funcat)

	print(head(c_funcat))

	funcat$funcat

	cog_symbol <- sort(setdiff(unique(funcat$funcat),"S"))
	COG_colors <- colorRampPalette(brewer.pal(9, "Set1"))((length(cog_symbol)))
	names(COG_colors) <- cog_symbol
	COG_colors <- c(COG_colors,c("S"="grey"))

	plot <- ggplot(data = c_funcat) + 
		geom_bar(aes_string(x = "factor(nb_org)", fill = "funcat")) +
		scale_fill_manual(values = COG_colors)+
		coord_flip() +
		xlab("# of organims")+
		ylab("# of MICFAM")

	ggsave(paste0(args[1],"/figures/","cog.pdf"), device = "pdf", height= (par("din")[2]*3), width= (par("din")[2]*3),plot)
}
#circular plot

classification = data.frame(classification = classification_vector)

classification$classification <- reorder(classification$classification, X = classification$classification, FUN = function(x) length(x))

at <- nrow(classification) - as.numeric(cumsum(sort(table(classification)))-0.5*sort(table(classification)))

label=paste0(round(sort(table(classification))/sum(table(classification)),2) * 100,"%")

plot <- ggplot(data = classification) +
		geom_bar(aes(x="", fill = classification), color='black', width = 1) +
		scale_fill_manual(values = color_chart) +
        coord_polar(theta = "y")+
		annotate(geom = "text", y = at, x = 1, label = label) +
		theme(panel.border = element_blank(),
		panel.background = element_blank(),
		axis.text = element_blank(),
		axis.ticks = element_blank(),
        panel.grid  = element_blank(),
		axis.title.x=element_blank(),
          axis.title.y=element_blank())

ggsave(paste0(args[1],"/figures/","classification.pdf"), device = "pdf", plot)
if (file.exists(cog)){
	cog_func = data.frame(classification = classification_vector, funcat = funcat)

	print(head(cog_func))

	#circular plot COG
	for (class in c("persistant","shell","cloud")){

		cog_func_class = cog_func[cog_func$classification == class,"funcat", drop=FALSE]

		cog_func_class$funcat <- reorder(cog_func_class$funcat, X = cog_func_class$funcat, FUN = function(x) length(x))
		at <- nrow(cog_func_class) - as.numeric(cumsum(sort(table(cog_func_class)))-0.5*sort(table(cog_func_class)))
		occ <- sort(table(cog_func_class))
		percent=round(occ/sum(occ),3) * 100
		label = paste0(names(occ)," = ",percent, "%")
		x = ifelse(percent>1, 1, 1.5)

		percent_label = data.frame(x= x,y = at, label = label, colour = names(occ))
	
		print(label)
		print(COG_colors)
	
		plot <- ggplot() +
			geom_bar(data = cog_func_class, aes(x="", fill = funcat), colour = "black", width = 1) +
			coord_polar(theta = "y")+
			geom_label_repel(data = percent_label, aes_string(x="x",y="y",label="label", colour="colour", size="(1/x)*2"),fill = "white", segment.color = 'black', show.legend=FALSE)+
			scale_fill_manual(values = COG_colors) +
			scale_color_manual(values = COG_colors) +
			theme(panel.border = element_blank(),
				panel.background = element_blank(),
				axis.text = element_blank(),
				axis.ticks = element_blank(),
        			panel.grid  = element_blank(),
				axis.title.x=element_blank(),
			        axis.title.y=element_blank())

		ggsave(paste0(args[1],"/figures/","cog_circular_",class,".pdf"), device = "pdf", height= (par("din")[2]*3), width= (par("din")[2]*3), plot)
	}
}


#tile_plot

organism_names <- unlist(strsplit(readLines(paste0(args[1],"/organisms.txt")), "\n"))

colnames(binary_matrix) <- organism_names
#binary_matrix = rbind(binary_matrix,colSums(binary_matrix))

#binary_matrix = binary_matrix[-nrow(binary_matrix),order(tail(binary_matrix,n=1))]

#binary_matrix = cbind(binary_matrix,classification = classification_vector)
#binary_matrix = binary_matrix[order(match(binary_matrix[,"classification"],c("persistant", "shell", "cloud") )),]

#binary_matrix_clust   <- hclust(dist(binary_matrix, method = "manhattan"))
#binary_matrix_clust_t <- hclust(dist(t(binary_matrix), method = "manhattan"))

#binary_matrix <- binary_matrix[,binary_matrix_clust_t$order]
binary_matrix = data.frame(binary_matrix,classification = classification_vector, check.names=FALSE)
binary_matrix = binary_matrix[order(match(binary_matrix$classification,c("persistant", "shell", "cloud") )),]
#binary_matrix <- binary_matrix[binary_matrix_clust$order,]

persistant_size = table(classification_vector)["persistant"]
print(persistant_size)
shell_size = table(classification_vector)["shell"]

col <- setdiff(colnames(binary_matrix),"classification")
total_shell = colSums(binary_matrix[binary_matrix$classification == "shell",col])
total_persistant = colSums(binary_matrix[binary_matrix$classification == "persistant",col])
total_non_cloud = colSums(binary_matrix[binary_matrix$classification != "cloud",col])
total_pan = colSums(binary_matrix[,col])

binary_matrix$familles <- seq(1,nrow(binary_matrix))
data = melt(binary_matrix, id.vars=c("familles"))

print(head(data))
colnames(data) = c("fam","org","value")

#ratio_shell = data.frame(x=col, y= round((total_shell/total_non_cloud)*shell_size+persistant_size))
ratio_persistant = data.frame(x=col, y= round((total_persistant/total_pan)*persistant_size))
print(ratio_persistant)

data$value <- factor(data$value, levels = c(TRUE,FALSE,"persistant", "shell", "cloud"), labels = c("presence","absence","persistant", "shell", "cloud"))
print("gg")
plot <- ggplot(data = data)+
        geom_raster(aes_string(x="org",y="fam", fill="value"))+
		scale_fill_manual(values = c("presence"="green","absence"="grey80",color_chart)) +
	geom_point(data = ratio_persistant,aes_string(x="x",y="y"), color = color_chart["strict_core"], size =0.5, shape=4)+
	#geom_point(data = ratio_shell,aes_string(x="x",y="y"), color = color_chart["shell"], size =0.5, shape=4)+
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.border = element_blank(), panel.background = element_blank())

ggsave(paste0(args[1],"/figures/","tile_plot.pdf"), device = "pdf", plot)

# MDS
coordMDS <- paste0(args[1],"/coordMDS_weigths.txt")
if (file.exists(coordMDS)){
	data = read.table(coordMDS, row.names=1, sep=",")

	print(data)
	colnames(data) <- c("coord1","coord2","weight")
	data$organism <- organism_names

	plot <- ggplot(data = data, aes_string(x="coord1",y="coord2", label="organism", size = "weight", color = "weight", fill = "weight"))+
	geom_point()+
	ggrepel::geom_text_repel()+
	scale_colour_gradient2()+
	scale_fill_gradient2()+
	scale_size_area()

	ggsave(paste0(args[1],"/figures/","MDS.pdf"), device = "pdf", plot)
}
# Mash tree

mash_distances <- paste0(args[1],"/mash_distance.csv")
if (file.exists(mash_distances)){
	library("ape")
	library("ggtree")

	dist <- read.table(mash_distances)
	dist <- dist[,2:ncol(dist)]
	colnames(dist) <- organism_names
	rownames(dist) <- organism_names

	nj_tree <- nj(as.dist(dist))
	ggsave(paste0(args[1],"/figures/","mash_nj_tree.pdf"), device = "pdf",ggplot(nj_tree, aes(x, y)) + geom_tree() + geom_tiplab(size=3, color="blue") + theme_tree2(), height = 20, width = 25)

	#clust <- hclust(as.dist(dist))
}
