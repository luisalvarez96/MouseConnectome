#All codes by Luis Alvarez luisalvarez.10.96@gmail.com git: https://github.com/luisalvarez96

library(tidyr)
library(dplyr)
library(igraph)
library(tictoc)
library(ggplot2)
library(reshape2)
library('dendextend')
library(stringr)
library(corrplot)

#### Net from raq colors ####


getColors('~/Downloads/engram_cormatrix_r_r2_p (1).csv', 0.1, 0.1)

NoColors <- NULL
NoColors <- rbind(NoColors, getColors('~/Downloads/engram_cormatrix_r_r2_p (1).csv', 0.1, 0.1)[,c(1,14,15)] %>% filter(ClusterSize == 1))
NoColors <- rbind(NoColors, getColors('~/Downloads/engram_cormatrix_r_r2_p (1).csv', 0.1, 0.2)[,c(1,14,15)] %>% filter(ClusterSize == 1))
NoColors <- rbind(NoColors, getColors('~/Downloads/engram_cormatrix_r_r2_p (1).csv', 0.05, 0.1)[,c(1,14,15)] %>% filter(ClusterSize == 1))
NoColors <- rbind(NoColors, getColors('~/Downloads/engram_cormatrix_r_r2_p (1).csv', 0.05, 0.2)[,c(1,14,15)] %>% filter(ClusterSize == 1))
NoColors <- NoColors[!duplicated(NoColors$acronym),]


NetfromColors <- function(file, pval, h){
  nodes <- getColors(file, pval, h)[,c(1,14,15)] %>% filter(ClusterSize > 1)
  edges <- conn_edges_connstr_rdx %>% filter(Source %in% nodes$acronym & Target %in% nodes$acronym)
  
  return(list(nodes, edges))
}


conn_nodes_p1h2colors <- NetfromColors('~/Downloads/engram_cormatrix_r_r2_p (1).csv', 0.1, 0.1)[[1]]
conn_edges_p1h2colors <- NetfromColors('~/Downloads/engram_cormatrix_r_r2_p (1).csv', 0.1, 0.1)[[2]]

tic()
perc_p1h2colors <- inv_percolation(conn_nodes_p1h2colors, conn_edges_p1h2colors)
toc()


ggplot(data.frame(perc_p1h2colors$pvals), aes(log10(p.value), MeanDegree)) + xlab('log10(CONNECTION STRENGTH)') +
  geom_point() + ggtitle('p1h2 Mean Degree') + geom_hline(yintercept = 10, linetype = "dashed", color = 'blue') + 
  geom_vline(xintercept = log10(2500), linetype = "dashed", color = 'red') + coord_cartesian(xlim = c(0, NA))

ggplot(data.frame(perc_p1h2colors$pvals), aes(log10(p.value), SCCSize)) + xlab('log10(CONNECTION STRENGTH)') +
  geom_point() + ggtitle('p1h2 REDUCED Biggest SCC') + 
  geom_vline(xintercept = log10(2500), linetype = "dashed", color = 'red') + coord_cartesian(xlim = c(0, NA))

ggplot(data.frame(perc_p1h2colors$pvals), aes(log10(p.value), WeakSize)) + xlab('log10(CONNECTION STRENGTH)') +
  geom_point() + ggtitle('p1h2 REDUCED WeakCC') +
  geom_vline(xintercept = log10(2500), linetype = "dashed", color = 'red') + coord_cartesian(xlim = c(0, NA))



nodes_info(conn_nodes_p1h2colors, conn_edges_p1h2colors %>% filter(ConnStr > 100)) %>% filter(!is.na(SCCId)) %>% group_by(SCCSize) %>% summarise(n())

nodes_info(conn_nodes_p1h2colors, conn_edges_p1h2colors %>% filter(ConnStr > 200)) %>% filter(!is.na(SCCId)) %>% group_by(SCCSize) %>% summarise(n())








conn_nodes_p1h1colors <- NetfromColors('~/Downloads/engram_cormatrix_r_r2_p (1).csv', 0.1, 0.2)[[1]]
conn_edges_p1h1colors <- NetfromColors('~/Downloads/engram_cormatrix_r_r2_p (1).csv', 0.1, 0.2)[[2]]

tic()
perc_p1h1colors <- inv_percolation(conn_nodes_p1h1colors, conn_edges_p1h1colors)
toc()


ggplot(data.frame(perc_p1h1colors$pvals), aes(log10(p.value), MeanDegree)) + xlab('log10(CONNECTION STRENGTH)') +
  geom_point() + ggtitle('p1h1 REDUCED Mean Degree') + geom_hline(yintercept = 10, linetype = "dashed", color = 'blue') + 
  geom_vline(xintercept = log10(2500), linetype = "dashed", color = 'red') + 
  geom_vline(xintercept = log10(conn_mat_connstr_rdx['CA1', 'PL']), linetype = "dashed") + 
  annotate('text', x = log10(conn_mat_connstr_rdx['CA1', 'PL']) - .1, y = 60, label = 'CA1 -> PL', angle = 90) + 
  geom_vline(xintercept = log10(conn_mat_connstr_rdx['ENTl', 'PL']), linetype = "dashed") + 
  annotate('text', x = log10(conn_mat_connstr_rdx['ENTl', 'PL']) - .1, y = 60, label = 'ENTl -> PL', angle = 90) +
  geom_vline(xintercept = log10(conn_mat_connstr_rdx['CA3', 'ENTl']), linetype = "dashed") + 
  annotate('text', x = log10(conn_mat_connstr_rdx['CA3', 'ENTl']) + .1, y = 75, label = 'CA3-> ENT', angle = 90) + coord_cartesian(xlim = c(0, NA))

ggplot(data.frame(perc_p1h1colors$pvals), aes(log10(p.value), SCCSize)) + xlab('log10(CONNECTION STRENGTH)') +
  geom_point() + ggtitle('p1h1 REDUCED Biggest SCC') + 
  geom_vline(xintercept = log10(2500), linetype = "dashed", color = 'red') + 
  geom_vline(xintercept = log10(conn_mat_connstr_rdx['CA1', 'PL']), linetype = "dashed") + 
  annotate('text', x = log10(conn_mat_connstr_rdx['CA1', 'PL']) - .1, y = 60, label = 'CA1 -> PL', angle = 90) + 
  geom_vline(xintercept = log10(conn_mat_connstr_rdx['ENTl', 'PL']), linetype = "dashed") + 
  annotate('text', x = log10(conn_mat_connstr_rdx['ENTl', 'PL']) - .1, y = 40, label = 'ENTl -> PL', angle = 90) +
  geom_vline(xintercept = log10(conn_mat_connstr_rdx['CA3', 'ENTl']), linetype = "dashed") + 
  annotate('text', x = log10(conn_mat_connstr_rdx['CA3', 'ENTl']) + .1, y = 95, label = 'CA3-> ENT', angle = 90) + coord_cartesian(xlim = c(0, NA))

ggplot(data.frame(perc_p1h1colors$pvals), aes(log10(p.value), WeakSize)) + xlab('log10(CONNECTION STRENGTH)') +
  geom_point() + ggtitle('p1h1 REDUCED WeakCC') +
  geom_vline(xintercept = log10(2500), linetype = "dashed", color = 'red') + 
  geom_vline(xintercept = log10(conn_mat_connstr_rdx['CA1', 'PL']), linetype = "dashed") + 
  annotate('text', x = log10(conn_mat_connstr_rdx['CA1', 'PL']) - .2, y = 60, label = 'CA1 -> PL', angle = 90) + 
  geom_vline(xintercept = log10(conn_mat_connstr_rdx['ENTl', 'PL']), linetype = "dashed") + 
  annotate('text', x = log10(conn_mat_connstr_rdx['ENTl', 'PL']) - .2, y = 60, label = 'ENTl -> PL', angle = 90) + coord_cartesian(xlim = c(0, NA))



nodes_info(conn_nodes_p1h1colors, conn_edges_p1h1colors %>% filter(ConnStr > 100)) %>% filter(!is.na(SCCId)) %>% group_by(SCCSize) %>% summarise(n())

nodes_info(conn_nodes_p1h1colors, conn_edges_p1h1colors %>% filter(ConnStr > 200)) %>% filter(!is.na(SCCId)) %>% group_by(SCCSize) %>% summarise(n())



#### story so far in terms of correlation matrices ####

##### just corr, nothing fancy ####
corr_raq <- read.csv("~/Downloads/engram_cormatrix_r_r2_p (1).csv", stringsAsFactors = F, header = T)
corr_raq <- dcast(corr_raq, row ~ column, value.var = "cor")
rownames(corr_raq) <- t(corr_raq[1])
corr_raq <- corr_raq[-1]
corr_raq <- as.matrix(corr_raq)

corrplot(corr_raq, method = 'color', order = 'hclust', tl.srt = 45, tl.cex = 0.5)
corrRect.hclust(as.matrix(corr_raq), k = 18)

hc <- hclust(as.dist(1-corr_raq))
plot(hc, cex = .6)
abline(h = hc$height[hc$order == 18], col = 'red', lty = 2)

#as.dendrogram(hclust(as.dist((1-corr_raq)/2))) %>% set("labels_cex", 0.6) %>% set("labels_col", "black") %>% plot(horiz = TRUE)


##### pval filter (initially) ####
corr_raq_pfilt <- read.csv("~/Downloads/engram_cormatrix_r_r2_p (1).csv", stringsAsFactors = F, header = T)
corr_raq_pfilt$cor[corr_raq_pfilt$p > .1] <- 0
corr_raq_pfilt <- dcast(corr_raq_pfilt, row ~ column, value.var = "cor")
rownames(corr_raq_pfilt) <- t(corr_raq_pfilt[1])
corr_raq_pfilt <- corr_raq_pfilt[-1]
corr_raq_pfilt <- as.matrix(corr_raq_pfilt)

corrplot(corr_raq_pfilt[hc$order, hc$order], method = 'color', order = 'original', tl.srt = 45, tl.cex = 0.5)
corrRect.hclust(corr_raq, k = 18)

corrplot(corr_raq_pfilt, method = 'color', order = 'hclust', tl.srt = 45, tl.cex = 0.5)
hc_pfilt <- hclust(as.dist(1-corr_raq_pfilt))
corrRect.hclust(corr_raq_pfilt, k = length(unique(cutree(hc_pfilt, h = 0.1))))

plot(hc_pfilt, cex = .6)
abline(h = 0.1, col = 'red', lty = 2)


##### corr filtration  ####
corr_raq_filt <- read.csv("~/Downloads/engram_cormatrix_r_r2_p (1).csv", stringsAsFactors = F, header = T)
corr_raq_filt <- dcast(corr_raq_filt, row ~ column, value.var = "cor")
rownames(corr_raq_filt) <- t(corr_raq_filt[1])
corr_raq_filt <- corr_raq_filt[-1]
corr_raq_filt <- as.matrix(corr_raq_filt)
corr_raq_filt[abs(corr_raq) < 0.6] <- 0

corr_raq_filt <- corr_raq
corr_raq_filt[abs(corr_raq) < (1-sort(hc_abs$height, decreasing = TRUE)[17])] <- 0

#corrplot(corr_raq_filt[hc$order, hc$order], method = 'color', order = 'original', tl.srt = 45, tl.cex = 0.5)
#corrRect.hclust(corr_raq, k = 18)

corrplot(corr_raq_filt, method = 'color', order = 'hclust', tl.srt = 45, tl.cex = 0.5)
corrRect.hclust(corr_raq_filt, k = 27)


# makse's choice: 18 original clusters 
# its just filtering with entire corr matrix, the point of making the correlation filtering is so that the clusters change and then some nodes get lost
nodes_corr_18orig <- data.frame(acronym = rownames(corr_raq_filt))
nodes_corr_18orig$hier_clus <- NA

hier_clus <- data.frame(clusters = cutree(hc, k=18))
for (i in 1:nrow(nodes_corr_18orig)) {
  cluster <- hier_clus$clusters[rownames(hier_clus) == nodes_corr_18orig$acronym[i]]
  if(length(cluster)){
    nodes_corr_18orig$hier_clus[i] <- cluster
  }
}
nodes_corr_18orig <- nodes_corr_18orig %>% group_by(hier_clus) %>% mutate(ClusterSize = n()) %>% filter(ClusterSize > 1)

edges_corr_18orig <- conn_edges_connstr_rdx %>% filter(Source %in% nodes_corr_18orig$acronym & Target %in% nodes_corr_18orig$acronym)
#edges_corr_18orig <- conn_edges_corr_18orig_connstr_rdx %>% filter(Source %in% nodes_corr_18orig$acronym & Target %in% nodes_corr_18orig$acronym)
#edges_corr_18orig <- edges_corr_18orig %>% filter(ConnStr >= 232)

nodes_corr_18orig <- nodes_info(nodes_corr_18orig, edges_corr_18orig)
nodes_corr_18orig <- nodes_corr_18orig[!is.na(nodes_corr_18orig$SCCId),]

write.table(nodes_corr_18orig, file = 'corr18orig_colors.txt', quote = F, row.names = F, col.names = T, sep = '\t')
write.table(edges_corr_18orig %>% filter(ConnStr >= 232), file = 'corr18orig_232ConnStrgraph.txt', quote = F, row.names = F, col.names = F, sep = '\t')
write.table(edges_corr_18orig %>% filter(ConnStr >= 606), file = 'corr18orig_606ConnStrgraph.txt', quote = F, row.names = F, col.names = F, sep = '\t')
write.table(edges_corr_18orig %>% filter(ConnStr >= 2500), file = 'corr18orig_2500ConnStrgraph.txt', quote = F, row.names = F, col.names = F, sep = '\t')


  #only connected nodes
nodes_corr_18orig <- nodes_info(nodes_corr_18orig, edges_corr_18orig %>% filter(ConnStr >= 232))
write.table(nodes_corr_18orig[!is.na(nodes_corr_18orig$SCCId),1:2], file = 'corr18_232_colors.txt', quote = F, row.names = F, col.names = F, sep = '\t')




# dealers choice: 27 new clusters
nodes_corr.6_27 <- data.frame(acronym = rownames(corr_raq_filt))
nodes_corr.6_27$hier_clus <- NA

hier_clus <- data.frame(clusters = cutree(hclust(as.dist(1-corr_raq_filt)), k=27))
for (i in 1:nrow(nodes_corr.6_27)) {
  cluster <- hier_clus$clusters[rownames(hier_clus) == nodes_corr.6_27$acronym[i]]
  if(length(cluster)){
    nodes_corr.6_27$hier_clus[i] <- cluster
  }
}
nodes_corr.6_27 <- nodes_corr.6_27 %>% group_by(hier_clus) %>% mutate(ClusterSize = n()) %>% filter(ClusterSize > 1)

edges_corr.6_27 <- conn_edges_connstr_rdx %>% filter(Source %in% nodes_corr.6_27$acronym & Target %in% nodes_corr.6_27$acronym)
#edges_corr.6_27 <- conn_edges_corr.6_27_connstr_rdx %>% filter(Source %in% nodes_corr.6_27$acronym & Target %in% nodes_corr.6_27$acronym)
edges_corr.6_27 <- edges_corr.6_27 %>% filter(ConnStr >= 232)

nodes_corr.6_27 <- nodes_info(nodes_corr.6_27, edges_corr.6_27)
nodes_corr.6_27 <- nodes_corr.6_27[!is.na(nodes_corr.6_27$SCCId),]

write.table(nodes_corr.6_27, file = 'corr6_27_colors.txt', quote = F, row.names = F, col.names = T, sep = '\t')
write.table(edges_corr.6_27 %>% filter(ConnStr >= 232), file = 'corr6_27_232ConnStr_graph.txt', quote = F, row.names = F, col.names = F, sep = '\t')
write.table(edges_corr.6_27 %>% filter(ConnStr >= 606), file = 'corr6_27_606ConnStr_graph.txt', quote = F, row.names = F, col.names = F, sep = '\t')
write.table(edges_corr.6_27 %>% filter(ConnStr >= 2500), file = 'corr6_27_2500ConnStr_graph.txt', quote = F, row.names = F, col.names = F, sep = '\t')



  #only connected nodes
nodes_corr.6_27 <- nodes_info(nodes_corr.6_27, edges_corr.6_27 %>% filter(ConnStr >= 232))
write.table(nodes_corr.6_27[!is.na(nodes_corr.6_27$SCCId),1:2], file = 'corr6_27_232_colors.txt', quote = F, row.names = F, col.names = F, sep = '\t')







#### Perc from corr ####

# hc <- hclust(as.dist(1-corr))
# 
# hier_clus <- data.frame(clusters = cutree(hc, k=clusters))
# 
# conn_nodes$hier_clus <- NA 
# for (i in 1:nrow(conn_nodes)) {
#   cluster <- hier_clus$clusters[conn_nodes$acronym[i] == rownames(hier_clus)]
#   if(length(cluster)){
#     conn_nodes$hier_clus[i] <- cluster
#   }
# }
# 
# 
# nodes <- get_hier_clus(h, dendrogram, conn_nodes_connstr_rdxCutoff) %>% group_by(hier_clus) %>% mutate(ClusterSize = n())
# 
# 
# nodes <- getColors(file, pval, h)[,c(1,14,15)] %>% filter(ClusterSize > 1)


CorrPerc <- function(file, NumClusters, ConnStrCutoff){
  corr <- read.csv(file, stringsAsFactors = F, header = T)
  corr <- dcast(corr, row ~ column, value.var = "cor")
  rownames(corr) <- t(corr[1])
  corr <- corr[-1]
  corr <- as.matrix(corr)
  
  nodes <- data.frame(acronym = rownames(corr))
  nodes$hier_clus <- NA
  edges <- conn_edges_connstr_rdx %>% filter(Source %in% nodes$acronym & Target %in% nodes$acronym)
  edges <- edges %>% filter(ConnStr >= ConnStrCutoff)
  
  cor_cut <- 0
  cor_cuts <- NULL
  while (cor_cut < 1) {
    print(cor_cut)
    corr[abs(corr) < cor_cut] <- 0
    
    hier_clus <- data.frame(clusters = cutree(hclust(as.dist(1-corr)), k=NumClusters))
    for (i in 1:nrow(nodes)) {
      cluster <- hier_clus$clusters[rownames(hier_clus) == nodes$acronym[i]]
      if(length(cluster)){
        nodes$hier_clus[i] <- cluster
      }
    }
    nodes <- nodes %>% group_by(hier_clus) %>% mutate(ClusterSize = n()) %>% filter(ClusterSize > 1)
    edges <- edges %>% filter(Source %in% nodes$acronym & Target %in% nodes$acronym)
    
    nodes <- nodes_info(nodes, edges)
    nodes <- nodes[!is.na(nodes$SCCId),]
    
    cor_cuts <- rbind(cor_cuts, c(cor_cut, max(unique(nodes$SCCSize)), max(unique(nodes$WeakSize)), mean(nodes$OutDegree, na.rm = T)))
    
    cor_cut <- cor_cut + .01
  }
  
  colnames(cor_cuts) <- c('CorrValue', 'SCCSize', 'WeakSize', 'MeanDegree')
  return(cor_cuts)
}

corr_perc <- CorrPerc("~/Downloads/engram_cormatrix_r_r2_p (1).csv", 18, 232)

ggplot(data.frame(corr_perc), aes(CorrValue, SCCSize)) + xlab('Correlation Values') + geom_point() + ggtitle('Correlation Percolation SCC Size') 




#### 

sum(!is.na(nodes_info(nodes_corr.6_27, edges_corr.6_27 %>% filter(ConnStr >= 232))$SCCId))
mean(nodes_info(nodes_corr.6_27, edges_corr.6_27 %>% filter(ConnStr >= 2500))$OutDegree, na.rm = T)
nodes_info(nodes_corr.6_27, edges_corr.6_27 %>% filter(ConnStr >= 2500)) %>% filter(!is.na(SCCId)) %>% group_by(SCCSize) %>% summarise(n())



sum(!is.na(nodes_info(nodes_corr_18orig, 
                      read.delim("davids_output_18connected/corr18orig_232ConnStrgraph.txt_corr18_232_colors.txt.ar.balanced.directed.out.graph.txt", sep = ' ', header = F))$SCCId))
mean(nodes_info(nodes_corr_18orig, 
                      read.delim("davids_output_18connected/corr18orig_232ConnStrgraph.txt_corr18_232_colors.txt.ar.balanced.directed.out.graph.txt", sep = ' ', header = F))$OutDegree, na.rm = T)
nodes_info(nodes_corr_18orig, 
           read.delim("davids_output_18connected/corr18orig_232ConnStrgraph.txt_corr18_232_colors.txt.ar.balanced.directed.out.graph.txt", sep = ' ', header = F)) %>% 
  filter(!is.na(SCCId)) %>% group_by(SCCSize) %>% summarise(n())




#### Post-david ####
addStrNAdded <- function(davids_output, colors_file, corr = corr_raq, output_name = '', rtrn = F, abs = T){
  output <- read.delim(davids_output, sep = ' ', header = F)
  print('read graph')
  colnames(output) <- c('Source', 'Target')
  output$ConnStr <- NA
  for (i in 1:nrow(output)) {
#    if(!i%%100){print(i)} 
    output$ConnStr[i] <- conn_edges_connstr_rdx$ConnStr[conn_edges_connstr_rdx$Source == output$Source[i] & conn_edges_connstr_rdx$Target == output$Target[i]]
  }
#  output$Added <- apply(output[1:2], 1, function(row) all(apply(edges_corr_18orig[edges_corr_18orig$ConnStr > Cutoff,1:2], 1, function(row2) !all(row == row2))))
#  print('added')
  
  nodes <- read.delim(colors_file, header = F)
  print('read nodes')
  colnames(nodes) <- c('acronym', 'Color')
  if(abs){
    hc <- hclust(as.dist(1-abs(corr)))
    nodes$SubColor <- NA
    for (i in unique(nodes$Color)) {
      nodesInColor <- nodes$acronym[nodes$Color == i]
      subcor <- corr[which(cutree(hc, length(unique(nodes$Color))) == i), which(cutree(hc, length(unique(nodes$Color))) == i)]
      if(length(nodesInColor) > 1 & any(subcor < 0)){
        hier_clus <- data.frame(clusters = cutree(hclust(as.dist(1-subcor)), k=2))
        print(paste(i, nodesInColor))
        for (j in 1:length(nodesInColor)) {
          print(nodesInColor[j])
          subcolor <- hier_clus$clusters[rownames(hier_clus) == nodesInColor[j]]
          nodes$SubColor[nodes$acronym == nodesInColor[j]] <- subcolor
        }
      }
    }
  }
  nodes <- nodes %>% group_by(Color) %>% mutate(ClusterSize = n())
  nodes <- nodes_info(nodes, output)
  
  print(paste('edges:', nrow(output)))
  print(paste('nodes:', sum(!is.na(nodes$SCCId))))
  print(paste('Mean Degree:', mean(nodes$OutDegree, na.rm = T)))
  print('SCCs')
  print(nodes %>% filter(!is.na(SCCId)) %>% group_by(SCCSize) %>% summarise(n()))
  
  
  if(nchar(output_name)){
    write.table(output, file = paste(output_name, '_graph.csv', sep = ''), quote = F, row.names = F, col.names = F, sep = '\t')
    write.table(nodes[!is.na(nodes$SCCId),], file = paste(output_name, '_nodes.csv', sep = ''), quote = F, row.names = F, col.names = T, sep = '\t')
  }
  if(rtrn){
    list <- NULL
    list$nodes <- nodes[!is.na(nodes$SCCId),]
    list$edges <- output
    return(list)
  }
}


#addStrNAdded("davids_output/abs232graph.txt_abs232colors.txt.ar.balanced.directed.out.graph.txt", 232, 'davids_output/abs232colors.txt')
#addStrNAdded("davids_output/retr232graph.txt_retr232colors.txt.ar.balanced.directed.out.graph.txt", 'davids/retr232colors.txt', 'retr232', corr_raq_retrieval, T)
abs232_25colors <- addStrNAdded("davids_output/abs232graph.txt_abs232_25colors.txt.ar.balanced.directed.out.graph.txt", 'abs(corr)/abs232_25colors.txt', rtrn = T)
abs232_25colors_edges <- abs232_25colors$edges
abs232_25colors_nodes <- abs232_25colors$nodes

  
addCollapsed <- function(nodes_file, edges_file, output_file, return = F){
  nodes <- read.delim(nodes_file, header = T)
  edges <- read.delim(edges_file, header = F, sep = '\t')
  
  col <- collapseColors(nodes, edges)
  
  write.table(col[[2]], file = paste(output_file, 'edges.csv', sep = '_'), quote = F, row.names = F, col.names = F, sep = '\t')
  write.table(col[[1]], file = paste(output_file, 'nodes.csv', sep = '_'), quote = F, row.names = F, col.names = T, sep = '\t')
  
  if(return){
    list <- NULL
    list$nodes <- col[[1]]
    list$edges <- col[[2]]
    
    return(list)
  }
}


Comparing <- function(dir, toCompare = abs232_edges){
  print(paste('ref has', nrow(toCompare), 'edges'))
  graph_files <- grep("graph\\.txt$", list.files(dir, full.names = TRUE), value = TRUE)
  
  for (graph in graph_files) {
#    print(graph)
    
    output <- read.delim(graph, sep = ' ', header = F)
#    print('read graph')
    colnames(output) <- c('Source', 'Target')
    # output$ConnStr <- NA
    # for (i in 1:nrow(output)) {
    #   #    if(!i%%100){print(i)} 
    #   output$ConnStr[i] <- conn_edges_connstr_rdx$ConnStr[conn_edges_connstr_rdx$Source == output$Source[i] & conn_edges_connstr_rdx$Target == output$Target[i]]
    # }
    
    print(paste(gsub("^[^0-9]*([0-9]+).*", "\\1", graph), 'shared', nrow(inner_join(output[1:2], toCompare[1:2]))))
#    print(anti_join(toCompare[1:2], output[1:2]))
    
  }
}

#### Dealing with negative Correlations ####

corrplot(abs(corr_raq), method = 'color', order = 'hclust', tl.srt = 45, tl.cex = 0.5)
corrRect.hclust(as.matrix(abs(corr_raq)), k = 18)

PreDavid<- function(correlation = corr_raq, colors, file_name, CutOffs = c(232, 606, 2500)){
  nodes <- data.frame(acronym = rownames(correlation))
  nodes$hier_clus <- NA
  
  hc <- hclust(as.dist(1-correlation))
  hier_clus <- data.frame(clusters = cutree(hc, k=colors))
  for (i in 1:nrow(nodes)) {
    cluster <- hier_clus$clusters[rownames(hier_clus) == nodes$acronym[i]]
    if(length(cluster)){
      nodes$hier_clus[i] <- cluster
    }
  }
  
  nodes <- nodes %>% group_by(hier_clus) %>% mutate(ClusterSize = n())# %>% filter(ClusterSize > 1)

  edges <- conn_edges_connstr_rdx %>% filter(Source %in% nodes$acronym & Target %in% nodes$acronym)
  
  for (cutoff in CutOffs) {
    write.table(edges %>% filter(ConnStr >= cutoff), file = paste(file_name, cutoff, 'graph.txt', sep = ''), quote = F, row.names = F, col.names = F, sep = '\t')
    nodes <- nodes_info(nodes, edges %>% filter(ConnStr >= cutoff))
    write.table(nodes[!is.na(nodes$SCCId),1:2], file = paste(file_name, cutoff, 'colors.txt', sep = ''), quote = F, row.names = F, col.names = F, sep = '\t')
  }
  # write.table(edges %>% filter(ConnStr >= 606), file = paste(file_name, '606graph.txt', sep = ''), quote = F, row.names = F, col.names = F, sep = '\t')
  # write.table(edges %>% filter(ConnStr >= 2500), file = paste(file_name, '2500graph.txt', sep = ''), quote = F, row.names = F, col.names = F, sep = '\t')
  
  # nodes <- nodes_info(nodes, edges %>% filter(ConnStr >= 606))
  # write.table(nodes[!is.na(nodes$SCCId),1:2], file = paste(file_name, '606colors.txt', sep = ''), quote = F, row.names = F, col.names = F, sep = '\t')
  # nodes <- nodes_info(nodes, edges %>% filter(ConnStr >= 2500))
  # write.table(nodes[!is.na(nodes$SCCId),1:2], file = paste(file_name, '2500colors.txt', sep = ''), quote = F, row.names = F, col.names = F, sep = '\t')
}







addStrNAdded("davids_output/abs232graph.txt_abs232colors.txt.ar.balanced.directed.out.graph.txt", 232)

#### CI ####

CollectiveInfluence <- function(nodes, edges, l){
  graph <- graph_from_edgelist(as.matrix(edges[, 1:2]), directed = TRUE)
  for (i in 1:nrow(nodes)) {
    print(nodes$acronym[i])
    neighbors <- neighborhood(graph, order = l, nodes = nodes$acronym[i], mode = 'out', mindist = l-1)
    nodes$CI[i] <- degree(graph, nodes$acronym[i], mode = 'out') * sum(degree(graph, unlist(neighbors), mode = 'out') - 1)
    
  }
  return(nodes)
}

#### Perc from Topology ####

ConnPerc <- function(nodes, edges, prohibitted = c()){
  nodes <- nodes_info(nodes, edges)
  cols <- c("Betweenness", "Closeness", "InCloseness", "OutCloseness", "CI")
  
  q <- 0
#  perc <- c(q, max(nodes$WeakSize))
  max <- max(nodes$WeakSize)
  perc <- NULL
  for (Column in cols) {
    print(Column)
    col <- NULL
    nodes_perc <- nodes
    edges_perc <- edges
    
    nodes_perc <- nodes_perc[!is.na(nodes_perc$SCCId),]
#    col <- max(nodes_perc$WeakSize)
    while (nrow(nodes_perc) > length(prohibitted)) {
      #    q <- max(nodes_perc[Column], na.rm = T)
      q <- q+1
#      print(q)
      nodes_perc <- nodes_info(nodes_perc, edges_perc)
#      nodes_perc <- nodes_perc[!is.na(nodes_perc$SCCId),]
      
      
      nodes_perc <- nodes_perc[order(nodes_perc[[Column]], decreasing = T),]
      j <- 1
      while (nodes_perc$acronym[j] %in% prohibitted){
        j <- j+1
        print(paste('j=',j))
        }
      col <- rbind(col, c(nodes_perc$acronym[j], max(nodes_perc$WeakSize)))
      nodes_perc <- nodes_perc[-j,]
      #    nodes_perc <- nodes_perc %>% filter(Column < max(Column)) #filtrate more than once node 
      edges_perc <- edges_perc %>% filter(Source %in% nodes_perc$acronym & Target %in% nodes_perc$acronym)
      
      nodes_perc <- nodes_perc %>% filter(acronym %in% unique(c(edges_perc$Source, edges_perc$Target)))
      print(paste('nodes:', nrow(nodes_perc)))
    }
#    perc <- cbind(perc, col)
    perc[[Column]] <- col
  }
  return(perc)
}


nodes <- abs232_nodes
edges <- abs232_edges


RmROIs <- function(ToRm, nodes = abs232_nodes, edges = abs232_edges, rtrn = F, prnt = T){
  nodes_rm <- nodes %>% filter(!acronym %in% ToRm)
  edges_rm <- edges %>% filter(Source %in% nodes_rm$acronym & Target %in% nodes_rm$acronym)
  
  nodes_rm <- nodes_info(nodes_rm, edges_rm)
  
  if(sum(is.na(nodes_rm$SCCId))){
    if(prnt){
      print('Following nodes get disconneted (original stats shown):')
      disc <- nodes_rm$acronym[is.na(nodes_rm$SCCId)]
      print(nodes %>% filter(acronym %in% disc))
    }
    nodes_rm <- nodes_rm[!is.na(nodes_rm$SCCId),]
  }
  
  if(prnt){
    if(nrow(nodes %>% group_by(SCCSize) %>% summarise()) != nrow(nodes_rm %>% group_by(SCCSize) %>% summarise())){
      print('# of SCCs have changed!')
      print('Before')
      print(unlist(nodes %>% group_by(SCCSize) %>% summarise()))
      print('After')
      print(unlist(nodes_rm %>% group_by(SCCSize) %>% summarise()))
    } else {
      if(any(nodes %>% group_by(SCCSize) %>% summarise() != nodes_rm %>% group_by(SCCSize) %>% summarise())){
        print('Size of SCCs have changed!')
        print('Before')
        print(unlist(nodes %>% group_by(SCCSize) %>% summarise()))
        print('After')
        print(unlist(nodes_rm %>% group_by(SCCSize) %>% summarise()))
      }
    }
  }
  
  if(rtrn){
    list <- NULL
    list$nodes <- nodes_rm
    list$edges <- edges_rm %>% filter(Source %in% nodes_rm$acronym & Target %in% nodes_rm$acronym)
    return(list)
  }
}

RmROIList <- function(List, nodes = abs232_nodes, edges = abs232_edges, biggest_group = 3){
  SCCs <- unlist(nodes_info(nodes, edges) %>% group_by(SCCSize) %>% summarise())
  
  list <- NULL
  for (i in 1:biggest_group) {
    groups <- combn(List, i)
    print(paste('going into groups of', i))   
    
    for (j in 1:ncol(groups)) {
      nodes_rm <- RmROIs(groups[,j], nodes, edges, T, F)$nodes
      
      sccs <- unlist(nodes_rm %>% group_by(SCCSize) %>% summarise())
      size <- sccs[!sccs %in% SCCs]
      
      if(length(SCCs)-length(sccs) > 0){
        size <- 0
#        print(c(node, length(sccs)-1, size))
      } 
      if(length(SCCs)-length(sccs) < 0){
        print(paste(groups[,j], 'splits SCC'))
      }
      
      list <- rbind(list, c(paste(groups[,j], collapse = ', '), length(groups[,j]), length(sccs)-1, as.integer(size)))
      
    }
  }


  colnames(list) <- c('RmROI(s)', '#ROIs', '#SCCs', 'MaxSCCSize')
  return(data.frame('RmROIs' = list[,1], 'NumbROIs' = list[,2], 'NumbSCCs' = list[,3], 'MaxSCCSize' = as.integer(list[,4])))
}

tic()
rms <- RmROIList(abs232_nodes$acronym[abs232_nodes$SCCSize == 27])
toc()
tic()
rms_retr <- RmROIList(abs232_retr_nodes$acronym[abs232_retr_nodes$SCCSize == 44], abs232_retr_nodes, abs232_retr_edges)
toc()

#### dendrogram height #### 

#the height or threshold cut here is not directly on the corr matrix but on the distance generated by 1-corr

plot(hc_abs, cex = .6)
rect(par("usr")[1], sort(hc_abs$height, decreasing = TRUE)[18], par("usr")[2], sort(hc_abs$height, decreasing = TRUE)[17], 
     col = "red", border = NA, density = 50)
abline(h = 0.4, col = 'red', lty = 2)

# in order of hclust for abs(correlations) to see the neg corr in clusters from abs
corrplot(corr_raq[hc_abs$order,hc_abs$order], method = 'color', order = 'original', tl.srt = 45, tl.cex = 0.5)
corrRect.hclust(as.matrix(abs(corr_raq)), k = 18)

#### running through clusters separetely ####

# subcolors for neg correlations
abs232_nodes$SubColor <- NA
for (i in unique(abs232_nodes$Color)) {
  nodesInColor <- abs232_nodes$acronym[abs232_nodes$Color == i]
  subcor <- corr_raq[which(cutree(hc_abs, 18) == i), which(cutree(hc_abs, 18) == i)]
  if(length(nodesInColor) > 1 & any(subcor < 0)){
    hier_clus <- data.frame(clusters = cutree(hclust(as.dist(1-subcor)), k=2))
    for (j in 1:length(nodesInColor)) {
      subcolor <- hier_clus$clusters[rownames(hier_clus) == nodesInColor[j]]
      abs232_nodes$SubColor[abs232_nodes$acronym == nodesInColor[j]] <- subcolor
    }
    
  }
}

# finding min in clusters and max outside
clusters <- matrix(0, nrow = nrow(corr_raq), ncol = nrow(corr_raq), dimnames = list(rownames(corr_raq), colnames(corr_raq))) 
for (i in unique(abs232_nodes$Color)) {
#  nodesInColor <- abs232_nodes$acronym[abs232_nodes$Color == i]
  clusters[which(cutree(hc_abs, 18) == i), which(cutree(hc_abs, 18) == i)] <- corr_raq[which(cutree(hc_abs, 18) == i), which(cutree(hc_abs, 18) == i)]
}
outside <- corr_raq[hc_abs$order,hc_abs$order]-clusters[hc_abs$order,hc_abs$order]

corrplot(clusters[hc_abs$order,hc_abs$order], method = 'color', order = 'original', tl.srt = 45, tl.cex = 0.5)
corrplot(outside, method = 'color', order = 'original', tl.srt = 45, tl.cex = 0.5)

outside_filt <- outside
outside_filt[abs(outside) < 0.6] <- 0
corrplot(outside_filt, method = 'color', order = 'original', tl.srt = 45, tl.cex = 0.5)
corrRect.hclust(as.matrix(abs(corr_raq)), k = 18)

clusters_filt <- clusters
clusters_filt[abs(clusters) > 0.6] <- 0
corrplot(clusters_filt[hc_abs$order,hc_abs$order], method = 'color', order = 'original', tl.srt = 45, tl.cex = 0.5)
corrRect.hclust(as.matrix(abs(corr_raq)), k = 18)


corr_raq_filt <- corr_raq
corr_raq_filt[abs(corr_raq) < (1-sort(hc_abs$height, decreasing = TRUE)[18])] <- 0
corr_raq_filt[abs(corr_raq) <= 0.67] <- 0
corrplot(corr_raq_filt[hc_abs$order,hc_abs$order], method = 'color', order = 'original', tl.srt = 45, tl.cex = 0.5)
corrRect.hclust(as.matrix(abs(corr_raq)), k = 18)



#### collapsing #### 
collapseColors <- function(nodes, edges){
  nodes_col <- c() 
  for(i in unlist(nodes %>% group_by(Color) %>% summarise())){
    fiber <- nodes %>% filter(Color == i)
    SCC <- fiber %>% filter(SCCSize > 1)
    if(nrow(SCC) > 1){
      print(paste("SCCs in fiber ", i, ": ", SCC$acronym))
    }
    
    if(nrow(SCC)){
      nodes_col <- rbind(nodes_col, SCC[1,])
    } else {
      nodes_col <- rbind(nodes_col, fiber[1,])
    }
  }
  #nodes_col <- nodes_col[complete.cases(nodes_col),]
  
  colnames(edges) <- c('Source', 'Target')
  edges_col <- edges
  edges_col <- edges_col[edges_col$Target %in% nodes_col$acronym,]
  for(i in 1:nrow(edges_col)){
    print(paste(i, edges_col$Source[i]))
    if(!edges_col$Source[i] %in% nodes_col$acronym){
      edges_col$Source[i] <- nodes_col$acronym[ nodes_col$Color == nodes$Color[ nodes$acronym == edges_col$Source[i] ] ]
    }
  }
  print('aqui')
  
  nodes_col <- nodes_info(nodes_col, edges_col) 
  
  list <- NULL
  list$nodes <- nodes_col
  list$edges <- edges_col
  return(list)
}




#### Cycles ####

getcycles_fromSCC <- function(graph, nodes, sccid){
  sccg <- induced.subgraph(graph, nodes)
  ids <- 1:vcount(sccg)
  Cycles <- NULL
  
  print(paste(vcount(sccg), 'nodes in SCC'))
  
  #Cycles starting in every node in ascending order of vertex id , eliminating each node as taken into account
  #  Cycles <- NULL
  while(vcount(sccg) != 1){
    #print(paste('nodes = ', vcount(sccg)))
    
    if(vcount(sccg) > 20 & vcount(sccg) %% 5 == 0){print(paste(vcount(sccg), 'nodes left in the SCC'))}
    
    #All cycles of smallest (vertexid) node 
    nodeid <- ids[1]
    
    neighbors <- neighbors(sccg, nodeid, mode = 'out')
    if(length(neighbors) > 10){print(paste(length(neighbors), 'neighbors'))}
    
    temp_cycles <- NULL
    for (j in neighbors) {
      #print(paste('node =', get.vertex.attribute(sccg, "name", nodeid), 'neighbor = ', get.vertex.attribute(sccg, "name", j)))
      temp_cycles <- c(temp_cycles, lapply(all_simple_paths(sccg, j, nodeid, mode="out"), function(p) c(nodeid,p)))
      #print(paste('temp_cycles =', length(temp_cycles)))
    }
    
    #add the new cycles to list of all cycles
    cycles <- NULL
    if(length(temp_cycles) != 0){
      for (j in 1:length(temp_cycles)) {
        cycles$SCCId <- sccid
        cycles$Length[j] <- length(temp_cycles[[j]]) - 1
        cycles$Path[j] <- paste(data.frame(c = get.vertex.attribute(sccg, "name", temp_cycles[[j]]))$c, collapse = ' ')
      }
    }
    Cycles <- rbind(Cycles, data.frame(cycles, stringsAsFactors = F))
    #print(paste('Cycles =', nrow(Cycles)))
    
    #Remove node i just taken into account
    sccg <- induced.subgraph(sccg, ids[-1])
    ids <- 1:vcount(sccg)
  }
  return(Cycles)
}

getcycles_fromEdges <- function(edges){
  #edges <- edges %>% filter(Source != Target)
  nodes <- as.data.frame(sort(unique(c(edges$Source, edges$Target))))
  colnames(nodes)[1] <- "Label"
  g <- graph_from_data_frame(d = edges, vertices = nodes, directed = 1)
  
  sc <- clusters(g, mode = "strong")
  if(length(sc$csize)){
    nodes$SCCId <- sc$membership
    SCCids <- NULL
    for (i in 1:length(sc$csize)) {
      if(sc$csize[i] > 1){SCCids <- c(SCCids, i)}
    }
    
    Cycles <- NULL
    #Find cycles in each SCC 
    for (sccid in SCCids){
      Cycles <- rbind(Cycles, getcycles_fromSCC(g, nodes$SCCId == sccid, sccid))
    }
    Cycles <- Cycles[!duplicated(Cycles),]
    return(Cycles)
  } else {
    return(NULL)
  }
  
}




#### percolation of corr matrix as adjacency ####

perc_corr_raq <- inv_percolation(data.frame(acronym = rownames(corr_raq)), read.csv("~/Downloads/engram_cormatrix_r_r2_p (1).csv", stringsAsFactors = F, header = T), T)
ggplot(data.frame(perc_corr_raq$pvals), aes((p.value), SCCSize)) + geom_point() + xlab('abs(Correlation)') + ylab('GCCSize') +
  geom_vline(xintercept = 0.8246, linetype = "dashed", color = 'blue') + 
  annotate('text', x = 0.8246 + .02, y = 60, label = 'Cluster #15', angle = 90) +
  geom_vline(xintercept = 0.7808, linetype = "dashed", color = 'blue') + 
  annotate('text', x = 0.7808 + .02, y = 30, label = 'Cluster #18', angle = 90) 
View(nodes_info(data.frame(acronym = rownames(corr_raq)), edges %>% filter(p_values >= 0.9547774)))


#which clusters gets filled first
while (condition) {
  
}


disconnectedThresholds_18_enc <- DiscThresh(18, )
for (i in 1:18) {
  print(paste(i, ';', min(abs(corr_raq[which(cutree(hc_abs, 18) == i), which(cutree(hc_abs, 18) == i)]))))
}

#### comparando mismos colores en retrieval ####
corr_raq_retrieval <- read.csv("~/Downloads/engram_retrieval_cormatrix_r_r2_p.csv", stringsAsFactors = F, header = T)
corr_raq_retrieval <- dcast(corr_raq_retrieval, row ~ column, value.var = "cor")
rownames(corr_raq_retrieval) <- t(corr_raq_retrieval[1])
corr_raq_retrieval <- corr_raq_retrieval[-1]
corr_raq_retrieval <- as.matrix(corr_raq_retrieval)

#corr_raq_retrieval <- corr_raq_retrieval[hc_abs$order,hc_abs$order]

corrplot(corr_raq_retrieval, method = 'color', order = 'original', tl.srt = 45, tl.cex = 0.5)
corrRect.hclust(as.matrix(abs(corr_raq)), k = 18)

corr_raq_retrieval_filt <- corr_raq_retrieval
corr_raq_retrieval_filt[abs(corr_raq_retrieval) < 0.3] <- 0
corrplot(corr_raq_retrieval_filt[hc_abs$order,hc_abs$order], method = 'color', order = 'original', tl.srt = 45, tl.cex = 0.5)
corrRect.hclust(as.matrix(abs(corr_raq)), k = 18)



#order of matrixes should be the same
SynchDec <- function(NumClusters, Cutoff, corr = corr_raq_retrieval, hc = hc_abs){
  cols <- c()
  for (i in 1:NumClusters) {
#    nodesInColor <- abs232_nodes$acronym[abs232_nodes$Color == i]
    subcor <- corr[which(cutree(hc, NumClusters) == i), which(cutree(hc, NumClusters) == i)]
    if(length(subcor) > 1){
      Dec <- (sum(abs(subcor) > Cutoff)-nrow(subcor))/(length(subcor)-nrow(subcor))
    } else {
      Dec <- NA
    }
#    print(mean(subcor))
    mn <- mean(abs(subcor[abs(subcor) > Cutoff & abs(subcor) < 1]))
    min <- min(abs(subcor))
    cols <- rbind(cols, c(paste(i, ':', rownames(corr)[cutree(hc_abs, NumClusters) == i][1], '(', nrow(subcor), 'nodes)'), Dec, mn, min))
#    cols <- rbind(cols, c(paste(i, ':', '(', ncol(subcor), 'nodes)'), Dec, mn))
  }
  colnames(cols) <- c('Color', 'Cij_remaining', 'Mean', 'FilledAt')
  return(data.frame(cols))
}




perc_corr_raq_retrieval <- inv_percolation(data.frame(acronym = rownames(corr_raq)), read.csv("~/Downloads/engram_retrieval_cormatrix_r_r2_p.csv", stringsAsFactors = F, header = T), T)
ggplot(data.frame(perc_corr_raq_retrieval$pvals), aes((p.value), SCCSize)) + geom_point() + xlab('abs(Correlation)')

#retrieval gets disc at ~0.8.... but I've already got disconnected from (ENCODING) clusters at .4 


DiscThresh <- function(NumClusters, InitialCutoff = 0, corr = corr_raq_retrieval, hc = hc_abs){
  corr_filt <- corr
  corr_filt[abs(corr_filt) < InitialCutoff] <- 0
  thresh <- InitialCutoff
  
  nodes <- data.frame(acronym = rownames(corr), Cluster = NA, ClusterSize = NA, DiscThreshold = NA)
  
  while (thresh < 1) {
    print(thresh)
    for (i in 1:NumClusters) {
      subcor <- corr_filt[which(cutree(hc, NumClusters) == i), which(cutree(hc, NumClusters) == i)]
      if(length(subcor) > 1){
        for (node in rownames(subcor)) {
#          print(paste(node, sum(subcor[node,] != 0)))
          if(sum(subcor[node,] != 0) == 1 & is.na(nodes$DiscThreshold[nodes$acronym == node])){
            nodes$DiscThreshold[nodes$acronym == node] <- thresh
          }
          nodes$Cluster[nodes$acronym == node] <- i
          nodes$ClusterSize[nodes$acronym == node] <- nrow(subcor)
        }
      }
    }
    thresh <-  thresh + 0.01
    corr_filt[abs(corr_filt) < thresh] <- 0
    
  }
  
  return(nodes)
}


for (node in rownames(subcor)) {
  print(nodes[nodes$acronym == node])
}





