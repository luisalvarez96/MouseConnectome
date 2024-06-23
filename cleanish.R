library(tidyr)
library(dplyr)
library(igraph)
library(tictoc)
library(ggplot2)
library(reshape2)
library('dendextend')




####
nodes_info <- function(nodes, edges, l = 2, directed = T){
  #doesn't work if there is a node 0  -> use .txt not _Results.csv
  #edges need to be label -> label, not ids
  graph <- graph_from_edgelist(as.matrix(edges[, 1:2]), directed = directed)
  
  # SCCIds <- data.frame(SCCId = components(graph, mode = "strong")$membership)
  # OutDegree <- data.frame(Out = degree(graph, mode = "out", loops = F, normalized = F))
  # InDegree <- data.frame(In = degree(graph, mode = "in", loops = F, normalized = F))
  # Coreness <- data.frame(Core = coreness(graph, mode = "out"))
  # Betweenness <- data.frame(Btwns = betweenness(graph, directed = T))
  # Closeness <- data.frame(Close = closeness(graph, mode = 'all'))
  # InCloseness <- data.frame(Close = closeness(graph, mode = 'in'))
  # OutCloseness <- data.frame(Close = closeness(graph, mode = 'out'))
  # WeakIds <- data.frame(WeakId = components(graph, mode = "weak")$membership)
  
  SCCIds <- components(graph, mode = "strong")$membership
  OutDegree <- degree(graph, mode = "out", loops = F, normalized = F)
  InDegree <- degree(graph, mode = "in", loops = F, normalized = F)
  Coreness <- coreness(graph, mode = "out")
  Betweenness <- betweenness(graph, directed = T)
  Closeness <- closeness(graph, mode = 'all')
  InCloseness <- closeness(graph, mode = 'in')
  OutCloseness <- closeness(graph, mode = 'out')
  WeakIds <- components(graph, mode = "weak")$membership
  
  for(i in 1:nrow(nodes)){
    # nodes$SCCId[i] <- SCCIds[nodes$acronym[i],]
    # nodes$OutDegree[i] <- OutDegree[nodes$acronym[i],]
    # nodes$InDegree[i] <- InDegree[nodes$acronym[i],]
    # nodes$Coreness[i] <- Coreness[nodes$acronym[i],]

    nodes$SCCId[i] <- SCCIds[nodes$acronym[i]]
    nodes$OutDegree[i] <- OutDegree[nodes$acronym[i]]
    nodes$InDegree[i] <- InDegree[nodes$acronym[i]]
    nodes$Coreness[i] <- Coreness[nodes$acronym[i]]
    
    if(nodes$acronym[i] %in% c(edges$Source, edges$Target)){
      neighbors <- neighborhood(graph, order = l, nodes = nodes$acronym[i], mode = 'out', mindist = l-1)
      nodes$CI[i] <- degree(graph, nodes$acronym[i], mode = 'out') * sum(degree(graph, unlist(neighbors), mode = 'out') - 1)
    } else {
      nodes$CI[i] <- NA
      }
    
    # nodes$Betweenness[i] <- Betweenness[nodes$acronym[i],]
    # nodes$Closeness[i] <- Closeness[nodes$acronym[i],]
    # nodes$InCloseness[i] <- InCloseness[nodes$acronym[i],]
    # nodes$OutCloseness[i] <- OutCloseness[nodes$acronym[i],]
    # nodes$WeakId[i] <- WeakIds[nodes$acronym[i],]
    
    nodes$Betweenness[i] <- Betweenness[nodes$acronym[i]]
    nodes$Closeness[i] <- Closeness[nodes$acronym[i]]
    nodes$InCloseness[i] <- InCloseness[nodes$acronym[i]]
    nodes$OutCloseness[i] <- OutCloseness[nodes$acronym[i]]
    nodes$WeakId[i] <- WeakIds[nodes$acronym[i]]
  }
  
  nodes <- nodes %>%
    group_by(SCCId) %>%
    mutate(SCCSize = n()) %>%
    ungroup() 
  
  nodes <- nodes %>%
    group_by(WeakId) %>%
    mutate(WeakSize = n()) %>%
    ungroup() 
  
  return(nodes)  
}
####

#### New Connectome


#read adjacency-pvalues matrix
conn_mat_newAllen <- read.csv('~/Downloads/normalized_connection_density.csv', stringsAsFactors = F, row.names = 1)
conn_mat_newAllen <- conn_mat_newAllen[grep('ipsi', colnames(conn_mat_newAllen))]
colnames(conn_mat_newAllen) <- conn_mat_newAllen[1,]
conn_mat_newAllen <- conn_mat_newAllen[-1,]
for (col in colnames(conn_mat_newAllen)){
  conn_mat_newAllen[, col] <- as.numeric(conn_mat_newAllen[, col])
}



#### collapse for raquel ####

ROIs <- data.frame('ROIs' = rownames(conn_mat_newAllen))
ROIs$name_allen <- NA
for (i in 1:nrow(ROIs)) {
  allen <- allen_dict$name[ROIs$ROIs[i] == allen_dict$acronym]
  if(length(allen)){ROIs$name_allen[i] <- allen}
}
exceptions <- c('Isocortex', 'sAMY', 'x', 'y', 'Eth', 'CTXsp', 'IntG', 'Pa4', 'Pa5', 'PDTg', 'PeF', 'ProS', 'PoT', 'Xi',
                'RAmb', 'Su3', 'SubG', 'TEa', 'VeCB', 'Acs5', 'ANcr1', 'ANcr2', 'PTLp', 'SSs')

#a mano: 'P-sen', 'P-mot', 'P-sat', PVa PVi PVp PVpo, SC, AMd AMv, MBsen MBmot MBsta, ENTl ENTm, EPd EPv, LSc LSr LSv, ICc ICd ICe, MY-sen MY-mot MY-sat, MOp, MOs, PMdm, PMv, AP, APr
#what with SC PV? that dont have 'full ROI' 
amano <- c('P-sen', 'P-mot', 'P-sat', 'PVa', 'PVi', 'PVp', 'PVpo', 'SCs', 'SCop', 'SCsg', 'SCzo', 'SCm', 'SCdg', 'SCdw', 'SCiw', 'SCig', 'AMd', 'AMv', 'MBsen', 'MBmot', 'MBsta', 
           'ENTl', 'ENTm', 'EPd', 'EPv', 'LSc', 'LSr', 'LSv', 'ICc', 'ICd', 'ICe', 'MY-sen', 'MY-mot', 'MY-sat', 'MOp', 'MOs', 'PMd', 'PMv', 'AP', 'APr',
           'VISal', 'VISam', 'VISl', 'VISp', 'VISpl', 'VISpm', 'VISli', 'VISpor', 'VISa', 'VISrl', 'SNr', 'SNc')


ROIs <- ROIs[!ROIs$ROIs %in% exceptions, ]
ROIs <- ROIs[!ROIs$ROIs %in% amano, ]
ROIs <- ROIs[grep('[a-z]', ROIs$ROIs),]
ROIs$short <- NA
ROIs$NumLayers <- NA
for (i in 1:nrow(ROIs)) {
  short <- gsub('[a-z]|-', '', ROIs$ROIs[i])
  num <- length(grep(short, rownames(conn_mat_newAllen), fixed = T))
  #  num <- length(grep(paste('^', short, sep = ''), rownames(conn_mat_newAllen), fixed = T))
  if(length(num)){
    ROIs$short[i] <- short
    ROIs$NumLayers[i] <- num
  }
}

ROIs$short[grep('SSp', ROIs$ROIs)] <- 'SSp'




Collapse <- function(ROI, conn_mat1){
  
  rows <- conn_mat1[grep(ROI, rownames(conn_mat1)), ]
  cols <- conn_mat1[, grep(ROI, colnames(conn_mat1))]
  newrow <- colSums(rows)
  newcol <- rowSums(cols)
  
  conn_mat1 <- rbind(conn_mat1, newrow)
  conn_mat1 <- cbind(conn_mat1, c(newcol, NA))
  
  conn_mat1 <- conn_mat1[grep(ROI, rownames(conn_mat1), invert = TRUE), ]
  conn_mat1 <- conn_mat1[, grep(ROI, colnames(conn_mat1), invert = TRUE)]
  
  rownames(conn_mat1)[length(conn_mat1)] <- ROI
  colnames(conn_mat1)[length(conn_mat1)] <- ROI
  
  return(conn_mat1) 
}


manual <- function(ROI, ROI_list, conn_mat1){
  
  rows <- conn_mat1[rownames(conn_mat1) %in% ROI_list, ]
  cols <- conn_mat1[, colnames(conn_mat1) %in% ROI_list]
  
  newrow <- colSums(rows)
  newcol <- rowSums(cols)
  
  conn_mat1 <- rbind(conn_mat1, newrow)
  conn_mat1 <- cbind(conn_mat1, c(newcol, NA))
  
  conn_mat1 <- conn_mat1[!rownames(conn_mat1) %in% ROI_list, ]
  conn_mat1 <- conn_mat1[, !colnames(conn_mat1) %in% ROI_list]
  
  rownames(conn_mat1)[length(conn_mat1)] <- ROI
  colnames(conn_mat1)[length(conn_mat1)] <- ROI
  
  return(conn_mat1) 
}



apata <- c()
apata$P <- c('P', 'P-sen', 'P-mot', 'P-sat')
apata$PV <-  c('PVa', 'PVi', 'PVp', 'PVpo')
apata$AM <- c('AM', 'AMd', 'AMv')
apata$MB <- c('MB', 'MBsen', 'MBmot', 'MBsta') 
apata$ENT <- c('ENT', 'ENTl', 'ENTm')
apata$EP <- c('EP', 'EPd', 'EPv')
apata$LS <- c('LS', 'LSc', 'LSr', 'LSv')
apata$IC <- c('IC', 'ICc', 'ICd', 'ICe')
apata$MY <- c('MY', 'MY-sen', 'MY-mot', 'MY-sat')
apata$MO <- c('MO', 'MOp', 'MOs')
apata$PM <- c('PMd', 'PMv')
apata$AP <- c('AP', 'APr')
apata$SC <- c('SCs', 'SCop', 'SCsg', 'SCzo', 'SCm', 'SCdg', 'SCdw', 'SCiw', 'SCig')
apata$VIS <- c('VIS', 'VISal', 'VISam', 'VISl', 'VISp', 'VISpl', 'VISpm', 'VISli', 'VISpor', 'VISa', 'VISrl')
apata$SN <- c('SNr', 'SNc')




for (ROI in unique(ROIs$short)) {
  print(ROI)
  conn_mat_newAllen <- Collapse(ROI, conn_mat_newAllen)
}

for (i in 1:length(apata)) {
  print(names(apata)[i])
  conn_mat_newAllen <- manual(names(apata)[i], apata[[i]], conn_mat_newAllen)
}


write.table(conn_mat_newAllen, file = 'new_allen.csv', quote = F, row.names = T, col.names = T, sep = '\t')
write.table(ROIs, file = 'grouped1.csv', quote = F, row.names = T, col.names = T, sep = '\t')













#nodes df
conn_nodes_newAllen_Full <- data.frame(acronym = row.names(conn_mat_newAllen_Full), stringsAsFactors = F)
conn_nodes_newAllen_Full$name <- NA 
for (i in 1:nrow(conn_nodes_newAllen_Full)) {
  if(conn_nodes_newAllen_Full$acronym[i] %in% allen_dict$acronym){
    conn_nodes_newAllen_Full$name[i] <- allen_dict$name[conn_nodes_newAllen_Full$acronym[i] == allen_dict$acronym]
  } 
}


#edges df
tic()
conn_edges_newAllen_Full <- NULL
for (i in 1:nrow(conn_mat_newAllen)) {
  for (j in 1:nrow(conn_mat_newAllen)) {
    if(!is.na(conn_mat_newAllen[i,j])){
      conn_edges_newAllen_Full <- rbind(conn_edges_newAllen_Full, c(rownames(conn_mat_newAllen)[i], colnames(conn_mat_newAllen)[j], conn_mat_newAllen[i,j]))
    }
  }
}
colnames(conn_edges_newAllen_Full) <- c('Source', 'Target', 'p_values')
conn_edges_newAllen_Full <- data.frame(Source = conn_edges_newAllen_Full[,1], Target = conn_edges_newAllen_Full[,2], p_values = as.numeric(conn_edges_newAllen_Full[,3]))
toc()


#### percolation ####  
tic('calc p-values')
perc_newAllen <- pval_percolation(conn_nodes_newAllen, conn_edges_newAllen_Full)
toc()

#plots
ggplot(data.frame(perc_newAllen$pvals), aes(log10(p.value), SCCSize)) + 
  geom_point() + ggtitle('Biggest SCC') + geom_vline(xintercept = -3.5, linetype = "dashed", color = 'red') +
  geom_vline(xintercept = conn_edges_newAllen['CA1', 'CA3',3], linetype = "dashed", color = 'blue')
  #annotate('text', x = )
 
ggplot(data.frame(perc_newAllen$pvals), aes(log10(p.value), MeanDegree)) + 
  geom_point() + ggtitle('Mean Degree') + geom_vline(xintercept = -3.5, linetype = "dashed", color = 'red') + 
  geom_hline(yintercept = 10, linetype = "dashed", color = 'blue')

ggplot(data.frame(perc_newAllen$pvals), aes(log10(p.value), WeakSize)) + 
  geom_point() + ggtitle('WeakCC') + geom_vline(xintercept = -3.5, linetype = "dashed", color = 'red') #+ 
  #geom_vline(xintercept = -8, linetype = "dashed", color = 'blue') + geom_vline(xintercept = -9, linetype = "dashed", color = 'blue') 
  






ggplot(data.frame(perc_newAllen_Full$pvals), aes(log10(p.value), SCCSize)) + 
  geom_point() + ggtitle('Full Allen Biggest SCC') + geom_vline(xintercept = -3.5, linetype = "dashed", color = 'red')

ggplot(data.frame(perc_newAllen_Full$pvals), aes(log10(p.value), MeanDegree)) + 
  geom_point() + ggtitle('Full Allen Mean Degree') + geom_hline(yintercept = 10, linetype = "dashed", color = 'blue') + 
  geom_vline(xintercept = -3.5, linetype = "dashed", color = 'red') 

ggplot(data.frame(perc_newAllen_Full$pvals), aes(log10(p.value), WeakSize)) + 
  geom_point() + ggtitle('Full Allen WeakCC') + ylim(0, NA) + geom_vline(xintercept = -3.5, linetype = "dashed", color = 'red') 










ggplot(data.frame(perc_connstr$pvals), aes(log10(p.value), SCCSize)) + 
  geom_point() + ggtitle('CONNECTION STRENGTH Biggest SCC')

ggplot(data.frame(perc_connstr$pvals), aes(log10(p.value), MeanDegree)) + 
  geom_point() + ggtitle('CONNECTION STRENGTH Mean Degree') + geom_hline(yintercept = 10, linetype = "dashed", color = 'blue') +
  geom_vline(xintercept = log10(2500), linetype = "dashed", color = 'red')

ggplot(data.frame(perc_connstr$pvals), aes(log10(p.value), WeakSize)) + 
  geom_point() + ggtitle('CONNECTION STRENGTH WeakCC') + ylim(0, NA)









ggplot(data.frame(perc_connstr_inv$pvals), aes(log10(p.value), SCCSize)) + xlab('log10(CONNECTION STRENGTH)') +
  geom_point() + ggtitle('CONNECTION STRENGTH Biggest SCC') +
  geom_vline(xintercept = log10(2500), linetype = "dashed", color = 'red') + #geom_vline(xintercept = log10(4000), linetype = "dashed", color = 'red') + 
  geom_vline(xintercept = log10(conn_mat_connstr['CA3', 'CA1']), linetype = "dashed") +
  geom_vline(xintercept = log10(conn_mat_connstr['ENT', 'DG']), linetype = "dashed") +
  geom_vline(xintercept = log10(conn_mat_connstr['CA3', 'ENT']), linetype = "dashed") #+
#  coord_cartesian(xlim = c(3, NA), ylim = c(NA, NA))

ggplot(data.frame(perc_connstr_inv$pvals), aes(log10(p.value), MeanDegree)) + xlab('log10(CONNECTION STRENGTH)') +
  geom_point() + ggtitle('CONNECTION STRENGTH Mean Degree') + geom_hline(yintercept = 10, linetype = "dashed", color = 'blue') +
  geom_vline(xintercept = log10(2500), linetype = "dashed", color = 'red') + #geom_vline(xintercept = log10(4000), linetype = "dashed", color = 'red') +
  geom_vline(xintercept = log10(conn_mat_connstr['CA3', 'CA1']), linetype = "dashed") +
  annotate('text', x = log10(conn_mat_connstr['CA3', 'CA1']) - .2, y = 100, label = 'CA3 -> CA1', angle = 90) +
  geom_vline(xintercept = log10(conn_mat_connstr['ENT', 'DG']), linetype = "dashed") +
  annotate('text', x = log10(conn_mat_connstr['ENT', 'DG']) + .2, y = 200, label = 'ENT-> DG', angle = 90) +
  geom_vline(xintercept = log10(conn_mat_connstr['CA3', 'ENT']), linetype = "dashed") +
  annotate('text', x = log10(conn_mat_connstr['CA3', 'ENT']) - .2, y = 300, label = 'CA3-> ENT', angle = 90) + coord_cartesian(xlim = c(3, NA), ylim = c(NA, 15))
  

ggplot(data.frame(perc_connstr_inv$pvals), aes(log10(p.value), WeakSize)) + xlab('log10(CONNECTION STRENGTH)') +
  geom_point() + ggtitle('CONNECTION STRENGTH WeakCC') + ylim(0, NA) +
  geom_vline(xintercept = log10(2500), linetype = "dashed", color = 'red') + #geom_vline(xintercept = log10(4000), linetype = "dashed", color = 'red') + 
  geom_vline(xintercept = log10(conn_mat_connstr['CA3', 'CA1']), linetype = "dashed") +
  geom_vline(xintercept = log10(conn_mat_connstr['ENT', 'DG']), linetype = "dashed") +
  geom_vline(xintercept = log10(conn_mat_connstr['CA3', 'ENT']), linetype = "dashed") #+
#  coord_cartesian(xlim = c(3, NA), ylim = c(NA, NA))









ggplot(data.frame(perc_newAllen_Full_inv$pvals), aes(log10(p.value), SCCSize)) + xlab('log10(NORM CONNECTION DENS)') +
  geom_point() + ggtitle('Full Allen Biggest SCC') + geom_vline(xintercept = -3.5, linetype = "dashed", color = 'red') +
  geom_vline(xintercept = log10(conn_mat_newAllen_Full['CA3', 'CA1']), linetype = "dashed") + 
  geom_vline(xintercept = log10(conn_mat_newAllen_Full['ENT', 'DG']), linetype = "dashed") + 
  geom_vline(xintercept = log10(conn_mat_newAllen_Full['CA3', 'ENT']), linetype = "dashed")  
  

ggplot(data.frame(perc_newAllen_Full_inv$pvals), aes(log10(p.value), MeanDegree)) + xlab('log10(NORM CONNECTION DENS)') +
  geom_point() + ggtitle('Full Allen Mean Degree') + geom_hline(yintercept = 10, linetype = "dashed", color = 'blue') + 
  geom_vline(xintercept = -3.5, linetype = "dashed", color = 'red') +
  geom_vline(xintercept = log10(conn_mat_newAllen_Full['CA3', 'CA1']), linetype = "dashed") + 
  annotate('text', x = log10(conn_mat_newAllen_Full['CA3', 'CA1']) - .2, y = 100, label = 'CA3 -> CA1', angle = 90) +
  geom_vline(xintercept = log10(conn_mat_newAllen_Full['ENT', 'DG']), linetype = "dashed") + 
  annotate('text', x = log10(conn_mat_newAllen_Full['ENT', 'DG']) + .2, y = 200, label = 'ENT-> DG', angle = 90) +
  geom_vline(xintercept = log10(conn_mat_newAllen_Full['CA3', 'ENT']), linetype = "dashed") + 
  annotate('text', x = log10(conn_mat_newAllen_Full['CA3', 'ENT']) - .2, y = 300, label = 'CA3-> ENT', angle = 90) +
  coord_cartesian(xlim = c(-4, NA), ylim = c(NA, 100))

ggplot(data.frame(perc_newAllen_Full_inv$pvals), aes(log10(p.value), WeakSize)) + xlab('log10(NORM CONNECTION DENS)') +
  geom_point() + ggtitle('Full Allen WeakCC') + ylim(0, NA) + geom_vline(xintercept = -3.5, linetype = "dashed", color = 'red') +
  geom_vline(xintercept = log10(conn_mat_newAllen_Full['CA3', 'CA1']), linetype = "dashed") + 
  geom_vline(xintercept = log10(conn_mat_newAllen_Full['ENT', 'DG']), linetype = "dashed") + 
  geom_vline(xintercept = log10(conn_mat_newAllen_Full['CA3', 'ENT']), linetype = "dashed")


#### Corr clustesr ####

corr_raq <- read.csv("~/Downloads/engram_retrieval_cormatrix_r_r2_p.csv", stringsAsFactors = F, header = T)
corr_raq <- corr_raq[1:3]
corr_raq <- dcast(corr_raq, row ~ column, value.var = "cor")
rownames(corr_raq) <- t(corr_raq[1])
corr_raq <- corr_raq[-1]

dendrogram_raq <- as.dendrogram(hclust(as.dist((1-corr_raq)/2)))

conn_nodes_newAllen <- get_hier_clus(10, dendrogram_raq, conn_nodes_newAllen)


#### Amended sinClusters #####
sinCluster_raq <- read.csv("ROIs_sinCluster_raq.csv", stringsAsFactors = F, header = T)
sinCluster_raq <- sinCluster_raq[c(1,2,6)]

chng_acr <- sinCluster_raq[!sinCluster_raq$Raq.data %in% rownames(conn_mat_newAllen),]
for (i in 1:nrow(chng_acr)) {
  colnames(conn_mat_newAllen)[colnames(conn_mat_newAllen) == chng_acr$acronym[i]] <- chng_acr$Raq.data[i]
  rownames(conn_mat_newAllen)[rownames(conn_mat_newAllen) == chng_acr$acronym[i]] <- chng_acr$Raq.data[i]
}







#Full new connectome, to ungroup some ROIs
conn_mat_newAllen_Full <- read.csv('~/Downloads/normalized_connection_density.csv', stringsAsFactors = F, row.names = 1)
conn_mat_newAllen_Full <- conn_mat_newAllen_Full[grep('ipsi', colnames(conn_mat_newAllen_Full))]
colnames(conn_mat_newAllen_Full) <- conn_mat_newAllen_Full[1,]
conn_mat_newAllen_Full <- conn_mat_newAllen_Full[-1,]
for (col in colnames(conn_mat_newAllen_Full)){
  conn_mat_newAllen_Full[, col] <- as.numeric(conn_mat_newAllen_Full[, col])
}

tic()
conn_edges_newAllen_Full <- NULL
for (i in 1:nrow(conn_mat_newAllen_Full)) {
  for (j in 1:nrow(conn_mat_newAllen_Full)) {
    if(!is.na(conn_mat_newAllen_Full[i,j])){
      conn_edges_newAllen_Full <- rbind(conn_edges_newAllen_Full, c(rownames(conn_mat_newAllen_Full)[i], colnames(conn_mat_newAllen_Full)[j], conn_mat_newAllen_Full[i,j]))
    }
  }
}
colnames(conn_edges_newAllen_Full) <- c('Source', 'Target', 'p_values')
conn_edges_newAllen_Full <- data.frame(Source = conn_edges_newAllen_Full[,1], Target = conn_edges_newAllen_Full[,2], p_values = as.numeric(conn_edges_newAllen_Full[,3]))
toc() #33min


#Full new connectome CONNECTION_STRENGTH
conn_mat_connstr <- read.csv('~/Downloads/connection_strength.csv', stringsAsFactors = F, row.names = 1)
conn_mat_connstr <- conn_mat_connstr[grep('ipsi', colnames(conn_mat_connstr))]
colnames(conn_mat_connstr) <- conn_mat_connstr[1,]
conn_mat_connstr <- conn_mat_connstr[-1,]
for (col in colnames(conn_mat_connstr)){
  conn_mat_connstr[, col] <- as.numeric(conn_mat_connstr[, col])
}

tic()
conn_edges_connstr <- NULL
for (i in 1:nrow(conn_mat_connstr)) {
  for (j in 1:nrow(conn_mat_connstr)) {
    if(!is.na(conn_mat_connstr[i,j])){
      conn_edges_connstr <- rbind(conn_edges_connstr, c(rownames(conn_mat_connstr)[i], colnames(conn_mat_connstr)[j], conn_mat_connstr[i,j]))
    }
  }
}
colnames(conn_edges_connstr) <- c('Source', 'Target', 'p_values')
conn_edges_connstr <- data.frame(Source = conn_edges_connstr[,1], Target = conn_edges_connstr[,2], p_values = as.numeric(conn_edges_connstr[,3]))
toc() #20min

conn_nodes_connstr <- data.frame(acronym = row.names(conn_mat_connstr), stringsAsFactors = F)
conn_nodes_connstr$name <- NA 
for (i in 1:nrow(conn_nodes_connstr)) {
  if(conn_nodes_connstr$acronym[i] %in% allen_dict$acronym){
    conn_nodes_connstr$name[i] <- allen_dict$name[conn_nodes_connstr$acronym[i] == allen_dict$acronym]
  } 
}


tic('calc p-values')
perc_connstr <- pval_percolation(conn_nodes_connstr, conn_edges_connstr)
toc()


#### inv perc ####
inv_percolation <- function(nodes, edges, is.corr = F) {
  colnames(edges)[3] <- 'p_values'
  if(is.corr){
    edges$p_values <- abs(edges$p_values)
    directed <- F
  } else {
    directed <- T
  }
  conn_nodes_pval <- nodes
  max <- max(edges$p_values)
  pval <- min(edges$p_values)
  print(pval)
  pvals <- NULL
  
  while (1 >= pval & pval <= max) {
    if((round(pval,3)*100)%%5 == 0){
      print(pval)
    }
    conn_nodes_pval <- nodes_info(conn_nodes_pval, edges %>% filter(p_values >= pval), directed = directed)
    conn_nodes_pval <- conn_nodes_pval[!is.na(conn_nodes_pval$SCCId),]
    pvals <- rbind(pvals, c(pval, max(unique(conn_nodes_pval$SCCSize)), max(unique(conn_nodes_pval$WeakSize)), mean(conn_nodes_pval$OutDegree, na.rm = T)))
    if(pval > .001){
      pval <- pval + .001
    } else {
      pval <- pval*10
    }
    if(pval == 0){
      pval <- 10**-10
    }
  }
  
  while (pval >= 1 & pval <= max) {
    print(pval)
    conn_nodes_pval <- nodes_info(conn_nodes_pval, edges %>% filter(p_values >= pval), directed = directed)
    conn_nodes_pval <- conn_nodes_pval[!is.na(conn_nodes_pval$SCCId),]
    pvals <- rbind(pvals, c(pval, max(unique(conn_nodes_pval$SCCSize)), max(unique(conn_nodes_pval$WeakSize)), mean(conn_nodes_pval$OutDegree, na.rm = T)))
    if(pval > 100){
      pval <- pval + 100
    }
    if(100 >= pval & pval > 10){
      pval <- pval + 10
    }
    if(10 >= pval & pval >= 1){
      pval <- pval + 1
    }
  }
  
  list <- NULL
  colnames(pvals) <- c('p-value', 'SCCSize', 'WeakSize', 'MeanDegree')
  list$pvals <- pvals
  list$nodes <- conn_nodes_pval
  return(list)
}



#### 

raq_ROIs <- data.frame(arc = colnames(corr_raq), inAllen = colnames(corr_raq) %in% conn_nodes_connstr$acronym)
raq_ROIs$Collapsed <- NA
for (i in 1:nrow(raq_ROIs)) {
  matches <- grep(paste('^', raq_ROIs$arc[i], sep = ''), conn_nodes_connstr$acronym, value = T) 
  lower <- grep('[a-z]', matches, value = T)
  if(length(lower) & length(matches) > 1){
    raq_ROIs$Collapsed[i] <- paste(lower, collapse = ', ')
  }
}

also_dup <- raq_ROIs[complete.cases(raq_ROIs),]
also_dup$Collapsed[1] <- 'AMd, AMv'
also_dup$Collapsed[6] <- 'COApl, COApm'
also_dup <- also_dup[1:9,]
also_dup <- also_dup[!also_dup$arc %in% c('LG', 'MD', 'PR'),] #LG not in allen

check <- function(ROI, ROI_list, conn_mat1){
  ROI_list <- unlist(strsplit(ROI_list, split = ', '))
  
  rows <- conn_mat1[rownames(conn_mat1) %in% ROI_list, ]
  cols <- conn_mat1[, colnames(conn_mat1) %in% ROI_list]
  
  newrow <- colSums(rows)
#  newcol <- rowSums(cols)
  
#  return(newrow)
  return(sum(newrow - conn_mat1[rownames(conn_mat1) == ROI, ]))
}


for (i in 1:nrow(also_dup)) {
  print(paste(also_dup$arc[i], check(also_dup$arc[i], also_dup$Collapsed[i], conn_mat_connstr)))
#  View(rbind(conn_mat_connstr[rownames(conn_mat_connstr) == also_dup$arc[i], ], check(also_dup$arc[i], also_dup$Collapsed[i], conn_mat_connstr)))
}


#### reducin conn str to only raq ROIs ####
conn_mat_connstr_rdx <- conn_mat_connstr
colnames(conn_mat_connstr_rdx)[colnames(conn_mat_connstr_rdx) == 'SUM'] <- 'SUMl'
rownames(conn_mat_connstr_rdx)[rownames(conn_mat_connstr_rdx) == 'SUM'] <- 'SUMl'

conn_mat_connstr_rdx <- manual('LG', c('LGd', 'LGv'), conn_mat_connstr_rdx)

conn_mat_connstr_rdx <- conn_mat_connstr_rdx[rownames(conn_mat_connstr_rdx) %in% raq_ROIs$arc, ]
conn_mat_connstr_rdx <- conn_mat_connstr_rdx[, colnames(conn_mat_connstr_rdx) %in% raq_ROIs$arc]


tic()
conn_edges_connstr_rdx <- NULL
for (i in 1:nrow(conn_mat_connstr_rdx)) {
  for (j in 1:nrow(conn_mat_connstr_rdx)) {
    if(!is.na(conn_mat_connstr_rdx[i,j])){
      conn_edges_connstr_rdx <- rbind(conn_edges_connstr_rdx, c(rownames(conn_mat_connstr_rdx)[i], colnames(conn_mat_connstr_rdx)[j], conn_mat_connstr_rdx[i,j]))
    }
  }
}
colnames(conn_edges_connstr_rdx) <- c('Source', 'Target', 'p_values')
conn_edges_connstr_rdx <- data.frame(Source = conn_edges_connstr_rdx[,1], Target = conn_edges_connstr_rdx[,2], ConnStr = as.numeric(conn_edges_connstr_rdx[,3]))
toc()

conn_nodes_connstr_rdx <- data.frame(acronym = row.names(conn_mat_connstr_rdx), stringsAsFactors = F)
conn_nodes_connstr_rdx$name <- NA 
for (i in 1:nrow(conn_nodes_connstr_rdx)) {
  if(conn_nodes_connstr_rdx$acronym[i] %in% allen_dict$acronym){
    conn_nodes_connstr_rdx$name[i] <- allen_dict$name[conn_nodes_connstr_rdx$acronym[i] == allen_dict$acronym]
  } 
}

tic()
perc_connstr_rdx_inv <- inv_percolation(conn_nodes_connstr_rdx, conn_edges_connstr_rdx)
toc()


ggplot(data.frame(perc_connstr_rdx_inv$pvals), aes(log10(p.value), MeanDegree)) + xlab('log10(CONNECTION STRENGTH)') +
  geom_point() + ggtitle('CONNECTION STRENGTH REDUCED Mean Degree') + geom_hline(yintercept = 10, linetype = "dashed", color = 'blue') + 
  geom_vline(xintercept = log10(2500), linetype = "dashed", color = 'red') + #geom_vline(xintercept = log10(4000), linetype = "dashed", color = 'red') +
  geom_vline(xintercept = log10(conn_mat_connstr_rdx['CA3', 'CA1']), linetype = "dashed") + 
  annotate('text', x = log10(conn_mat_connstr_rdx['CA3', 'CA1']) - .2, y = 25, label = 'CA3 -> CA1', angle = 90) +
  geom_vline(xintercept = log10(conn_mat_connstr_rdx['ENTl', 'DG']), linetype = "dashed") + 
  annotate('text', x = log10(conn_mat_connstr_rdx['ENTl', 'DG']) + .2, y = 50, label = 'ENT-> DG', angle = 90) +
  geom_vline(xintercept = log10(conn_mat_connstr_rdx['CA3', 'ENTl']), linetype = "dashed") + 
  annotate('text', x = log10(conn_mat_connstr_rdx['CA3', 'ENTl']) - .2, y = 75, label = 'CA3-> ENT', angle = 90) #+ coord_cartesian(xlim = c(2.5, NA), ylim = c(NA, 10))

ggplot(data.frame(perc_connstr_rdx_inv$pvals), aes(log10(p.value), SCCSize)) + xlab('log10(CONNECTION STRENGTH)') +
  geom_point() + ggtitle('CONNECTION STRENGTH REDUCED Biggest SCC') + #geom_hline(yintercept = 10, linetype = "dashed", color = 'blue') + 
  geom_vline(xintercept = log10(2500), linetype = "dashed", color = 'red') + #geom_vline(xintercept = log10(4000), linetype = "dashed", color = 'red') +
  geom_vline(xintercept = log10(conn_mat_connstr_rdx['CA3', 'CA1']), linetype = "dashed") + 
  geom_vline(xintercept = log10(conn_mat_connstr_rdx['ENTl', 'DG']), linetype = "dashed") + 
  geom_vline(xintercept = log10(conn_mat_connstr_rdx['CA3', 'ENTl']), linetype = "dashed") + coord_cartesian(xlim = c(2.5, NA))

ggplot(data.frame(perc_connstr_rdx_inv$pvals), aes(log10(p.value), WeakSize)) + xlab('log10(CONNECTION STRENGTH)') +
  geom_point() + ggtitle('CONNECTION STRENGTH REDUCED WeakCC') + geom_hline(yintercept = 10, linetype = "dashed", color = 'blue') + 
  geom_vline(xintercept = log10(2500), linetype = "dashed", color = 'red') + #geom_vline(xintercept = log10(4000), linetype = "dashed", color = 'red') +
  geom_vline(xintercept = log10(conn_mat_connstr_rdx['CA3', 'CA1']), linetype = "dashed") + 
  geom_vline(xintercept = log10(conn_mat_connstr_rdx['ENTl', 'DG']), linetype = "dashed") + 
  geom_vline(xintercept = log10(conn_mat_connstr_rdx['CA3', 'ENTl']), linetype = "dashed") + coord_cartesian(xlim = c(2.5, NA))


conn_edges_connstr_rdxCutoff <- conn_edges_connstr_rdx %>% filter(ConnStr >= 2500)
conn_nodes_connstr_rdxCutoff <- nodes_info(conn_nodes_connstr_rdx, conn_edges_connstr_rdxCutoff)

write.table(conn_nodes_connstr_rdxCutoff, file = 'conn_str_rdx_nodes.csv', quote = F, row.names = F, col.names = T, sep = '\t')
write.table(conn_edges_connstr_rdxCutoff, file = 'conn_str_rdx_edges.csv', quote = F, row.names = F, col.names = F, sep = '\t')


conn_str_rdx <- main('~/Dropbox (Graduate Center)/Fibers/Mouse_brain/rat/conn_str_rdx_edges.csv', '\t', Weighted = 0)
conn_str_rdx$Nodes <- conn_str_rdx$Nodes %>% group_by(FiberId) %>% mutate(FiberSize = n()) %>% ungroup()
conn_str_rdx$Nodes$FiberId[conn_str_rdx$Nodes$FiberSize == 1] <- NA

conn_nodes_connstr_rdxCutoff$FiberId <- NA
for (i in 1:nrow(conn_nodes_connstr_rdxCutoff)) {
  if(conn_nodes_connstr_rdxCutoff$acronym[i] %in% conn_str_rdx$Nodes$Label){
    conn_nodes_connstr_rdxCutoff$FiberId[i] <- conn_str_rdx$Nodes$FiberId[conn_nodes_connstr_rdxCutoff$acronym[i] == conn_str_rdx$Nodes$Label]
  }
}
conn_nodes_connstr_rdxCutoff <- conn_nodes_connstr_rdxCutoff %>% group_by(FiberId) %>% mutate(FiberSize = n()) %>% ungroup()
conn_nodes_connstr_rdxCutoff$FiberSize[is.na(conn_nodes_connstr_rdxCutoff$FiberId)] <- NA


#### 
conn_mat_newAllen_rdx <- conn_mat_newAllen_Full
colnames(conn_mat_newAllen_rdx)[colnames(conn_mat_newAllen_rdx) == 'SUM'] <- 'SUMl'
rownames(conn_mat_newAllen_rdx)[rownames(conn_mat_newAllen_rdx) == 'SUM'] <- 'SUMl'

conn_mat_newAllen_rdx <- manual('LG', c('LGd', 'LGv'), conn_mat_newAllen_rdx)

conn_mat_newAllen_rdx <- conn_mat_newAllen_rdx[rownames(conn_mat_newAllen_rdx) %in% raq_ROIs$arc, ]
conn_mat_newAllen_rdx <- conn_mat_newAllen_rdx[, colnames(conn_mat_newAllen_rdx) %in% raq_ROIs$arc]


tic()
conn_edges_newAllen_rdx <- NULL
for (i in 1:nrow(conn_mat_newAllen_rdx)) {
  for (j in 1:nrow(conn_mat_newAllen_rdx)) {
    if(!is.na(conn_mat_newAllen_rdx[i,j])){
      conn_edges_newAllen_rdx <- rbind(conn_edges_newAllen_rdx, c(rownames(conn_mat_newAllen_rdx)[i], colnames(conn_mat_newAllen_rdx)[j], conn_mat_newAllen_rdx[i,j]))
    }
  }
}
colnames(conn_edges_newAllen_rdx) <- c('Source', 'Target', 'p_values')
conn_edges_newAllen_rdx <- data.frame(Source = conn_edges_newAllen_rdx[,1], Target = conn_edges_newAllen_rdx[,2], NormDens = as.numeric(conn_edges_newAllen_rdx[,3]))
toc()

write.table(conn_edges_newAllen_rdx, file = 'norm_dens_rdx_edges.csv', quote = F, row.names = F, col.names = F, sep = '\t')
