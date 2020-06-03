source("cactus_clone_assignment.R")
if (!require(cardelino)) install.packages('cardelino')
library(cardelino)
args = commandArgs(trailingOnly=TRUE)
#args = c("~/Projects/GECCO/intersection/wes_inter.rds","~/Projects/GECCO/preprocess/sc/sc_data.rds","~/Projects/GECCO/cardelino/","~/Projects/GECCO/Canopy/", "S12118", "~/Projects/GECCO/Clustering/cell_cluster.rds","~/Projects/GECCO/code/gecco_clone_id.R","1")
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
CANOPY <- args[4]
CACTUSDIR <- args[3]
wes_inter <- readRDS(args[1])
sc <- readRDS(args[2])
nTree = 30
cell_cluster <- readRDS(args[6])
#whichCluster could be a arg for more than one  clusters 
whichCluster <- 1
relax_rate_prior <- c(0.5,9.5)#c(as.numeric(args[7], digits=15),as.numeric(args[8], digits=15))
dir.create(CACTUSDIR)
subject <- args[5]

Z <- readRDS(paste0(CANOPY,subject,"/Z.rds"))
RESULTS <- paste0(CACTUSDIR,subject,"/")
dir.create(RESULTS)
sample <- wes_inter[[subject]]
sample_sc <- sc[[subject]][sc[[subject]]$MUTATION_ID %in% sample$MUT_ID_minus,]


if (!(file.exists(paste0(RESULTS,"A_clone.rds")) & file.exists(paste0(RESULTS,"D_clone.rds"))) ) {
  A_clone <- data.frame(unique(sample_sc$cell))
  D_clone <- data.frame(unique(sample_sc$cell))
  rownames( A_clone) <- unique(sample_sc$cell)
  rownames( D_clone) <- unique(sample_sc$cell)
  A_clone <- cbind(A_clone,matrix(0, nrow(A_clone), length(sample$MUT_ID_minus)))
  D_clone <- cbind(D_clone,matrix(0, nrow(D_clone), length(sample$MUT_ID_minus)))
  colnames(A_clone) <- c("Cell",sample$MUT_ID_minus)
  colnames(D_clone) <- c("Cell",sample$MUT_ID_minus)
  
  for(read in 1:nrow(sample_sc)){
    d <- D_clone[D_clone$Cell==sample_sc[read,]$cell,][[sample_sc[read,]$MUTATION_ID]]
    D_clone[D_clone$Cell==sample_sc[read,]$cell,][[sample_sc[read,]$MUTATION_ID]] <- d+1
    if(sample_sc[read,]$refAllele!=sample_sc[read,]$base){
      a <- A_clone[A_clone$Cell==sample_sc[read,]$cell,][[sample_sc[read,]$MUTATION_ID]]
      A_clone[A_clone$Cell==sample_sc[read,]$cell,][[sample_sc[read,]$MUTATION_ID]] <- a+1
    }
  }
  saveRDS(A_clone,paste0(RESULTS,"A_clone.rds"))
  saveRDS(D_clone,paste0(RESULTS,"D_clone.rds"))
}else{
  A_clone <- readRDS(paste0(RESULTS,"A_clone.rds"))
  D_clone <- readRDS(paste0(RESULTS,"D_clone.rds"))
  rownames( A_clone) <- unique(sample_sc$cell)
  rownames( D_clone) <- unique(sample_sc$cell)
}
for (j in 1:nTree){
  print(paste0("tree ",j))
  cell_cluster_selected <- cell_cluster[[subject]][,c(1,whichCluster)]
  colnames(cell_cluster_selected) <- c("cell","cluster")
  assignments <- cactus_clone_assignment(t(A_clone[,2:ncol(A_clone)]), t(D_clone[,2:ncol(D_clone)]), Config = Z[[j]],cell_cluster_selected,n_chain = 10,n_proc = 30, relax_rate_prior = relax_rate_prior)

  dffg <- assign_cells_to_clones(assignments$prob,0.2)
  head(dffg)
  saveRDS(assignments,paste0(RESULTS,"tree_",j,"_assignments.rds"))
  saveRDS(dffg,paste0(RESULTS,"tree_",j,"_cellToClone.rds"))
  REPORTS   <-  paste0(RESULTS,"reports/")
  dir.create(REPORTS)
  pdf(file=paste0(REPORTS,"tree_",j,"_cellassignment.pdf"))
  print(prob_heatmap(assignments$prob))
  dev.off()       
}






