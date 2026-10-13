CreateNetworkBasis <- function(net,
                               mm_target,
                               N_pop
){
  sampled_nodes<-net$ego
  strata_names <- sort(unique(c(net$ego$strata,net$alter$strata)))
  sampled_nodes$net.degree<-table(factor(net$alter$.egoID,levels = sampled_nodes$.egoID))

  ## Create empty matrix
  mixing_matrix <- matrix(0,
                          nrow = length(strata_names),
                          ncol = length(strata_names),
                          dimnames = list(as.character(strata_names),
                                          as.character(strata_names)))
  ## Populate each element of the matrix with the corresponding target
  for(strata_row in strata_names){
    for(strata_col in strata_names){
      strata_string <- glue("mm[strata={strata_row},strata={strata_col}]")
      mixing_matrix[as.character(strata_row), as.character(strata_col)] <- mm_target[strata_string]
    }
  }
  mixing_matrix[lower.tri(mixing_matrix)]<-t(mixing_matrix[upper.tri(mixing_matrix)])

  stdf<-expand.grid(strata_names,strata_names) #set potential tie df
  stel<-stdf[sample(1:nrow(stdf),size=sum(mixing_matrix),replace=T,prob=mixing_matrix+.001),] #sample strata/strata ties
  sttb<-table(net$ego$strata) #get marginal sums for stratas
  popegos<-cbind(1:N_pop,as.numeric(sample(names(sttb),N_pop, replace=T,prob=sttb))) #sample potential set of ids per strata
  colnames(popegos)<-c("egoID","strata")
  popegos_split<-split(popegos[,1],popegos[,2])
  stel$f_id<-unlist(lapply(stel[,1],function(x) sample(popegos_split[[x]],1))) #sample list of from ids per strata
  stel$t_id<-unlist(lapply(1:nrow(stel),function(x) { #sample "to" ids that do not include "from" id
    sample_pool<-popegos_split[[stel[x,2]]]
    sample_pool<-sample_pool[which(sample_pool!=stel$f_id[x])]
    sample(sample_pool,size=1)}))

  #convert to igraph, then network
  basis_igraph <- graph_from_data_frame(d = stel[, c("f_id", "t_id")],
                                        ## Only keep the node characteristics we need, for speed
                                        vertices = popegos[,c("egoID", "strata")],
                                        directed = FALSE)
  basis_network <- asNetwork(basis_igraph)
  return(basis_network)
}
