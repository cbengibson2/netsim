library(ergm.ego)
library(intergraph)
#library(igraph)

citysim <- function(egonet=NULL, #initial, egocentric network (egor class)
                    egos=NULL, #initial list of egos (with 'xnames' as columns)
                    samstats=NULL,
                    N_pop=NULL, #population N
                    xnames=NULL, #node variables to match mixing statistics
                    upper_esps=1000, #absolute upper bound on the expected number of edgewise shared partners
                    iteration_max= 50,
                    target_tolerance = 0.02,
                    change_tolerance = 0.01,
                    san_iter = 5, #how many times max to run SAN
                    nsteps = N * 2^8, #steps in SAN
                    verbose = TRUE) {

    if((is.null(egonet) & is.null(samstats)) | (!is.null(egonet) & !is.null(samstats))) stop("Please supply either 'egonet' in egor format, or 'egos' AND 'samstats'")

    if(is.null(egos)) egos<-egonet$ego



    #build strata from "xnames", if supplied
    if(!("strata"%in%names(egos))){
      yz<-rep(1,nrow(egos))
      form1 <- as.formula(paste("yz ~ ", paste(xnames, collapse= "+")))
      strata_assignment<-aggregate(form1,egos,sum)
      strata_assignment$strata<-1:nrow(strata_assignment)
      strata_assignment<-strata_assignment[,-which(names(strata_assignment)%in%'yz')]
      if(!is.null(egonet)){
        egonet$ego<-merge(strata_assignment,egonet$ego,all.y=T)
        egonet$alter<-merge(strata_assignment,egonet$alter,all.y=T)
      }
      egos<-merge(strata_assignment,egos,all.y=T)
    }

    #build a simulated population of nodes
    sttb<-table(egos$strata)
    strata_names <- names(sttb)
    popegos<-cbind(1:N_pop,as.numeric(sample(strata_names,N_pop, replace=T,prob=sttb))) #sample potential set of ids per strata
    colnames(popegos)<-c("egoID","strata")

    ##########################################
    ###Produce sufficient statistic targets

    population_adjustment<-N_pop/nrow(egos) #adjustment factor

    #adjust stats for population size if supplied, format
    if(!is.null(samstats)){
      target<-lapply(samstats,function(x) x*population_adjustment)
      names(target)<-c("mm","triangle")
      mixing_matrix<-as.matrix(target[["mm"]])
    }

    #get target stats from egonet, if supplied
    if(!is.null(egonet)){
      mixing_matrix <- matrix(0,
                              nrow = length(strata_names),
                              ncol = length(strata_names),
                              dimnames = list(as.character(strata_names),
                                              as.character(strata_names)))
      #get counts of mixing by strata
      target <- list()
      target[["mm"]] <- summary(egonet ~ mm("strata", levels2 = TRUE),scaleto=N_pop)
      #find esps for triangle counting
      esps<-summary(egonet ~ esp(1:upper_esps),scaleto=N_pop)
      max_esp<-max(which(esps[1:upper_esps]>0))+5
      esps <- summary(egonet ~ esp(1:max_esp),scaleto=N_pop)
      target[["triangle"]] <- sum(as.numeric(gsub("esp","",names(esps)))*esps[1:length(target$esps)])/3
      target[["triangle"]]<-target[["triangle"]]*population_adjustment #adjust mm for population size

      ## Populate each element of the matrix with the corresponding target
      mm_target<-target[["mm"]]*population_adjustment #adjust mm for population size
      names(mm_target)<-gsub("\\]","",gsub("strata=","",gsub("mm\\[","",names(mm_target))))
      for(strata_row in strata_names){
        for(strata_col in strata_names){
          strata_string <- paste0(strata_row,",",strata_col)
          mixing_matrix[as.character(strata_row), as.character(strata_col)] <- mm_target[strata_string]
        }
      }
      mixing_matrix[lower.tri(mixing_matrix)]<-t(mixing_matrix[upper.tri(mixing_matrix)])
    }

  ######
  ###Produce basis network

  stdf<-expand.grid(strata_names,strata_names) #set potential tie df
  stel<-stdf[sample(1:nrow(stdf),size=sum(mixing_matrix),replace=T,prob=mixing_matrix+.001),] #sample strata/strata ties
  popegos_split<-split(popegos[,1],popegos[,2])
  stel$f_id<-unlist(lapply(stel[,1],function(x) sample(popegos_split[[x]],1))) #sample list of from ids per strata
  stel$t_id<-unlist(lapply(1:nrow(stel),function(x) { #sample "to" ids that do not include "from" id
    sample_pool<-popegos_split[[stel[x,2]]]
    sample_pool<-sample_pool[which(sample_pool!=stel$f_id[x])]
    sample(sample_pool,size=1)
    }
  ))

  #convert to igraph, then network
  basis_igraph <- igraph::graph_from_data_frame(d = stel[, c("f_id", "t_id")],
                                        ## Only keep the node characteristics we need, for speed
                                        vertices = popegos[,c("egoID", "strata")],
                                        directed = FALSE)
  basis_network <- intergraph::asNetwork(basis_igraph)


  ################################ Fit Networks
  #rewire triangles, holding constant degree and number of edges
  ## Fit triangle target

  sim_tri_network <- basis_network
  triangle_target_ratio <- 0
  iteration_count <- 0
  previous_triangles <- summary(basis_network ~ triangles)
  triangle_change_ratio <- 1
  N<-length(basis_network$val)

  ## Create formula for fitting triangles
  form <- sim_tri_network ~ triangles
  ## Initialize target statistic vector in ergm format
  target_stats <- summary(form)
  ## Replace target statistic with observed statistic
  target_stats[1] <- target$triangle

  while((triangle_target_ratio < (1-target_tolerance))|(triangle_target_ratio > (1 + target_tolerance))){

    ## Run the SAN
    sim_tri_network <- san(form,
                   target.stats = target_stats,
                   only.last = TRUE,
                   constraints = ~ edges + degreedist,
                   control = list(
                     ## Use low temperature because greedy searches seem to work fine for triangles
                     SAN.tau = 0.1,
                     SAN.maxit = san_iter,
                     SAN.invcov.diag = T,
                     SAN.nsteps = nsteps,
                     ## Evenly divide the steps
                     SAN.nsteps.alloc = function(nsim) 1^seq_len(nsim)
                   ),
                   verbose = verbose
    )

    observed_triangles <- summary(sim_tri_network ~ triangles)
    triangle_target_ratio <- observed_triangles/target$triangle

    ## Stop if there has been little change in the number of triangles
    triangle_change_ratio <- observed_triangles/previous_triangles
    triangle_change_abs <- abs(triangle_change_ratio-1)
    if(triangle_change_abs < change_tolerance){
      cat("Stopping due to small change in triangles\n")
      break
    }

    previous_triangles <- observed_triangles

    ## Stop if going beyond max iterations
    iteration_count <- iteration_count + 1

    cat(paste0("iteration count ", iteration_count, ", change since last = ",round(triangle_change_abs*100, 2),"%"))
    if(iteration_count >= iteration_max){
      cat("Stopping at max iterations\n")
      break
    }
  }

  #add back xnames
  xnames_replace<-xnames[!(xnames%in%network::list.vertex.attributes(sim_tri_network))]
  for(xname in xnames_replace) {
    vals<-strata_assignment[match(network::get.vertex.attribute(sim_tri_network,"strata"),strata_assignment$strata),xname]
    network::set.vertex.attribute(sim_tri_network,xname,vals)
  }

  return(sim_tri_network)
}

