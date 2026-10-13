library(ergm.ego)
library(intergraph)
library(data.table)

simulate_SIR <- function(edges, beta = 0.3, gamma = 0.1, nsteps = 50,
                         initial_infected = 1, verbose = TRUE) {
  # Ensure data.table
  edges <- as.data.table(edges)
  setnames(edges, c("Var1", "Var2"))

  # Build undirected adjacency list
  adj <- rbind(
    edges[, .(node = Var1, neighbor = Var2)],
    edges[, .(node = Var2, neighbor = Var1)]
  )[, .(neighbors = list(unique(neighbor))), by = node]

  nodes <- adj$node

  # Initialize states
  state <- data.table(node = nodes, status = "S")
  state[sample(.N, initial_infected), status := "I"]

  # Record time series
  ts <- data.table(step = 0,
                   S = sum(state$status == "S"),
                   I = sum(state$status == "I"),
                   R = sum(state$status == "R"))

  # Simulation loop
  for (t in 1:nsteps) {
    infected <- state[status == "I", node]
    susceptible <- state[status == "S", node]

    new_infected <- character(0)
    for (i in infected) {
      neigh <- adj[node == i, neighbors][[1]]
      sus_neigh <- neigh[neigh %in% susceptible]
      newly <- sus_neigh[runif(length(sus_neigh)) < beta]
      new_infected <- c(new_infected, newly)
    }

    recovered <- infected[runif(length(infected)) < gamma]

    state[node %in% new_infected, status := "I"]
    state[node %in% recovered, status := "R"]

    ts <- rbind(ts, data.table(
      step = t,
      S = sum(state$status == "S"),
      I = sum(state$status == "I"),
      R = sum(state$status == "R")
    ))

    if (verbose)
      cat(sprintf("Step %d: S=%d I=%d R=%d\n",
                  t, ts$S[t+1], ts$I[t+1], ts$R[t+1]))

    if (sum(state$status == "I") == 0) break
  }

  # --- Summary measures ---
  peak_I <- max(ts$I)
  time_to_peak <- ts$I %>% which.max() - 1
  final_size <- max(ts$R)
  duration <- max(ts$step)
  total_nodes <- nrow(state)

  summary <- data.table(
    peak_I = peak_I,
    peak_I_frac = peak_I / total_nodes,
    time_to_peak = time_to_peak,
    final_size = final_size,
    final_size_frac = final_size / total_nodes,
    duration = duration
  )

  return(list(timeseries = ts, summary = summary))
}


citysim_findtransitivity <- function(egonet=NULL, #initial, egocentric network (egor class)
                    egos=NULL, #initial list of egos (with 'xnames' as columns)
                    samstats=NULL,
                    N_pop=NULL, #population N
                    xnames=NULL, #node variables to match mixing statistics
                    upper_esps=1000, #absolute upper bound on the expected number of edgewise shared partners
                    iteration_max= 50,
                    target_tolerance = 0.02,
                    change_tolerance = 0.01,
                    san_iter = 5, #how many times max to run SAN
                    nsteps = N_pop * 2^8, #steps in SAN
                    verbose = TRUE,
                    transitivity_value) {

  if((is.null(egonet) & is.null(samstats)) | (!is.null(egonet) & !is.null(samstats))) stop("Please supply either 'egonet' in egor format, or 'egos' AND 'samstats'")
  if(is.null(egos)) egos<-egonet$ego
  sim_tri_networkl<-stats_iter<-list()
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

  ## Replace simulated statistic with target statistic
  deg<-table(c(stel[, c("f_id")], stel[, c("t_id")]))
  triples <- sum(deg * (deg - 1) / 2)
  target_triangles<-((transitivity_value*triples)/3)

  if(target_triangles > max(0,target_stats[1])){
    target_stats[1]<-target_triangles
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
        cat("Stopping due to small change in triangles \n")
        break
      }

      previous_triangles <- observed_triangles

      ## Stop if going beyond max iterations
      iteration_count <- iteration_count + 1
      #stats_iter[[iteration_count]]<- summary(sim_tri_network~triangles+degree(0:100)+esp(1:200)+density)
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
    stats_iter<-c(mode=mode,N_pop=N_pop,target_trans=transitivity_value,summary(sim_tri_network~triangles+degree(0:100)+esp(1:200)+density))

    #run SIR
    edges <- data.table(
      Var1 = intergraph::asDF(sim_tri_network)$edges[,1],
      Var2 = intergraph::asDF(sim_tri_network)$edges[,2]
    )
    res <- simulate_SIR(edges, beta = 0.3, gamma = 0.1, nsteps = 30)
    stats_iter<-c(stats_iter,unlist(res$summary))
    return(list(stats_iter,sim_tri_network))
  }
}

#use case 2: sufficient statistics + set of sampled egos
dir<-getwd()
load(paste0(dir,"/basis_inputs.Rdata"))
netsl<-list()
ctr<-0
for(N_pop in c(1000,2000)){ #,10000,20000,50000,100000
  for(transitivity_value in c(.00001,.0001,.001,.01,.1,.25)){
   nets<-list()
   samstats_iter<-samstats
   for(mode in names(samstats)){
    #solve for # of triangles
     nets[[mode]]<-citysim_findtransitivity(egos=sampled_nodes, #initial list of egos (with 'xnames' as columns)
                          samstats=samstats_iter[[mode]], #statistics
                          N_pop=N_pop,
                          transitivity_value=transitivity_value)
   }
   ctr<-ctr+1
   netsl[[ctr]]<-nets
  }
}

lapply(netsl,function(x) x[[1]])
