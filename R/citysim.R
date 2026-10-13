citysim <- function(egonet=NULL, #initial, egocentric network (egor class)
                    N_pop=NULL, #population N
                    xnames=NULL, #node variables to match mixing statistics
                    max_esp=50) {

  #build strata from "xnames"
   yz<-rep(1,nrow(egonet$ego))
   form1 <- as.formula(paste("yz ~ ", paste(xnames, collapse= "+")))
   marginals<-aggregate(form1,egonet$ego,sum)
   marginals$strata<-1:nrow(marginals)
   marginals<-marginals[,-which(names(marginals)%in%'yz')]
   egonet$ego<-merge(marginals,egonet$ego,all.y=T)
   egonet$alter<-merge(marginals,egonet$alter,all.y=T)

   #build target list
   target <- list()
   target[["mm"]] <- summary(egonet ~ mm("strata", levels2 = TRUE),scaleto=N_pop)
   target[["esp"]] <- summary(egonet ~ esp(1:max_esp),scaleto=N_pop)
   target[["triangle"]] <- sum(as.numeric(gsub("esp","",names(target[["esp"]])))*target$esp[1:length(target$esp)])/3

   #create basis network with attribute mixing
   basis_network <- CreateNetworkBasis(egonet,mm_target = target$mm,N_pop=N_pop)

   #rewire triangles, holding constant degree number of edges
   net<-FitSan(basis_network, target)

  return(net)
}

