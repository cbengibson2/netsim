source("city_sim_singlescript.R")
N_pop<-10000

#use case 1: fully observed egonet
data(faux.mesa.high)
egonet <- as.egor(faux.mesa.high)
net<-citysim(egonet=egonet,
             N_pop=N_pop,
             xnames=c("Grade","Sex"))


#use case 2: sufficient statistics + set of sampled egos
load("basis_inputs.Rdata")
tm<-Sys.time()
N_pop<-100000
nets<-list()
for(mode in names(samstats)){
  nets[[mode]]<-citysim(egos=sampled_nodes, #initial list of egos (with 'xnames' as columns)
                        samstats=samstats[[mode]], #population targets
                        N_pop=N_pop)
}
tm-Sys.time()
