rm(list=ls())
library(readxl)
path<-"prem_mixing_matrices.xlsx"
stats<-lapply(excel_sheets(path),function(x) as.data.frame(read_excel(path,sheet=x)))
names(stats)<-gsub("m_","",excel_sheets(path))

#adjustment 1. mixing matrix
load("pre_basis_inputs.Rdata")
load("eldemog2.Rdata")

eldemog2$mode<-as.character(eldemog2$Purpose)
eldemog2$mode[eldemog2$mode%in%c("Visit","Social/Recreation","College","Pick up or drop off a passenger","Shop")]<-"Other"
modes<-split(eldemog2,eldemog2$mode)
names(modes)<-tolower(names(modes))
edge.counts<-lapply(modes,nrow)
names(edge.counts)<-tolower(names(edge.counts))

adj<-(nrow(sampled_nodes)/1601330)
targets<-list()
target<-list()
for(mode in unique(names(edge.counts))){
  target$mm<-((t(stats$pop.prop[mode])*stats[[mode]][,-1])*sum(edge.counts[[mode]]))*adj
  #testnet<-(intergraph::asNetwork(igraph::graph_from_data_frame(d = modes[[mode]][, c("Person.Id.2", "Person.Id.1")],directed = FALSE)))
  #summary(testnet~triangles)
  target$triangle<-(length(igraph::triangles(igraph::graph_from_data_frame(d = modes[[mode]][, c("Person.Id.2", "Person.Id.1")],directed = FALSE)))/3)*adj
  targets[[mode]]<-target
}

#mm[upper.tri(mm)]<-mm[lower.tri(mm)]<-(mm[upper.tri(mm)]+mm[lower.tri(mm)])
sampled_nodes$agecat<-1
sampled_nodes$agecat[sampled_nodes$Age>=20]<-2
sampled_nodes$agecat[sampled_nodes$Age>=30]<-3
sampled_nodes$agecat[sampled_nodes$Age>=40]<-4
sampled_nodes$agecat[sampled_nodes$Age>=50]<-5
sampled_nodes$agecat[sampled_nodes$Age>=60]<-6
sampled_nodes$agecat[sampled_nodes$Age>=70]<-7
sampled_nodes$strata<-sampled_nodes$agecat

#need triangles per mode
#make "basis_inputs" include mode x c(mm,triangles) only
samstats<-targets
save(sampled_nodes,samstats,file="basis_inputs.Rdata")

