library(phytools)
library(svglite)
#source("http://www.phytools.org/cophylo/v0.1/cophylo.R")

#Load and prepare trees
tr_A <- ape::read.tree("pipolbs_firmictues_feb2022_out_mafft-linsi_trimAl.aln.contree")
tr_B <- ape::read.tree("Srec_firmicutes_feb2022_out_mafft-linsi_TrimAl.fa.contree")
tr_A_root <- reroot(tr_A, 1, 0.5)
tr_B_root <- reroot(tr_B, 1, 0.5)
plot(tr_A_root)
plot(tr_B_root)

#Read and prepare links between tips
AB_links_df <- read.table("pol_Srec_assoc.txt", header=FALSE)
AB_links_df_noColname <- AB_links_df; colnames(AB_links_df_noColname) <- NULL
AB_links_matrix <- as.matrix(AB_links_df_noColname)

#Check links (in 16SrRNA-piPolB setdiff should be 0)
list_A_labels <- tr_A_root$tip.label
list_links_A <- AB_links_df$V1
setdiff(list_A_labels, list_links_A)
setdiff(list_links_A, list_A_labels)

#Create cophylogeny
obj<-cophylo(tr_A_root, tr_B_root, rotate=TRUE, assoc = AB_links_matrix)
svg("piPolBs_vs_MOBP.svg")
plot(obj, fsize=0.1, lwd = 0.25, link.col=make.transparent("blue",0.6),
     link.lwd = 0.5, link.lty = 1, link.type="curved", scale.bar=c(0.5, 0.5))
dev.off()



cophyloplot(tr_pol_root, tr_rel_root, show.tip.label = FALSE,
            fsize=0.01, assoc = pol_rel_assoc) #use.edge.length = TRUE


###Test añadir ramas
tr_rel <- ape::read.tree("profile_MOBM_Firmicutes_relaxases_clean_out_mafft-linsi_trimal.faa.contree")
plot(tr_rel)
tr_rel_root <- reroot(tr_rel, 1, 0.5)
plot(tr_rel_root)
tr_pol_root_extra <- bind.tip(tr_rel_root, "attatch", edge.length = 1, where = "root", position = 0)
plot(tr_pol_root_extra)
tr_pol_root_extra_2 <- reroot(tr_pol_root_extra, 18, 0.5) # 18 length
plot(tr_pol_root_extra_2)


### Test unión
tr1<-rtree(10)
tr2<-rtree(10)
tr1$tip.label[1]<-"NA"
tr2$root.edge<-0
tr3<-paste.tree(tr1,tr2)

tr_pol_root_extra_2$tip.label[18] <- "NA"
plot(tr_pol_root_extra_2)

tr_pol_root_extra_2_2 <- tr_pol_root_extra_2
tr_pol_root_extra_2_2$root.edge<-0

tr3<-paste.tree(tr_pol_root_extra_2,tr_pol_root_extra_2_2)
plot(tr3)


### Full process
tr_MOBP <- ape::read.tree("T4SS_MOBP1_Firmicutes_relaxases_clean_out_mafft-linsi_trimal.faa.contree")
tr_MOBP_root <- reroot(tr_MOBP, 1, 0.5)
tr_MOBP_root_att <- bind.tip(tr_MOBP_root, "ATTATCH", edge.length = 0.1, where = "root", position = 0)
tr_MOBP_reroot_att <- reroot(tr_MOBP_root_att, length(tr_MOBP_root_att$tip.label), 0.05)
plot(tr_MOBP_reroot_att, cex = 0.2)

tr_MOBT <- ape::read.tree("profile_MOBT_Firmicutes_relaxases_clean_out_mafft-linsi_trimal.faa.contree")
tr_MOBT_root <- reroot(tr_MOBT, 1, 0.5)
tr_MOBT_root_att <- bind.tip(tr_MOBT_root, "ATTATCH", edge.length = 0.1, where = "root", position = 0)
tr_MOBT_reroot_att <- reroot(tr_MOBT_root_att, length(tr_MOBT_root_att$tip.label), 0.05) # 18 length
plot(tr_MOBT_reroot_att, cex = 0.4)

tr_MOBP_reroot_att$tip.label[length(tr_MOBP_root_att$tip.label)]<-"NA"
plot(tr_MOBP_reroot_att, cex = 0.4)
tr_MOBT_reroot_att$root.edge<-0
tr_MOBP_MOBT <-paste.tree(tr_MOBP_reroot_att,tr_MOBT_reroot_att)
plot(tr_MOBP_MOBT, cex = 0.2)

tr_MOBM <- ape::read.tree("profile_MOBM_Firmicutes_relaxases_clean_out_mafft-linsi_trimal.faa.contree")
tr_MOBM_root <- reroot(tr_MOBM, 1, 0.5)


tr_MOBP_MOBT$tip.label[length(tr_MOBP_MOBT$tip.label)]<-"NA"
plot(tr_MOBP_MOBT, cex = 0.2)
tr_MOBM_root$root.edge<-0
tr_MOBP_MOBT_MOBM <-paste.tree(tr_MOBP_MOBT,tr_MOBM_root)
plot(tr_MOBP_MOBT_MOBM, cex = 0.2)



###COphylogeny
pol_rel_assoc_df <- read.table("pol_rel_assoc_MOB_ALL.txt")
colnames(pol_rel_assoc_df) <- NULL
pol_rel_assoc <- as.matrix(pol_rel_assoc_df)

obj<-cophylo(tr_pol_root, tr_MOBP_MOBT_MOBM, rotate=TRUE, assoc = pol_rel_assoc)
svg("piPolBs_vs_all_MOB.svg")
plot(obj, fsize=0.1, lwd = 0.25, link.col=make.transparent("blue",0.6),
     link.lwd = 0.5, link.lty = 1, link.type="curved")
dev.off()


