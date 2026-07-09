library(ape)
library(phytools)
library(tidyverse)
library(RColorBrewer)
library(phylotools)
library(diversitree)
# library(corHMM)
library(ggbiplot)
library(factoextra) # extracts eingen values from pca objects
library(janitor) # row_to_names function
library(ggrepel)
library(viridis)


# setwd("/mnt/data/Dropbox/Pracovne/Venca/Project2022/HGT_tree_mapping")

setwd("/home/kk/Dropbox/Pracovne/Venca/Project2022/HGT_tree_mapping")

# setwd("/home/karol/Dropbox/Pracovne/Venca/Project2022/HGT_tree_mapping")

# read the HGT_tree_mapping

nexus_tree <- readNexus(file="hordeum_4dtv_time_ordered.nex", format="raxml")

tr2<-rotateNodes(nexus_tree,nodes="all")

plotTree(tr2, fsize = 0.5, direction="rightwards", lwd = 1)

# remove polyploids, multiple accessions or phased chromosomes
# removed Hjub1,2 - polyploid; Hmur1,2 - polyploid; Hsec1,2 - polyploid; HmarB - individual with no hit; Hbul2 - phased one haplotype omited

tips_torm <- c("Hjub1", "Hjub2", "Hmur1", "Hmur2","Hsec1", "Hsec2", "HmarB", "Hbul2", "Sorgh", "Oryza", "Brach", "Triti")

tr3 <- drop.tip(tr2, tips_torm)

plotTree(tr3, fsize = 0.5, direction="rightwards", lwd = 1)



# read table with the HGT data

data <- read_delim("241018_confirmed_HGT_data_edKK.csv", delim = "\t", col_names = T)

# remove data of polyploids and phased data

data_torm <- c("Hordeum_jubatum", "Hordeum_murinum_poly", "Hordeum_secalinum", "Hordeum_bulbosum2")

data <- data %>% filter(!hordeum_species %in% data_torm) %>% 
        mutate(hordeum_chromosome = if_else(grepl("contig",hordeum_chromosome),"other", hordeum_chromosome)) %>%
        mutate(hordeum_chromosome = replace(hordeum_chromosome, hordeum_chromosome=="chr1H_1","chr1H"))
        

data_hits <- data %>% unite(hit, orthogroup,hordeum_loc_start, sep="_", remove=FALSE) %>% 
             select(hordeum_species, hit, hordeum_chromosome) %>%
             mutate(row=row_number())

data_hits_w <- spread(data_hits, key=hordeum_chromosome, value=hit, fill=NA) %>% select(-row) %>%
               mutate(across(.cols = starts_with("chr"), ~replace(.,!is.na(.),1))) %>%
               mutate(across(.cols = starts_with("chr"), ~replace(.,is.na(.),0))) %>%
               mutate(other=replace(other, !is.na(other), 1)) %>%
               mutate(other=replace(other, is.na(other), 0)) %>%
               mutate_at(c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H", "other"),as.numeric) %>%
               group_by(hordeum_species) %>% summarize_all(sum) %>%
               add_row(hordeum_species = "Hordeum_vulgare", chr1H = 0, chr2H = 0, chr3H = 0, chr4H = 0, chr5H = 0, chr6H = 0,
                        chr7H = 0, other = 0) %>%
               add_row(hordeum_species = "Hordeum_gussoneanum", chr1H = 0, chr2H = 0, chr3H = 0, chr4H = 0, chr5H = 0, chr6H = 0,
                        chr7H = 0, other = 0) %>%
               select(-other)
               
# data genes

data_genes <- data %>% select(hordeum_species, orthogroup) %>% group_by(hordeum_species, orthogroup) %>% summarize(hits=n())

data_genes_w <- data_genes %>% spread(key=orthogroup, value= hits, fill=NA) %>%
                mutate(across(.cols = contains("HOG"), ~replace(.,is.na(.),0)))

# save a table showing the number of hits for each othogroup in each hordeum_species

write.table(data_genes_w, "genes_in_species.csv", sep ="\t", col.names = T, row.names = F)

# replace any number higher than zero by 1

data_genes_w <- data_genes_w %>% mutate(across(.cols = contains("HOG"), ~replace(.,. > 0,1)))

# add rows for Hordeum_vulgare and Hordeum_gussoneanum 

data_genes_w <- map_df(data_genes_w, ~ c(0,0, .x))

data_genes_w[1,1]= "Hordeum_vulgare"
data_genes_w[2,1]= "Hordeum_gussoneanum"


# rename tips to fit the HGT data

# create table to rename tips

new_names <- str_sort(unique(data_genes_w$hordeum_species))
tips <- str_sort(tr3$tip.label)

torename <- data.frame(tips, new_names)
colnames(torename) <- c("tip_label", "new_name")

# rename labels
tr_final <- sub.taxa.label(tr3,torename)

# plot tree with node labels

pdf("tree_for_fig.pdf", paper  = "a4r", width = 20)
plot(tr_final,no.margin=TRUE,edge.width=2,cex=0.7, fsize=12)
dev.off()
nodelabels(text=1:tr_final$Nnode,node=1:tr_final$Nnode+Ntip(tr_final))

write.tree(tr_final, "tree_for_fig.nwk")

## map the occurence of the genes individually..

# prepare data

aaa <- data.frame(select(data_genes_w, N0.HOG0007888))
rownames(aaa) <- c(data_genes_w$hordeum_species)
aaa <- mutate_if(aaa, is.numeric, as.factor)

# fit_joint_test <- corHMM(tr_final, aaa, node.states="joint", rate.cat=1, rate.mat=matrix(c(NA,1,1,NA),2,2,)) # from the book

bbb<-setNames(aaa$N0.HOG0007888,rownames(aaa))
levels(bbb)


# fitting the model

fitER<-fitMk(tr_final,bbb,model="ER")
fitER

# reconstruction of the ancestral states
# joint ancestral states

fit.joint<-ancr(fitER,type="joint")
fit.joint

# plot the tree
# cols<-setNames(viridisLite::viridis(n=2),levels(feed.mode)) # according to Revell

cols <- setNames(brewer.pal(n = 2, name = "Dark2"), levels(bbb))

plotTree.datamatrix(tr_final,as.data.frame(bbb),
  colors=list(cols),header=FALSE, fsize=1)
legend("topright",legend=levels(bbb),pch=22,title="N0.HOG0007888",
  pt.cex=2.5,pt.bg=cols,bty="n",cex=1.5)
nodelabels(pie=to.matrix(fit.joint$ace,levels(bbb)),
  piecol=cols,cex=0.4)


# reconstruction of the ancestral states
# joint ancestral states

fit.marginal<-ancr(fitER,type="marginal")
fit.marginal

# save the node states

ccc <- as.data.frame(fit.marginal$ace)
colnames(ccc) <- c("zero", "one")
rownames(ccc) <- c(rownames(fit.marginal$ace))

ddd <- data.frame(ccc$zero)
colnames(ddd) <- c("N0.HOG0007888")
rownames(ddd) <- rownames(ccc)



# plot the tree

# cols<-setNames(viridisLite::viridis(n=2),levels(feed.mode))

cols <- setNames(brewer.pal(n = 2, name = "Dark2"), levels(bbb))

plotTree.datamatrix(tr_final,as.data.frame(bbb),
  colors=list(cols),header=FALSE, fsize=1)
legend("topright",legend=levels(bbb),pch=22,title="N0.HOG0007888",
  pt.cex=2.5,pt.bg=cols,bty="n",cex=1.5)
nodelabels(pie=fit.marginal$ace,piecol=cols,cex=0.4)


# run it in a loop
JOINT <- "241119_joint"
MARGINAL <- "241119_marginal"
dir.create(file.path(JOINT))
dir.create(file.path(MARGINAL))

df_joint <- data.frame(tr_final$node.label)
df_marginal <- data.frame(tr_final$node.label)



genelist <- c(colnames(data_genes_w)[-1])

for (x in genelist) {

    # prepare data

    aaa <- data.frame(select(data_genes_w, as.name(x)))
    rownames(aaa) <- c(data_genes_w$hordeum_species)
    aaa <- mutate_if(aaa, is.numeric, as.factor)

    bbb<-setNames(aaa[,x],rownames(aaa))


    # fit the model
    # fitting the model

    fitER<-fitMk(tr_final,bbb,model="ARD")


    # reconstruction of the ancestral states
    # joint ancestral states

    fit.joint<-ancr(fitER,type="joint")

    # get the ancestral state probabilities

    ccc <- data.frame(fit.joint$ace)
    colnames(ccc) <- c(as.name(x))
    rownames(ccc) <- c(rownames(fit.joint$ace))


    df_joint <- cbind(df_joint, ccc)
    rownames(df_joint) <- c(rownames(ccc))


    # plot the tree

    cols <- setNames(brewer.pal(n = 2, name = "Dark2"), levels(bbb))

    pdf(paste0(JOINT,"/",x,"_joint.pdf"), paper = "a4r")
    plotTree.datamatrix(tr_final,as.data.frame(bbb),
        colors=list(cols),header=FALSE, fsize=1)
    legend("topleft",legend=levels(bbb),pch=22,title=x,
        pt.cex=2.5,pt.bg=cols,bty="n",cex=1.5)
    nodelabels(pie=to.matrix(fit.joint$ace,levels(bbb)),
    piecol=cols,cex=0.6)
    dev.off()



    # reconstruction of the ancestral states
    # marginal ancestral states

    fit.marginal<-ancr(fitER,type="marginal")

    # get the ancestral state probabilities

    eee <- as.data.frame(fit.marginal$ace)
    colnames(eee) <- c("zero", "one")
    rownames(eee) <- c(rownames(fit.marginal$ace))

    ddd <- data.frame(eee$zero)
    colnames(ddd) <- c(as.name(x))

    df_marginal <- cbind(df_marginal, ddd)
    rownames(df_marginal) <- c(rownames(eee))


    # plot the tree

    cols <- setNames(brewer.pal(n = 2, name = "Dark2"), levels(bbb))

    pdf(paste0(MARGINAL,"/",x,"_marginal.pdf"), paper = "a4r")
    plotTree.datamatrix(tr_final,as.data.frame(bbb),
        colors=list(cols),header=FALSE, fsize=1)
    legend("topleft",legend=levels(bbb),pch=22,title=x,
        pt.cex=2.5,pt.bg=cols,bty="n",cex=1.5)
    nodelabels(pie=fit.marginal$ace,piecol=cols,cex=0.6)
    dev.off()
    }

# export tables with the node probabilities

df_joint_ed <- df_joint[,-1]
df_joint_ed$node <- rownames(df_joint)

df_marginal_ed <- df_marginal[,-1]
df_marginal_ed$node <- rownames(df_marginal)


write.table(df_joint_ed, "241119_node_prob_joint.csv", sep ="\t", col.names = T, row.names = F)

write.table(df_marginal_ed, "241119_node_prob_marginal.csv", sep ="\t", col.names = T, row.names = F)







