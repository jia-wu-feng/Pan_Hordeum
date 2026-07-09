library(data.table)
library(ggplot2)
setwd("/Volumes/agruppen/dg7/fengj/panh2/chr4H")

read_paf<-function(file, primary_only=T, save=F){
  fread(head=F, cmd=paste("cat", file, "| cut -f -13"),
        col.names=c("query", "query_length", "query_start", "query_end",
                    "orientation", "reference", "reference_length",
                    "reference_start", "reference_end", "matches", "alnlen", "mapq", "type"))->z
  z[, c("query_start", "query_end", "reference_start", "reference_end") := list(query_start + 1, query_end + 1, reference_start + 1, reference_end + 1)]
  z[, orientation := as.integer(ifelse(orientation == "+", 1, -1))]
  z[, type := sub("tp:A:", "", type)]
  if(primary_only){
    z["P", on="type"]->z
  }
  if(save){
    saveRDS(z, file=sub("paf.gz$", "Rds", file))
  }
  z[]
}

list<-c("Hordeum_bogdanii","Hordeum_brevisubulatum","Hordeum_bulbosum_Hap1",
        "Hordeum_californicum","Hordeum_chilense","Hordeum_comosum","Hordeum_cordobense",
        "Hordeum_erectifolium","Hordeum_euclaston","Hordeum_flexuosum","Hordeum_gussoneanum",
        "Hordeum_intercedens","Hordeum_jubatum_Sub1","Hordeum_jubatum_Sub2",
        "Hordeum_marinum_BCC2001","Hordeum_murinum_BCC2017","Hordeum_murinum_Sub1",
        "Hordeum_murinum_Sub2","Hordeum_muticum","Hordeum_patagonicum","Hordeum_pubiflorum",
        "Hordeum_pusillum","Hordeum_roshevitzii","Hordeum_secalinum_Sub1","Hordeum_secalinum_Sub2",
        "Hordeum_stenostachys","Hordeum_vulgare","Hordeum_gussoneanum_Sub1","Hordeum_gussoneanum_Sub2",
        "Hordeum_parodii_Sub1","Hordeum_parodii_Sub2","Hordeum_parodii_Sub3")

### Read all alignment files based on the reference genome and store them

for(ref_specie in list){
  print(paste("Ref:",ref_specie))
  for(qry_specie in list){
    print(paste("Qry:",qry_specie))
    paf<-read_paf(paste("paf2/",ref_specie,"..",qry_specie,".paf",sep = ""),primary_only=T)
    mappingq=quantile(paf$mapq,0.80)
    mappinga=quantile(paf$alnlen,0.80)
    setorder(paf,reference,reference_start)-> pafl
    pafl[mapq>=mappingq & alnlen>=mappinga & reference==query]-> pafl
    pafl[,refsp:=ref_specie]
    pafl[,qrysp:=qry_specie]
    if(qry_specie=="Hordeum_bogdanii"){
      copy(pafl)->finalpaf
    }else{
      rbind(finalpaf,pafl)->finalpaf
    }
  }
  
  finalpaf[,timeid:=paste(refsp,qrysp)]
  
  fread("/Users/Feng/Nextcloud2/Hordeum Pangenome/3 review/tree fig/divergence_times.csv")-> divtime
  divtime[,timeid:=paste(Species1,Species2)]
  
  divtime[finalpaf,on=.(timeid)][refsp!=qrysp]-> finalpaf
  
  saveRDS(finalpaf,file=paste("/Volumes/agruppen/dg7/fengj/panh2/chr4H/alignment2/",ref_specie,"_alignment.RDS",sep = ""))
}

###Read coverage of reference/query genome alignments

for(ref_specie in list){
    print(paste("Ref:",ref_specie))
    for(qry_specie in list){
    print(paste("Qry:",qry_specie))
    cov1<-fread(paste("paf2/",ref_specie,"..",qry_specie,".ref.cov",sep = ""),header = F)
    cov1[,refsp:=ref_specie]
    cov1[,qrysp:=qry_specie]
    cov2<-fread(paste("paf2/",ref_specie,"..",qry_specie,".qry.cov",sep = ""),header = F)
    cov2[,refsp:=ref_specie]
    cov2[,qrysp:=qry_specie]
    cov1[,color:="b1"]
    if(qry_specie=="Hordeum_bogdanii"){
      copy(cov1)->finalcov1
      copy(cov2)->finalcov2
    }else{
      rbind(finalcov1,cov1)->finalcov1
      rbind(finalcov2,cov2)->finalcov2
    }
  }
  ### Mark species separation time
  fread("/Users/Feng/Nextcloud2/Hordeum Pangenome/3 review/tree fig/divergence_times.csv")-> divtime
  divtime[,timeid:=paste(Species1,Species2)]
  ### Normailzed reference coverage based on average coverage
  finalcov1[,timeid:=paste(refsp,qrysp)]
  divtime[finalcov1,on=.(timeid)][refsp!=qrysp]-> finalcov1
  finalcov1[,mean(V5),key=.(timeid)]->sumqry
  colnames(sumqry)<-c("timeid","meancov")
  sumqry[finalcov1,on=.(timeid)]-> finalcov1
  finalcov1[,ncov:=V5/meancov]
  finalcov1[,color:=">0.5"]
  finalcov1[ncov<=0.5,color:="<=0.5"]
  ### Normailzed query coverage based on average coverage
  finalcov2[,timeid:=paste(refsp,qrysp)]
  divtime[finalcov2,on=.(timeid)][refsp!=qrysp]-> finalcov2
  finalcov2[,mean(V5),key=.(timeid)]->sumqry
  colnames(sumqry)<-c("timeid","meancov")
  sumqry[finalcov2,on=.(timeid)]-> finalcov2
  finalcov2[,ncov:=V5/meancov]
  finalcov2[,color:=">0.5"]
  finalcov2[ncov<=0.5,color:="<=0.5"]
  
  saveRDS(finalcov1,file=paste("/Volumes/agruppen/dg7/fengj/panh2/chr4H/alignment2/",ref_specie,"_cov1.RDS",sep = ""))
  saveRDS(finalcov2,file=paste("/Volumes/agruppen/dg7/fengj/panh2/chr4H/alignment2/",ref_specie,"_cov2.RDS",sep = ""))
}
