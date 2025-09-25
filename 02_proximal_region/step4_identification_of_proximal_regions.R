library(data.table)
library(ggplot2)
library(cowplot)
library(GenomicRanges)

setwd("/Volumes/agruppen/dg7/fengj/panh2/chr4H")

list<-c("Hordeum_bogdanii","Hordeum_brevisubulatum","Hordeum_bulbosum_Hap1",
        "Hordeum_californicum","Hordeum_chilense","Hordeum_comosum","Hordeum_cordobense",
        "Hordeum_erectifolium","Hordeum_euclaston","Hordeum_flexuosum","Hordeum_gussoneanum",
        "Hordeum_intercedens","Hordeum_jubatum_Sub1","Hordeum_jubatum_Sub2",
        "Hordeum_marinum_BCC2001","Hordeum_murinum_BCC2017","Hordeum_murinum_Sub1",
        "Hordeum_murinum_Sub2","Hordeum_muticum","Hordeum_patagonicum","Hordeum_pubiflorum",
        "Hordeum_pusillum","Hordeum_roshevitzii","Hordeum_secalinum_Sub1","Hordeum_secalinum_Sub2",
        "Hordeum_stenostachys","Hordeum_vulgare")

### Determine the threshold value as 0.5

pdf("Normalized_coverage_distribution1.pdf",width = 5, height = 10)
for(ref_specie in list){
  finalcov1<-readRDS(paste("/Volumes/agruppen/dg7/fengj/panh2/chr4H/alignment/",ref_specie,"_cov1.RDS",sep = ""))
  print(ggplot(finalcov1, aes(x= ncov, y = qrysp)) +
          xlab("Coverage/Average Coverage")+ylab(ref_specie)+facet_grid(round(DivergenceTime,2)~.,scales = "free",space = "free")+
          geom_violin()+theme_minimal())
}
dev.off()

pdf("Normalized_coverage_distribution2.pdf",width = 20, height = 25)
for(ref_specie in list){
  finalcov1<-readRDS(paste("/Volumes/agruppen/dg7/fengj/panh2/chr4H/alignment/",ref_species,"_cov1.RDS",sep = ""))
  print(ggplot(finalcov1, aes(x= V2/1e6, y = ncov,fill=color)) + labs(fill="Cut-off")+
          xlab(paste(ref_species,"Position"))+ylab("Coverage/Average Coverage")+facet_grid(round(DivergenceTime,2)+gsub("_","\n",qrysp)~V1)+
          geom_bar(stat = "identity")+theme_minimal())
}
dev.off()


### identification of proximal regions boundary

setwd("/Volumes/agruppen/dg7/fengj/panh2/chr4H/paf/Zidentity_boundary")

### Select species separation time comparison 9 million years ago.

list_s9<-c("Hordeum_bogdanii..Hordeum_vulgare","Hordeum_brevisubulatum..Hordeum_vulgare",
        "Hordeum_bulbosum_Hap1..Hordeum_euclaston","Hordeum_californicum..Hordeum_vulgare",
        "Hordeum_chilense..Hordeum_vulgare","Hordeum_comosum..Hordeum_vulgare",
        "Hordeum_cordobense..Hordeum_vulgare","Hordeum_erectifolium..Hordeum_vulgare",
        "Hordeum_euclaston..Hordeum_vulgare","Hordeum_flexuosum..Hordeum_vulgare",
        "Hordeum_gussoneanum..Hordeum_vulgare","Hordeum_intercedens..Hordeum_vulgare",
        "Hordeum_jubatum_Sub1..Hordeum_vulgare","Hordeum_jubatum_Sub2..Hordeum_vulgare",
        "Hordeum_marinum_BCC2001..Hordeum_vulgare","Hordeum_murinum_BCC2017..Hordeum_vulgare",
        "Hordeum_murinum_Sub1..Hordeum_vulgare","Hordeum_murinum_Sub2..Hordeum_vulgare",
        "Hordeum_muticum..Hordeum_vulgare","Hordeum_patagonicum..Hordeum_vulgare",
        "Hordeum_pubiflorum..Hordeum_vulgare","Hordeum_pusillum..Hordeum_vulgare",
        "Hordeum_roshevitzii..Hordeum_vulgare","Hordeum_secalinum_Sub1..Hordeum_vulgare",
        "Hordeum_secalinum_Sub2..Hordeum_vulgare","Hordeum_stenostachys..Hordeum_vulgare",
        "Hordeum_vulgare..Hordeum_euclaston")
  

for(pair_l in list_s9){
  ref_species=unlist(strsplit(pair_l,"[..]"))[1]
  qry_species=unlist(strsplit(pair_l,"[..]"))[3]
  print(paste(ref_species,qry_species))
  cov<-fread(paste(pair_l,".ref.cov",sep = ""),header = F)
  cov[,ref:=ref_species]
  cov[,qry:=qry_species]
  meancov=mean(cov$V5)
  cov[,ncov:=V5/meancov]
  cov[,color:=">0.5"]
  cov[ncov<=0.5,color:="<=0.5"]
  bed_data <- cov
  threshold <- 0.5       # Threshold
  # Filter target area (using data.table for efficient filtering)
  filtered_data <- bed_data[ncov < threshold]
  # Handling empty result cases
  if (nrow(filtered_data) == 0) {
    stop("No region satisfying the conditions was identified.")
  }
  
  # Conversion to a GenomicRanges object (automatically handling coordinate transformation)
  gr <- GRanges(
    seqnames = filtered_data$V1,
    ranges = IRanges(
      start = filtered_data$V2 + 1,  
      end = filtered_data$V3
    )
  )
  
  # Merge intervals (while retaining chromosome grouping)
  merged_gr <- reduce(gr,min.gapwidth = 6e6)
  
  # Convert back to a data.table and adjust coordinates
  result <- data.table(
    chrom = as.character(seqnames(merged_gr)),
    start = start(merged_gr) - 1,  
    end = end(merged_gr)
  )
  
  # Sort by the original chromosome order
  setorder(result, chrom, start)->result
  result[end-start>10e6]-> result
  result[,species:=ref_species]
  
  paf<-read_paf(paste(pair_l,".paf",sep = ""),primary_only=T)
  setorder(paf,reference,reference_start)-> pafl
  pafl[mapq>=20 & reference==query]-> pafl
  pafl[,refsp:=ref_species]
  pafl[,qrysp:=qry_species]
  
  pdf(paste(ref_species,"_boundary.pdf",sep = ""),6,6)
  for(chr in c("chr1H","chr2H","chr3H","chr4H","chr5H","chr6H","chr7H")){
    # Reference coverage bar plot
    x_bar <- ggplot(cov[V1==chr], aes(x = V2/1e6, y = ncov,fill = color)) +
      geom_bar(stat = "identity") + xlab("")+ylab("NCoverage")+
      theme_minimal() + theme(legend.position = "null")
    # Main Figure: Scatter Plot
    main_plot <- ggplot(pafl[mapq>=20&reference==chr], aes(x = reference_start/1e6, y = query_start/1e6, xend = reference_end/1e6, yend = query_end/1e6)) +
      xlab("")+ylab("Alignment")+
      geom_vline(data = result[chrom==chr], aes(xintercept = start/1e6) ,linetype=4,linewidth=1,color="blue")+
      geom_vline(data = result[chrom==chr], aes(xintercept = end/1e6),linetype=4,linewidth=1,color="blue")+
      geom_segment(arrow = arrow(length = unit(0.01, "inches")))+
      theme_minimal()+theme(legend.position = "null")
    # Combined Figure
    print(plot_grid(
      x_bar,main_plot, 
      ncol = 1,rel_heights = c(2, 6)
    ))
  }
  dev.off()
  if(pair_l=="Hordeum_bogdanii..Hordeum_vulgare"){
    copy(result)->finalresult
  }else{
    rbind(finalresult,result)->finalresult
  }
  if(pair_l=="Hordeum_bogdanii..Hordeum_vulgare"){
    copy(pafl)->finalpafl
  }else{
    rbind(finalpafl,pafl)->finalpafl
  }
  if(pair_l=="Hordeum_bogdanii..Hordeum_vulgare"){
    copy(cov)->finalcov
  }else{
    rbind(finalcov,cov)->finalcov
  }
}

finalresult[,start:=start/1e6]
finalresult[,end:=end/1e6]
write.table(finalresult,"pericen_all.txt",sep = "\t",quote = F,row.names = F,col.names = F)

### Genome alignment plot and normalized coverage bar plot, with proximal region boundaries identified by manual inspection.

