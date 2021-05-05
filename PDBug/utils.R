scatter_color <- c("#66C2A5","#FC8D62")
heatmap_color <- c("#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58")
sample_color <- c("#543005","#BF812D","#DFC27D","#F6E8C3","#35978F","#40004B","#9970AB","#E7D4E8","#5AAE61","#00441B","gray75" )

dummy.code <- function(fac){
  fac <- as.factor(fac)
  level <- levels(fac)
  
  D <- matrix(0,nrow = length(fac),ncol = length(level))
  
  fac_num <- as.numeric(fac)
  
  for(i in 1:length(fac)){
    
    D[i, fac_num[i]] <- 1
  }
  colnames(D) <- level
  D <- as.data.frame(D)
  
  return(D)
}

OTU_info <- function(ps){
  OTU_name <- matrix(paste("ASV",1:ncol(ps@otu_table),sep = "_"),ncol = 1)
  OTU_name <- cbind(OTU_name,colnames(ps@otu_table))
  OTU_name <- cbind(OTU_name,ps@tax_table)
  colnames(OTU_name) <- c("ASV_idx","Sequence",colnames(OTU_name)[3:8])
  rownames(OTU_name) <- 1:ncol(ps@otu_table)
  OTU_name <- OTU_name[,-2]
  return(OTU_name)
}

cutoff <- function(R, r0, mode){
  col_id <- colnames(R)
  row_id <- rownames(R)
  
  R_cut <- R
  R_cut[abs(R_cut) >= r0] = 1
  R_cut[abs(R_cut) < r0] = 0
  R_cut = R_cut * sign(R)
  
  idx_row <- which(rowSums(R_cut) > 0)
  idx_col <- which(colSums(R_cut) > 0)
  
  if(mode == "-1-0-1"){
    R <- R_cut[idx_row, idx_col]
    colnames(R) <- col_id[idx_col]
    rownames(R) <- row_id[idx_row]
  }else{
    R <- R[idx_row, idx_col]
    colnames(R) <- col_id[idx_col]
    rownames(R) <- row_id[idx_row]
  }
  return(R)
}

microtrans <- function(ps, trans, Level, sub_data, PD){
  if (trans == "Composition"){
    ps <- transform_sample_counts(ps, function(x){x / sum(x)})
  }else if(trans == "Log"){
    ps <- transform_sample_counts(ps, function(x) log(1 + x))
  }
  
  if(Level == "ASV"){
    Microbiome <- as.data.frame(as.matrix(ps@otu_table))
    temp_name <- paste("ASV",1:ncol(Microbiome),sep = "_")
    colnames(Microbiome) <- temp_name
  }else if(Level == "Genus"){
    ps <- tax_glom(ps,taxrank = "Genus")
    Microbiome <- as.data.frame(as.matrix(ps@otu_table))
    colnames(Microbiome) <- ps@tax_table[,6]
  }else if(Level == "Family"){
    ps <- tax_glom(ps,taxrank = "Family")
    Microbiome <- as.data.frame(as.matrix(ps@otu_table))
    colnames(Microbiome) <- ps@tax_table[,5]
  }else if(Level == "Order"){
    ps <- tax_glom(ps,taxrank = "Order")
    Microbiome <- as.data.frame(as.matrix(ps@otu_table))
    colnames(Microbiome) <- ps@tax_table[,4]
  }else if(Level == "Phylum"){
    ps <- tax_glom(ps,taxrank = "Phylum")
    Microbiome <- as.data.frame(as.matrix(ps@otu_table))
    colnames(Microbiome) <- ps@tax_table[,2]
  }
  
  if (sub_data == "PD"){
    Microbiome <- Microbiome[PD == "PD",]
  }
  return(Microbiome)
}

clinicaltrans <- function(clinical.df, vartype, sub_data, PD){
  
  Clinical <- matrix(0,nrow = nrow(clinical.df))
  
  for(i in 1:ncol(clinical.df)){
    if ((vartype[i] == "discrete")|(vartype[i] == "PD discrete")){
      if(length(levels(clinical.df[,i])) == 2){
        temp <- matrix(as.numeric(as.factor(clinical.df[,i])) - 1,ncol = 1)
        colnames(temp) <- names(vartype[i])
      }else{
        temp <- as.matrix(dummy.code(clinical.df[,i]))
        colnames(temp) <- paste(names(vartype)[i], colnames(temp),sep = ".")
      }
      
    }else if ((vartype[i] == "continuous")|(vartype[i] == "PD continuous")){
      temp <- matrix(as.numeric(as.character(clinical.df[,i])),ncol = 1)
      colnames(temp) <- names(vartype)[i]
    }
    Clinical <- cbind(Clinical,temp)
  }
  
  Clinical <- as.data.frame(Clinical)
  rownames(Clinical) <- rownames(clinical.df)
  Clinical <- Clinical[,-(vartype == "missing")]

  
  if (sub_data == "PD"){
    Clinical <- Clinical[PD == "PD",]
  }
  
  
  
  return(Clinical)
}

cal_corr <- function(Microbiome, Clinical, vartype, D1, D2, trans, corr.method, sub_data, PD){
 # if (sub_data == "PD"){
 #   Clinical <- Clinical[PD == 1,]
 #   Microbiome <- Microbiome[PD == 1,]
 #   Clinical <- Clinical[,-1]
 # }
  
  if ((D1 == "Microbiome")&(D2 == "Microbiome")){
    if(trans == "Composition"){
      R <- propr::perb(Microbiome)@matrix
    }else{
      R <- cor(Microbiome, method = corr.method, use = "pairwise.complete.obs")
    }
  }else if ((D1 == "Clinical")&(D2 == "Clinical")){
    print("1")
    R <- cor(Clinical, method = corr.method, use = "pairwise.complete.obs")
  }else if((D1 == "Microbiome")&(D2 == "Clinical")){
    R <- cor(Microbiome, Clinical, method = corr.method, use = "pairwise.complete.obs")
  }else if((D1 == "Clinical")&(D2 == "Microbiome")){
    R <- cor(Clinical, Microbiome, method = corr.method, use = "pairwise.complete.obs")
  }else{
    R <- cor(Microbiome, method = corr.method, use = "pairwise.complete.obs")
  }
  return(R)
}

plot_heat <- function(R, r0, mode){
  R[is.na(R)] = 0
  R0 <- cutoff(R,r0,mode)
  
  if(length(R0) == 0){
    return(ggplotly(ggplot()))
  }else{
    R0_melt <- reshape2::melt(R0)
    p_R0 <- ggplot(R0_melt) + 
      geom_tile(aes(x = Var1, y = Var2, fill = value)) + 
      scale_fill_gradientn(colors = brewer.pal(9,"YlGnBu")) + 
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
    
    return(ggplotly(p_R0, source = "heatmap"))
  }
}

plot_scatter <- function(Microbiome, Clinical, R, r0, mode, pt1, pt2, D1, D2, sub_data, PD){
  R[is.na(R)] = 0
  R0 <- cutoff(R,r0,mode)
  pt1 <- rownames(R0)[pt1]
  pt2 <- colnames(R0)[pt2]
  if (sub_data == "PD"){
    PD = PD[PD == "PD"]
  }
  
  if((D1 == "Microbiome") &(D2 == "Microbiome")){
    temp <- Microbiome[,c(pt1,pt2)]
    temp <- as.data.frame(temp)
    pvalue <- cor.test(temp[,1], temp[,2], method = "spearman")[["p.value"]]
    p_scatter <- ggplot(temp,aes(x = temp[,1], y = temp[,2], color = PD)) + 
      geom_point() + 
      labs(x = pt1, y = pt2, title = paste("r = ",round(R[pt1,pt2],digits = 2),"(",round(pvalue,digits = 2),")")) + 
      theme_bw() +
      scale_color_manual(values = scatter_color)
    return(ggplotly(p_scatter))
  }else if((D1 == "Clinical") &(D2 == "Microbiome")){
    temp <- cbind(Clinical[,pt1],Microbiome[,pt2])
    temp <- as.data.frame(temp)
    pvalue <- cor.test(temp[,1], temp[,2], method = "spearman")[["p.value"]]
    p_scatter <- ggplot(temp,aes(x = temp[,1], y = temp[,2], color = PD)) + 
      geom_point() + 
      labs(x = pt1, y = pt2, title = paste("r = ",round(R[pt1,pt2],digits = 2),"(",round(pvalue,digits = 2),")")) + 
      theme_bw() +
      scale_color_manual(values = scatter_color)
    return(ggplotly(p_scatter))
  }else if((D1 == "Clinical") &(D2 == "Clinical")){
    temp <- Clinical[,c(pt1,pt2)]
    temp <- as.data.frame(temp)
    pvalue <- cor.test(temp[,1], temp[,2], method = "spearman")[["p.value"]]
    p_scatter <- ggplot(temp,aes(x = temp[,1], y = temp[,2], color = PD)) + 
      geom_point() + 
      labs(x = pt1, y = pt2, title = paste("r = ",round(R[pt1,pt2],digits = 2),"(",round(pvalue,digits = 2),")")) + 
      theme_bw() +
      scale_color_manual(values = scatter_color)
    return(ggplotly(p_scatter))
  }else if((D1 == "Microbiome") &(D2 == "Clinical")){
    temp <- cbind(Microbiome[,pt1],Clinical[,pt2])
    temp <- as.data.frame(temp)
    pvalue <- cor.test(temp[,1], temp[,2], method = "spearman")[["p.value"]]
    p_scatter <- ggplot(temp,aes(x = temp[,1], y = temp[,2], color = PD)) + 
      geom_point() + 
      labs(x = pt1, y = pt2, title = paste("r = ",round(R[pt1,pt2],digits = 2),"(",round(pvalue,digits = 2),")")) + 
      theme_bw() +
      scale_color_manual(values = scatter_color)
    return(ggplotly(p_scatter))
  }
}

RelAbundChart <- function(psra, level = "Phylum", maxAbd = 10, minprop = 0.05, mode = 1){
  
  psra <- transform_sample_counts(psra, function(x){x / sum(x)})
  
  psra_glom <- tax_glom(psra, taxrank = level)
  psra_melt <- psmelt(psra_glom)
  level1 <- rlang::sym(level)
  psra_melt1 <- psra_melt %>% group_by(!!level1, PD, Sample)  %>% summarise(mean = mean(Abundance))
  psra_melt0 <- psra_melt %>% group_by(!!level1) %>% summarise(mean = mean(Abundance))
  
  if(mode == 2){
    level_sel <- as.vector(as.matrix(psra_melt0[psra_melt0[,"mean"]>=minprop,1]))
    psra_melt2 <- psra_melt1[as.character(as.vector(as.matrix(psra_melt1[,1]))) %in% level_sel,]
    psra_check <- psra_melt2 %>% group_by(PD) %>% summarise(mean = sum(mean))
    if(nrow(psra_melt2) <= nrow(psra_melt1)){
      psra_check1 <- psra_melt2 %>% group_by(Sample, PD) %>% summarise(mean = sum(mean))
      other <- psra_check1 %>% mutate(!!level1 := "Other taxa") %>% select(!!level1, everything())
      other[,"mean"] = 1 - other[,"mean"]
      other <- data.frame(other)
      psra_melt2 <- data.frame(psra_melt2)
      psra_melt1 <- rbind(psra_melt2,other)
    }
    
    
  }else if(mode == 1){
    psra_melt0 <- psra_melt0 %>% arrange(desc(mean))
    level_sel <- as.vector(as.matrix(psra_melt0[1:min(maxAbd,nrow(psra_melt0)),1]))
    psra_melt2 <- psra_melt1[as.character(as.vector(as.matrix(psra_melt1[,1]))) %in% level_sel,]
    psra_check <- psra_melt2 %>% group_by(PD) %>% summarise(mean = sum(mean))
    if(nrow(psra_melt2) <= nrow(psra_melt1)){
      psra_check1 <- psra_melt2 %>% group_by(Sample, PD) %>% summarise(mean = sum(mean))
      other <- psra_check1 %>% mutate(!!level1 := "Other taxa")
      other[,"mean"] = 1 - other[,"mean"]
      other <- data.frame(other)
      psra_melt2 <- data.frame(psra_melt2)
      psra_melt1 <- rbind(psra_melt2,other)
    }
    
  }
  else{
    psra_melt0 <- psra_melt0 %>% arrange(desc(mean))
    level_ind <- intersect(1:min(maxAbd,nrow(psra_melt0)),which(psra_melt0[,"mean"]>=minprop))
    level_sel <- as.vector(as.matrix(psra_melt0[level_ind,1]))
    psra_melt2 <- psra_melt1[as.character(as.vector(as.matrix(psra_melt1[,1]))) %in% level_sel,]
    psra_check <- psra_melt2 %>% group_by(PD) %>% summarise(mean = sum(mean))
    if(nrow(psra_melt2) <= nrow(psra_melt1)){
      psra_check1 <- psra_melt2 %>% group_by(Sample, PD) %>% summarise(mean = sum(mean))
      other <- psra_check1 %>% mutate(!!level1 := "Other taxa") %>% select(!!level1, everything())
      other[,"mean"] = 1 - other[,"mean"]
      other <- data.frame(other)
      psra_melt2 <- data.frame(psra_melt2)
      psra_melt1 <- rbind(psra_melt2,other)
    }
  }
  psra_melt1[,level] = factor(as.character(as.vector((as.matrix(psra_melt1[,level])))),levels = c(unique(as.character(as.vector((as.matrix(psra_melt2[,level]))))),"Other taxa"))
  return(psra_melt1) 
}

sample_plot <- function(data, title, order_tax, level = "Family"){
  level1 <- rlang::sym(level)
  data_order <- data[data[,1] == order_tax,]
  data_order$PD <- as.factor(data_order$PD)
  data_order <- data_order %>% group_by(PD) %>% dplyr::arrange(desc(mean),.by_group = T)
  sample_ord <- data_order$Sample
  data$Sample <- factor(data$Sample,levels = sample_ord)
  lev <- levels(factor(data[,1]))
  t <- lev[1]
  lev[lev == order_tax] = lev[1]
  lev[1] <- order_tax
  
  data[,1] <- factor(data[,1], levels = lev)
  
  p <- ggplot(data) +
    geom_bar(aes(x = Sample, y = mean, fill = !!level1), position = "stack", stat = "identity") +
    xlab("Sample") +
    ylab("Abundance") +
    scale_fill_manual(values = sample_color)  +
    ggtitle(title) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_blank(),
          strip.background = element_rect(fill = "white"),
          panel.grid = element_blank(),
          legend.text.align = 0,
          legend.key.height = unit(4, "mm"),
          legend.key.width = unit(4, "mm"),
          legend.text = element_text(size = 9),
          legend.margin = ggplot2::margin(r = 7, unit = "mm"))
  return(p)
}

