library(ggplot2)
library(factoextra)

draw_pca <- function(input_spectra){
  #draw pca for all samples
  combined_df <- as.data.frame(t(input_spectra))
  sbs96.pca <- prcomp(combined_df)
  df_sbs96.pca <- data.frame(names = row.names(sbs96.pca$x),sbs96.pca$x)
  rownames(df_sbs96.pca) <- NULL
  
  #add new columns
  df_sbs96.pca <- mutate(df_sbs96.pca,
                         cell_category = case_when(
                           startsWith(names,'HepG2') ~ 'HepG2',
                           startsWith(names,'MCF10A') ~ 'MCF10A')) %>%
    mutate(S9 = case_when(
      !grepl('S9',names) ~ 'FALSE',
      #since mcf10a is only treated with hamsterS9, S9 category is simply
      #set as with S9 or without
      #grepl('HamsterS9',names) ~ 'HamsterS9',
      TRUE ~ 'S9')) %>%
    mutate(Nitrosamine_type = case_when(
      grepl('NDEA',names) ~ 'NDEA',
      grepl('NDMA',names) ~ 'NDMA',
      grepl('NPIP',names) ~ 'NPIP',
      grepl('NPYR',names) ~ 'NPYR'
    ))
  
  df_sbs96.pca$CellTypeS9Status <- paste(df_sbs96.pca$S9,df_sbs96.pca$cell_category)
  
  p_combined <- ggplot(data = df_sbs96.pca) +
    geom_point(mapping = aes(x=PC1,y=PC2,color=Nitrosamine_type,shape = CellTypeS9Status,size = 1))+
    #shape 1: HepG2 w/o S9;
    #shape 2: MCF10A w/o S9;
    #shape 16: HepG2 w/ S9;
    #shape 17: MCF10A w/ S9
    scale_shape_manual(values=c(1,2,16,17))+
    ggtitle('Combined_PCA')+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    guides(size='none')#hide legend 'size'
  return(p_combined)
}
if(F){
  draw_pca(nitrosamines)
}

check_pca_quality <- function(input_spectra){
  combined_df <- as.data.frame(t(input_spectra))
  sbs96.pca <- prcomp(combined_df)
  bar_plot <- fviz_eig(sbs96.pca)
  rep_plot <- fviz_pca_var(sbs96.pca,
                           col.var = "contrib", # Color by contributions to the PC
                           gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                           repel = TRUE     # Avoid text overlapping
                           )
  return(list(bar_plot,rep_plot))
}

if(F){#draw pca for all samples
  combined_df <- as.data.frame(t(nitrosamines));dim(combined_df)
  sbs96.pca <- prcomp(combined_df);
  df_sbs96.pca <- data.frame(names = row.names(sbs96.pca$x),sbs96.pca$x); df_sbs96.pca;
  rownames(df_sbs96.pca) <- NULL;
  
  #add new columns
  df_sbs96.pca <- mutate(df_sbs96.pca,
                         cell_category = case_when(
                           startsWith(names,'HepG2') ~ 'HepG2',
                           startsWith(names,'MCF10A') ~ 'MCF10A')) %>%
    mutate(S9 = case_when(
      !grepl('S9',names) ~ 'FALSE',
      #since mcf10a is only treated with hamsterS9, S9 category is simply
      #set as with S9 or without
      #grepl('HamsterS9',names) ~ 'HamsterS9',
      TRUE ~ 'S9')) %>%
    mutate(Nitrosamine_type = case_when(
      grepl('NDEA',names) ~ 'NDEA',
      grepl('NDMA',names) ~ 'NDMA',
      grepl('NPIP',names) ~ 'NPIP',
      grepl('NPYR',names) ~ 'NPYR'
    ))
  
  df_sbs96.pca$CellTypeS9Status <- paste(df_sbs96.pca$S9,df_sbs96.pca$cell_category)
  
  p_combined <- ggplot(data = df_sbs96.pca) +
    geom_point(mapping = aes(x=PC1,y=PC2,color=Nitrosamine_type,shape = CellTypeS9Status,size = 1))+
    #shape 1: HepG2 w/o S9;
    #shape 2: MCF10A w/o S9;
    #shape 16: HepG2 w/ S9;
    #shape 17: MCF10A w/ S9
    scale_shape_manual(values=c(1,2,16,17))+
    ggtitle('Combined_PCA')+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    guides(size='none')#hide legend 'size'
  p_combined
  
  #check pca quality
  fviz_eig(sbs96.pca)
  fviz_pca_var(sbs96.pca,
               col.var = "contrib", # Color by contributions to the PC
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  )
  #draw results in 3D
  library(plotly)
  
  fig <- plot_ly(df_sbs96.pca, x = ~PC1, y = ~PC2, z = ~Nitrosamine_type, 
                 color = ~cell_category,colors = c('#BF382A', '#0C4B8E'),
                 symbol = ~S9)
  fig <- fig %>% add_markers()
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                     yaxis = list(title = 'PC2'),
                                     7,                       zaxis = list(title = 'Nitrosamine_type')))
  
  fig
  
  
  #pca for individual nitrosamines
  
  #draw pca for ndma, modify the dataframe first
  ndma_df <- as.data.frame(t(nitrosamines[,grep('NDMA',colnames(nitrosamines))]));
  sbs96NDMA.pca <- prcomp(ndma_df);sbs96NDMA.pca;
  df_sbs96NDMA.pca <- data.frame(names = row.names(sbs96NDMA.pca$x),sbs96NDMA.pca$x); df_sbs96.pca;
  rownames(df_sbs96NDMA.pca) <- NULL;
  df_sbs96NDMA.pca <- mutate(df_sbs96NDMA.pca,
                             cell_category = case_when(
                               startsWith(names,'HepG2') ~ 'HepG2',
                               startsWith(names,'MCF10A') ~ 'MCF10A')) %>%
    mutate(S9 = case_when(
      !grepl('S9',names) ~ 'FALSE',
      #grepl('HamsterS9',names) ~ 'HamsterS9',
      TRUE ~ 'S9'
    ))
  
  p <- ggplot(data = df_sbs96NDMA.pca) +
    geom_point(mapping = aes(x=PC1,y=PC2,color = cell_category,shape = S9,size = 5))+
    ggtitle('NDMA_PCA')+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    guides(size='none')#hide legend 'size'
  p
  
  #draw pca for ndea, modify the dataframe first
  NDEA_df <- as.data.frame(t(nitrosamines[,grep('NDEA',colnames(nitrosamines))]));
  sbs96NDEA.pca <- prcomp(NDEA_df);
  df_sbs96NDEA.pca <- data.frame(names = row.names(sbs96NDEA.pca$x),sbs96NDEA.pca$x);
  rownames(df_sbs96NDEA.pca) <- NULL;
  df_sbs96NDEA.pca <- mutate(df_sbs96NDEA.pca,
                             cell_category = case_when(
                               startsWith(names,'HepG2') ~ 'HepG2',
                               startsWith(names,'MCF10A') ~ 'MCF10A')) %>%
    mutate(S9 = case_when(
      !grepl('S9',names) ~ 'FALSE',
      #grepl('HamsterS9',names) ~ 'HamsterS9',
      TRUE ~ 'S9'
    ))
  
  p <- ggplot(data = df_sbs96NDEA.pca) +
    geom_point(mapping = aes(x=PC1,y=PC2,color = cell_category,shape = S9,size = 5))+
    ggtitle('NDEA_PCA')+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    guides(size='none')#hide legend 'size'
  p
  
  #draw pca for npip, modify the dataframe first
  NPIP_df <- as.data.frame(t(nitrosamines[,grep('NPIP',colnames(nitrosamines))]));
  sbs96NPIP.pca <- prcomp(NPIP_df);
  df_sbs96NPIP.pca <- data.frame(names = row.names(sbs96NPIP.pca$x),sbs96NPIP.pca$x);
  rownames(df_sbs96NPIP.pca) <- NULL;
  df_sbs96NPIP.pca <- mutate(df_sbs96NPIP.pca,
                             cell_category = case_when(
                               startsWith(names,'HepG2') ~ 'HepG2',
                               startsWith(names,'MCF10A') ~ 'MCF10A')) %>%
    mutate(S9 = case_when(
      !grepl('S9',names) ~ 'FALSE',
      #grepl('HamsterS9',names) ~ 'HamsterS9',
      TRUE ~ 'S9'
    ))
  
  p <- ggplot(data = df_sbs96NPIP.pca) +
    geom_point(mapping = aes(x=PC1,y=PC2,color = cell_category,shape = S9,size = 5))+
    ggtitle('NPIP_PCA')+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    guides(size='none')#hide legend 'size'
  p
  
  #draw pca for npyr, modify the dataframe first
  NPYR_df <- as.data.frame(t(nitrosamines[,grep('NPYR',colnames(nitrosamines))]));
  sbs96NPYR.pca <- prcomp(NPYR_df);
  df_sbs96NPYR.pca <- data.frame(names = row.names(sbs96NPYR.pca$x),sbs96NPYR.pca$x);
  rownames(df_sbs96NPYR.pca) <- NULL;
  df_sbs96NPYR.pca <- mutate(df_sbs96NPYR.pca,
                             cell_category = case_when(
                               startsWith(names,'HepG2') ~ 'HepG2',
                               startsWith(names,'MCF10A') ~ 'MCF10A')) %>%
    mutate(S9 = case_when(
      !grepl('S9',names) ~ 'FALSE',
      #grepl('HamsterS9',names) ~ 'HamsterS9',
      TRUE ~ 'S9'
    ))
  
  p <- ggplot(data = df_sbs96NPYR.pca) +
    geom_point(mapping = aes(x=PC1,y=PC2,color = cell_category,shape = S9,size = 5))+
    ggtitle('NPYR_PCA')+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    guides(size='none')#hide legend 'size'
  p
}

