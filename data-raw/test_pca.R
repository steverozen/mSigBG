library(ggplot2)
n_df <- as.data.frame(t(nitrosamines))
cc <- sub(".*HepG2.*", "HepG2", rownames(n_df))
cc <- sub(".*MCF10A.*", "MCF10A", cc)
S9 <- grepl("S9", rownames(n_df))
pch <- rep(2, nrow(n_df))
pch[cc == "HepG2" & !S9] <- 1
pch[cc == "HepG2" & S9] <- 16
pch[cc == "MCF10A" & S9] <- 17
nn <- cbind(cc, pch, S9, n_df)

#draw pca for all nitrosamines, modify the dataframe first
all.pca <- prcomp(n_df)
plot(all.pca$x, pch = pch)
library(factoextra)
fviz_eig(all.pca)
viz_pca_var(all.pca,
                          col.var = "contrib", # Color by contributions to the PC
                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                          repel = TRUE     # Avoid text overlapping
            )

p <- ggplot(data = as.data.frame(all.pca)) +
  geom_point(mapping = aes(x=PC1,y=PC2,
                           # shape = pch,
                           size = 5))+
  scale_shape_manual(values = pch)+
  ggtitle('All nitrosmaines')+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(size='none')#hide legend 'size'


#theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
##      legend.key.height = unit(0.5, 'cm'), #change legend key height
#      legend.key.width = unit(0.5, 'cm'), #change legend key width
#      legend.text = element_text(size=5)) #change legend text font size
p




df_sbs96NDMA.pca <- data.frame(names = row.names(sbs96NDMA.pca$x),sbs96NDMA.pca$x); df_sbs96.pca;
rownames(df_sbs96NDMA.pca) <- NULL;
df_sbs96NDMA.pca <- mutate(df_sbs96NDMA.pca,
                           cell_category = case_when(
                             startsWith(names,'HepG2') ~ 'HepG2',
                             startsWith(names,'MCF10A') ~ 'MCF10A')) %>%
  mutate(S9 = case_when(
    !grepl('S9',names) ~ 'FALSE',
    grepl('HamsterS9',names) ~ 'HamsterS9',
    TRUE ~ 'S9'
  ))



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
    grepl('HamsterS9',names) ~ 'HamsterS9',
    TRUE ~ 'S9'
  ))

p <- ggplot(data = df_sbs96NDMA.pca) +
  geom_point(mapping = aes(x=PC1,y=PC2,color = cell_category,shape = S9,size = 5))+
  ggtitle('NDMA_PCA')+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(size='none')#hide legend 'size'
  #theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
  ##      legend.key.height = unit(0.5, 'cm'), #change legend key height
  #      legend.key.width = unit(0.5, 'cm'), #change legend key width
  #      legend.text = element_text(size=5)) #change legend text font size
p

#draw pca for all samples
combined_df <- as.data.frame(t(nitrosamines));dim(combined_df)
sbs96.pca <- prcomp(combined_df);
df_sbs96.pca <- data.frame(names = row.names(sbs96.pca$x),sbs96.pca$x); df_sbs96.pca;
rownames(df_sbs96.pca) <- NULL;

df_sbs96.pca$new_cate <- paste(S9,cc)
df_sbs96.pca <- mutate(df_sbs96.pca,
                       cell_category = case_when(
                         startsWith(names,'HepG2') ~ 'HepG2',
                         startsWith(names,'MCF10A') ~ 'MCF10A')) %>%
  mutate(S9 = case_when(
    !grepl('S9',names) ~ 'FALSE',
    grepl('HamsterS9',names) ~ 'HamsterS9',
    TRUE ~ 'S9')) %>%
  mutate(Nitrosamine_type = case_when(
    grepl('NDEA',names) ~ 'NDEA',
    grepl('NDMA',names) ~ 'NDMA',
    grepl('NPIP',names) ~ 'NPIP',
    grepl('NPYR',names) ~ 'NPYR'
  ))
library(plotly)

fig <- plot_ly(df_sbs96.pca, x = ~PC1, y = ~PC2, z = ~Nitrosamine_type, 
               color = ~cell_category,colors = c('#BF382A', '#0C4B8E'),
               symbol = ~S9)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                   yaxis = list(title = 'PC2'),
            7,                       zaxis = list(title = 'Nitrosamine_type')))

fig

p3d <- ggplot(data = df_sbs96.pca) +
  geom_point(mapping = aes(x=PC1,y=PC2,color=Nitrosamine_type,shape = new_cate,size = 1))+
  scale_shape_manual(values=c(1,2,16,17))+
  ggtitle('Combined_PCA')+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(size='none')#hide legend 'size'
#theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
##      legend.key.height = unit(0.5, 'cm'), #change legend key height
#      legend.key.width = unit(0.5, 'cm'), #change legend key width
#      legend.text = element_text(size=5)) #change legend text font size
p3d


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
    grepl('HamsterS9',names) ~ 'HamsterS9',
    TRUE ~ 'S9'
  ))

p <- ggplot(data = df_sbs96NDEA.pca) +
  geom_point(mapping = aes(x=PC1,y=PC2,color = cell_category,shape = S9,size = 5))+
  ggtitle('NDEA_PCA')+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(size='none')#hide legend 'size'
#theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
##      legend.key.height = unit(0.5, 'cm'), #change legend key height
#      legend.key.width = unit(0.5, 'cm'), #change legend key width
#      legend.text = element_text(size=5)) #change legend text font size
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
    grepl('HamsterS9',names) ~ 'HamsterS9',
    TRUE ~ 'S9'
  ))

p <- ggplot(data = df_sbs96NPIP.pca) +
  geom_point(mapping = aes(x=PC1,y=PC2,color = cell_category,shape = S9,size = 5))+
  ggtitle('NPIP_PCA')+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(size='none')#hide legend 'size'
#theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
##      legend.key.height = unit(0.5, 'cm'), #change legend key height
#      legend.key.width = unit(0.5, 'cm'), #change legend key width
#      legend.text = element_text(size=5)) #change legend text font size
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
    grepl('HamsterS9',names) ~ 'HamsterS9',
    TRUE ~ 'S9'
  ))

p <- ggplot(data = df_sbs96NPYR.pca) +
  geom_point(mapping = aes(x=PC1,y=PC2,color = cell_category,shape = S9,size = 5))+
  ggtitle('NPYR_PCA')+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(size='none')#hide legend 'size'
#theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
##      legend.key.height = unit(0.5, 'cm'), #change legend key height
#      legend.key.width = unit(0.5, 'cm'), #change legend key width
#      legend.text = element_text(size=5)) #change legend text font size
p

