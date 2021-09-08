library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(ggplot2)
data.dir <- file.path('data-raw','package.variable.source.data','nitrosamine-vcfs');
AllFiles <- dir(path=data.dir, pattern = 'SNV',full.names = TRUE); AllFiles;
chr37_files <- c("data-raw/package.variable.source.data/nitrosamine-vcfs/HepG2_NDEA_cl1_SNVintersect.vcf",                 
                "data-raw/package.variable.source.data/nitrosamine-vcfs/HepG2_NDEA_cl2_SNVintersect.vcf",                 
                "data-raw/package.variable.source.data/nitrosamine-vcfs/HepG2_NDEA_cl3_SNVintersect.vcf",                       
                "data-raw/package.variable.source.data/nitrosamine-vcfs/HepG2_NDMA_cl1_SNVintersect.vcf",                 
                "data-raw/package.variable.source.data/nitrosamine-vcfs/HepG2_NDMA_cl2_SNVintersect.vcf",                 
                "data-raw/package.variable.source.data/nitrosamine-vcfs/HepG2_NDMA_cl3_SNVintersect.vcf",                         
                "data-raw/package.variable.source.data/nitrosamine-vcfs/HepG2_NPIP_cl1_SNVintersect.vcf",                 
                "data-raw/package.variable.source.data/nitrosamine-vcfs/HepG2_NPIP_cl2_SNVintersect.vcf",                
                "data-raw/package.variable.source.data/nitrosamine-vcfs/HepG2_NPIP_cl3_SNVintersect.vcf",                        
                "data-raw/package.variable.source.data/nitrosamine-vcfs/HepG2_NPYR_cl1_SNVintersect.vcf",               
                "data-raw/package.variable.source.data/nitrosamine-vcfs/HepG2_NPYR_cl2_SNVintersect.vcf",                      
                "data-raw/package.variable.source.data/nitrosamine-vcfs/MCF10A_NDEA+HamsterS9_10uM_cl3_SNVintersect.vcf", 
                "data-raw/package.variable.source.data/nitrosamine-vcfs/MCF10A_NDEA+HamsterS9_10uM_cl5_SNVintersect.vcf", 
                "data-raw/package.variable.source.data/nitrosamine-vcfs/MCF10A_NDMA+HamsterS9_10uM_cl2_SNVintersect.vcf", 
                "data-raw/package.variable.source.data/nitrosamine-vcfs/MCF10A_NDMA+HamsterS9_10uM_cl3_SNVintersect.vcf",          
                "data-raw/package.variable.source.data/nitrosamine-vcfs/MCF10A_NPIP+HamsterS9_10uM_cl1_SNVintersect.vcf", 
                "data-raw/package.variable.source.data/nitrosamine-vcfs/MCF10A_NPIP+HamsterS9_10uM_cl10_SNVintersect.vcf",
                "data-raw/package.variable.source.data/nitrosamine-vcfs/MCF10A_NPYR+HamsterS9_10uM_cl1_SNVintersect.vcf", 
                "data-raw/package.variable.source.data/nitrosamine-vcfs/MCF10A_NPYR+HamsterS9_10uM_cl8_SNVintersect.vcf")
chr38_files <- c("data-raw/package.variable.source.data/nitrosamine-vcfs/HepG2_NDEA+S9_10uM_cl1_SNVintersect.vcf",
                 "data-raw/package.variable.source.data/nitrosamine-vcfs/HepG2_NDMA+S9_10uM_cl1_SNVintersect.vcf",
                 "data-raw/package.variable.source.data/nitrosamine-vcfs/HepG2_NPIP+S9_10uM_cl1_SNVintersect.vcf",
                 "data-raw/package.variable.source.data/nitrosamine-vcfs/HepG2_NPYR+S9_10uM_cl2_SNVintersect.vcf",
                 "data-raw/package.variable.source.data/nitrosamine-vcfs/MCF10A_NDEA_50uM_cl2_SNVintersect.vcf",   
                 "data-raw/package.variable.source.data/nitrosamine-vcfs/MCF10A_NPIP_50uM_cl1_SNVintersect.vcf");
names37 <- sub(".*((MCF10A|HepG2).(NDEA|NDMA|NPIP|NPYR).(.*S9)?.*([1-9])).*", "\\2_\\3_\\4_cl\\5", chr37_files, perl = TRUE);
names38 <- sub(".*((MCF10A|HepG2).(NDEA|NDMA|NPIP|NPYR).(.*S9)?.*([1-9])).*", "\\2_\\3_\\4_cl\\5", chr38_files, perl = TRUE);

sbs.cat.37 <- ICAMS::VCFsToCatalogs(files= chr37_files,
                                    variant.caller = 'strelka', 
                                    ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5,
                                    region = 'genome',
                                    names.of.VCFs = names37,
                                    output.file = file.path(data.dir, 'check'));

 
sbs.cat.38 <- ICAMS::VCFsToCatalogs(files= chr38_files,
                                    variant.caller = 'strelka', 
                                    ref.genome = BSgenome.Hsapiens.UCSC.hg38,
                                    region = 'genome',
                                    names.of.VCFs = names38,
                                    output.file = file.path(data.dir, 'check'));

nitrosamines <- cbind(sbs.cat.37$catSBS96,sbs.cat.38$catSBS96);

write.csv(nitrosamines, file='CombinedChr37_38vcfs.csv');

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
                                   zaxis = list(title = 'Nitrosamine_type')))

fig

p3d <- ggplot(data = df_sbs96.pca) +
  geom_point(mapping = aes(x=PC1,y=PC2,color=Nitrosamine_type,fill = cell_category,shape = S9,size = 3))+
  scale_fill_manual(values = c('red','black'))+
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

