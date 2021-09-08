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
                                    output.file = file.path(data.dir, 'check'))

nitrosamines <- cbind(sbs.cat.37$catSBS96,sbs.cat.38$catSBS96)

write.csv(nitrosamines, file='CombinedChr37_38vcfs.csv')
