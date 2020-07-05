

plotCorrelTable<-function(correltable,fileName="CorrelTable.pdf",
                          
                          main="No title was provided",Colv=NULL,Rowv=NULL,
                          
                          cellnote=F,key=F){
  
  cairo_pdf(fileName, width=11, height=11, onefile=T)
  
  par(mfrow=c(1,1), oma=c(3,3,3,3), mar=c(3,3,3,3),cex.main=0.7)
  
  library(gplots)
  
  cellNote<-round(correltable,3)
  
  if(cellnote == F) cellNote<-matrix(data=NA,ncol=ncol(correltable),nrow=nrow(correltable))
  
  heatmap.2(x = correltable,Rowv = Rowv,Colv = Colv,trace = "none",key=key,
            
            cellnote=cellNote,notecol="black",notecex = 0.8,
            
            margins = c(6, 6),
            
            lmat=rbind( c(0, 3, 4), c(2,1,0) ),
            
            lwid=c(1, 12, 1.5),
            
            lhei=c(1, 10),
            
            main=main,
            
            cexRow= 0.7, cexCol = 0.7, col=redgreen(75),dendrogram = "none")
  
  dev.off()
  
}

