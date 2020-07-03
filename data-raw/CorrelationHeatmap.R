correlTable<-function(df,nonSmpCol="",type="cosine"){
  
  
  
  ## creates a matrix with cosine similarities between all samples in a catalog
  
  ## type can be either cosine, or any of the methods used by cor()
  
  
  
  if(!nonSmpCol[1] == ""){
    
    output<-matrix(data=NA,nrow=ncol(df)-length(nonSmpCol),ncol=ncol(df)-length(nonSmpCol))
    
    colnames(output)<-rownames(output)<-colnames(df)[-nonSmpCol]
    
    items<-colnames(df)[-nonSmpCol]
    
  } else {
    
    output<-matrix(data=NA,nrow=ncol(df),ncol=ncol(df))
    
    colnames(output)<-rownames(output)<-colnames(df)
    
    items<-colnames(df)
    
  }
  
  
  
  for(row in items){
    
    x<-df[,colnames(df) == row]
    
    if(!is.null(ncol(x))){x<-x[,1]}
    
    for(col in items){
      
      y<-df[,colnames(df) == col]
      
      if(!is.null(ncol(y))){y<-y[,1]}
      
      if(type=="cosine"){a<-cosine(x,y)} else {a<-cor(x,y,method = type)}
      
      output[row,col]<-a
      
    }
    
  }
  
  return(output)
  
}



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

