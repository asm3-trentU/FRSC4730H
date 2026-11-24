spec2vec <- function(mz,int,max_mz){
  vec = numeric(max_mz)
  for(i in 1:length(mz)){
    vec[mz[i]] = int[i] 
  }
  return(vec)
}

euclideanDistance <-function(x,y){
  if(length(x)!=length(y)){
    cat("\n\nVectors are not of equal length; cannot compute euclidean distance.\n\n")
    return(NA)
  }
  return(round(sqrt(sum((x-y)^2)),3))
}

cosineSim <-function(x,y){
  if(length(x)!=length(y)){
    cat("\n\nVectors are not of equal length; cannot compute cosine similarity.\n\n")
    return(NA)
  }
  
  num = sum(x*y)
  den1 = sqrt(sum(x*x))
  den2 = sqrt(sum(y*y))
  sim = num/(den1*den2+1e-8)
  return(round(sim,3))
}

LibrarySearch <- function(Unknown,Library,Metric,nHits){
  
  Unknown[,2] = Unknown[,2]/max(Unknown[,2]) # normalizes the relative intensity of the unknown mass spectrum
  nSearch = dim(Library$IndexData)[1]        # determines the number of library entries (i.e., number of comparisons)
  
  Name = character(nSearch)
  Score = numeric(nSearch)

  mzU = Unknown[,1]
  relU = Unknown[,2]
  
  stime <- system.time({
  #cat("Library Search in progress... ")
  for(i in 1:nSearch){
    index = Library$IndexData$ID[i]
    indices = which(Library$SpectralData$SpecID==index)
    mzL = Library$SpectralData$mz[indices]
    relL = Library$SpectralData$int[indices]
    relL = relL/max(relL) # make sure intensities are normalized
    
    max_mz = max_mz = max(max(mzU),max(mzL))
    vecU = spec2vec(mzU,relU,max_mz)
    vecL = spec2vec(mzL,relL,max_mz)
    
    Name[i] = Library$IndexData$Name[i]
    if(Metric=="Euclidean"){
      Score[i] = euclideanDistance(vecU,vecL)
    } else if (Metric=="CosineSim"){
      Score[i] = cosineSim(vecU,vecL)
    } else {
      cat("ERROR: Inappropriate Metric selection.\n\n")
      return(NULL)
    }
  }
  })[3]  
  #cat(paste("Library search complete! Elapsed time: ",round(stime,3)," seconds."))  
  
  if(Metric=="Euclidean"){
    Sort = order(Score,decreasing=FALSE)  
  } else if (Metric=="CosineSim") {
    Sort = order(Score,decreasing=TRUE)
  }
  
  Result = data.frame(Name = Name[Sort],Score = Score[Sort])
  return(list(HitList = Result[1:nHits,],
              RunTime = paste0("Elapsed time: ",round(stime,3)," seconds.")))
  
}

HT_plot <- function(spec1_mz,spec1_int,spec2_mz,spec2_int){
  xlim = c(min(spec1_mz,spec2_mz),max(spec1_mz,spec2_mz))
  plot(spec1_mz,spec1_int,col="blue",type="h",xlim=xlim,ylim=c(-1,1),
       yaxt="n",main = "Spec 1 (top) vs Spec 2 (bottom)",
       xlab = "mass-to-charge ratio",
       ylab = "relative intensity")
  matplot(spec2_mz,-spec2_int,col="red",type="h",add=TRUE)  

}