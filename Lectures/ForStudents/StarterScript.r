## FRSC 4730H - Programming Assignment 2 - Helper Script
##
## Developed by Arun Moorthy (arunmoorthy@trentu.ca)
## Last updated: September 7, 2026
## =============================================================================

## Step 0. Clear up the work space, load custom functions, load the reference mass spectral library
rm(list=ls()) # clears up the work space
source("CustomFunctions/LibrarySearchFunctions.R") # loads all the functions in the file
Library = readRDS("asm_ReferenceMSLibrary_SWDRUG314.RDS") # loads the reference library


## Example 1. A standard library search without any pre-processing of the unknown spectrum

Unknown = read.table("Unknowns/U1.txt")  # reads in the unknown mass spectrum
Metric = "Euclidean"                     # this can take one of two values "Euclidean" or "CosineSim"
nHits = 10                               # this is the number of most similar library entries we want returned

Results = LibrarySearch(Unknown,Library,Metric,nHits) # This returns the library search HitList and RunTime. 

View(Results$HitList)
cat(Results$RunTime) 

## Example 1a. Visualize a mass spectral comparison -

spec1_mz = Unknown[,1]
spec1_int = Unknown[,2]/max(Unknown[,2]) #normalize

index = which(Library$IndexData$Name=="N,N-diacetylnorephedrine")
indices = which(Library$SpectralData$SpecID==index)
spec2_mz = Library$SpectralData$mz[indices]
spec2_int = Library$SpectralData$int[indices]
spec2_int = spec2_int/max(spec2_int)

HT_plot(spec1_mz,spec1_int,spec2_mz,spec2_int)


## Example 2. A library search where the low intensity peaks from the unknown have been removed

Unknown = read.table("Unknowns/U2.txt")  # reads in the unknown mass spectrum
Metric = "CosineSim"                     # this can take one of two values "Euclidean" or "CosineSim"
nHits = 4                               # this is the number of most similar library entries we want returned

  # Normalize Unknown
  Unknown[,2] = Unknown[,2]/max(Unknown[,2])
  
  # Remove peaks with relative intensity less than 0.15 (15%)
  IndicesOfInterest = which(Unknown[,2]<0.15) # finds indices of peaks with relative intensity less than 0.15
  if(length(IndicesOfInterest)>0){
    Unknown = Unknown[-IndicesOfInterest,] # removes peaks in IndicesOfInterest
  }
  
Results = LibrarySearch(Unknown,Library,Metric,nHits) # This returns the library search HitList and RunTime. 

View(Results$HitList)
cat(Results$RunTime) 


## Example 3. A library search where only a subset of the library is searched against
Unknown = read.table("Unknowns/U3.txt")  # reads in the unknown mass spectrum
Metric = "Euclidean"                     # this can take one of two values "Euclidean" or "CosineSim"
nHits = 10                               # this is the number of most similar library entries we want returned
LibraryTWO = Library                    # this is a local copy of the library

  # Identify Library entries with molecular weight between 200 and 300, inclusive
  a1 = which(LibraryTWO$IndexData$MW>=200) # indices for entries greater than or equal 200
  a2 = which(LibraryTWO$IndexData$MW<=300) # indices for entries less than or equal to 300
  IndicesOfInterest = intersect(a1,a2) # intersection of indices
  
  if(length(IndicesOfInterest)>0){
    LibraryTWO$IndexData = LibraryTWO$IndexData[IndicesOfInterest,] # we're only keeping the library entries that meet the conditions
  }
  
Results = LibrarySearch(Unknown,LibraryTWO,Metric,nHits) # This returns the library search HitList and RunTime. 

View(Results$HitList)
cat(Results$RunTime) 



## Example 4. Pull a subset of library entries and construct a distance matrix, then visualize with MDS and cluster with agglomerative hierarchical clustering
U = seq(1,30)
D = array(0,dim=c(length(U),length(U)))

for(i in 1:length(U)){
  indices_i = which(Library$SpectralData$SpecID==U[i])
  spec_mz_i = Library$SpectralData$mz[indices_i]
  spec_int = Library$SpectralData$int[indices_i]
  spec_int_i = spec_int/max(spec_int)
  
  for(j in 1:length(U)){
    
    indices_j = which(Library$SpectralData$SpecID==U[j])
    spec_mz_j = Library$SpectralData$mz[indices_j]
    spec_int = Library$SpectralData$int[indices_j]
    spec_int_j = spec_int/max(spec_int)

    maxmz = max(spec_mz_i,spec_mz_j)
    A = spec2vec(spec_mz_i,spec_int_i,maxmz)
    B = spec2vec(spec_mz_j,spec_int_j,maxmz)
    
    D[i,j] = 1-cosineSim(A,B)
  }

}

dendrogram = hclust(as.dist(D))

mds = cmdscale(D,k=2)
rownames(mds) = U
plot(mds,pch=19,cex=0)
text(mds[,1],mds[,2],U)


U = seq(1,30)
data = array(0,dim=c(length(U),1000))
for(i in 1:length(U)){
  indices = which(Library$SpectralData$SpecID==i)
  mz = Library$SpectralData$mz[indices]
  int = Library$SpectralData$int[indices]
  data[i,] = spec2vec(mz,int,1000)
}

zeros = which(colSums(data)==0)
data = data[,-zeros]
pca_results = prcomp(data)
scores = pca_results$x[,1:2,drop=FALSE]
plot(scores,cex=0)
text(scores[,1],scores[,2],U)
