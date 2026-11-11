# Programming assignment 1
# Arun 
# 2026/10/26
## =============================================================================

# Question 1 U2412
spec = read.table("Unknowns/U2412.txt")
mz = spec[,1]
int = spec[,2]
int = int/max(int)
# plot(mz,int,type="h")
bp_index = which(int==max(int))
bp_mz = mz[bp_index]
cat(paste0("A1. The mass-to-charge ratio for the base peak in the spectrum is ",bp_mz,".\n"))

# Question 2 U27629.txt
spec = read.table("Unknowns/U27629.txt")
mz = spec[,1]
int = spec[,2]
int = int/max(int)
plot(mz,int,type="h")
abline(h=0.25,lty=3)
POI_index = which(int >= 0.25)
POI_length = length(POI_index)
cat(paste0("A2. The number of peaks of interest is ", POI_length,".\n"))  


# Question 3 U14499.txt
spec = read.table("Unknowns/U14499.txt")
mz = spec[,1]
int = spec[,2]
int = int/max(int)
plot(mz,int,type="h")
abline(v=142,lty=3, col="red")
POI_index = which( mz >= 260)
POI_length = length(POI_index)
cat(paste0("A3. The number of peaks of interest is ", POI_length,".\n"))  
  
# Question 4 U11651.txt
spec = read.table("Unknowns/U11651.txt")
mz = spec[,1]
int = spec[,2]
int = int/max(int)
plot(mz,int,type="h")
abline(h=0.1,lty=3)
POI_index = which(int >= 0.1)
maxMz_POI = max(mz[POI_index])
cat(paste0("A4. The max mz of the peaks of interest is ", maxMz_POI,".\n"))  
  
# Question 5 - U21136.txt
spec = read.table("Unknowns/U21136.txt")
mz = spec[,1]
int = spec[,2]
int = int/max(int)
plot(mz,int,type="h")
abline(h=0.1,lty=3)
POI_index = which(int >= 0.1)
minMz_POI = min(mz[POI_index])
cat(paste0("A5. The min mz of the peaks of interest is ", minMz_POI,".\n"))  

# Question 6 - U35603.txt and U24255.txt
specA = read.table("Unknowns/U35603.txt")
mzA = specA[,1]
intA = specA[,2]
intA = intA/max(intA)

specB = read.table("Unknowns/U24255.txt")
mzB = specB[,1]
intB = specB[,2]/max(specB[,2])

plot(mzA,intA,type="h",ylim=c(-1,1),xlim=c(min(c(mzA,mzB)),max(c(mzA,mzB))),col="purple")
matplot(mzB,-intB,type="h",col="blue",add=TRUE)
abline(h=.10,lty=3)
abline(h=-.10,lty=3)

POI_A_index = which(intA>=0.1)
POI_B_index = which(intB>=0.1)

mz_POI_A = mzA[POI_A_index]
mz_POI_B = mzB[POI_B_index]

overlap_mz = intersect(mz_POI_A,mz_POI_B)
max_overlap_mz = max(overlap_mz)
cat(paste0("A6. The max mz of common peaks is ",max_overlap_mz,".\n"))

# Question 7 - U9200.txt and U2356.txt
specA = read.table("Unknowns/U9200.txt")
mzA = specA[,1]
intA = specA[,2]
intA = intA/max(intA)

specB = read.table("Unknowns/U2356.txt")
mzB = specB[,1]
intB = specB[,2]/max(specB[,2])

plot(mzA,intA,type="h",ylim=c(-1,1),xlim=c(min(c(mzA,mzB)),max(c(mzA,mzB))),col="purple")
matplot(mzB,-intB,type="h",col="blue",add=TRUE)
abline(h=.10,lty=3)
abline(h=-.10,lty=3)

POI_A_index = which(intA>=0.1)
POI_B_index = which(intB>=0.1)

mz_POI_A = mzA[POI_A_index]
mz_POI_B = mzB[POI_B_index]

overlap_mz = intersect(mz_POI_A,mz_POI_B)
min_overlap_mz = min(overlap_mz)
cat(paste0("A7. The min mz of common peaks is ",min_overlap_mz,".\n"))


# Question 8 U32380.txt and U22062.txt Euclidean Distance = sqrt(sum((A-B)^2))
specA = read.table("Unknowns/U32380.txt")
mzA = specA[,1]
intA = specA[,2]/max(specA[,2])

specB = read.table("Unknowns/U22062.txt")
mzB = specB[,1]
intB = specB[,2]/max(specB[,2])

A = numeric(max(c(mzA,mzB)))
B = numeric(max(c(mzA,mzB)))

for(i in 1:length(mzA)){
  A[mzA[i]] = intA[i]
}

for(i in 1:length(mzB)){
  B[mzB[i]] = intB[i]
}

EucDist = sqrt(sum((A-B)^2))
EucDist_rounded = round(EucDist,3)

cat(paste0("A8. The euclidean distance between the spectra is ",EucDist_rounded,".\n"))


# Question 9. U17894.txt and U35582.txt.
specA = read.table("Unknowns/U17894.txt")
mzA = specA[,1]
intA = specA[,2]/max(specA[,2])

specB = read.table("Unknowns/U35582.txt")
mzB = specB[,1]
intB = specB[,2]/max(specB[,2])

A = numeric(max(c(mzA,mzB)))
B = numeric(max(c(mzA,mzB)))

for(i in 1:length(mzA)){
  if(intA[i]>=.10){
    A[mzA[i]] = intA[i]
  }
}

for(i in 1:length(mzB)){
  if(intB[i]>=.10){
    B[mzB[i]] = intB[i]  
  }
}

EucDist = sqrt(sum((A-B)^2))
EucDist_rounded = round(EucDist,3)

cat(paste0("A9. The euclidean distance between the spectra using only significant peaks is ",EucDist_rounded,".\n"))


# Question 10. U1037.txt and U25688.txt
specA = read.table("Unknowns/U1037.txt")
mzA = specA[,1]
intA = specA[,2]/max(specA[,2])

specB = read.table("Unknowns/U25688.txt")
mzB = specB[,1]
intB = specB[,2]/max(specB[,2])

A = numeric(max(c(mzA,mzB)))
B = numeric(max(c(mzA,mzB)))

for(i in 1:length(mzA)){
  if(intA[i]>=.10){
    A[mzA[i]] = intA[i]
  }
}

for(i in 1:length(mzB)){
  if(intB[i]>=.10){
    B[mzB[i]] = intB[i]  
  }
}

CS_num = sum(A*B)
CS_den1 = sqrt(sum(A*A))
CS_den2 = sqrt(sum(B*B))

CS = CS_num/(CS_den1*CS_den2)
CS_rounded = round(CS,3)
cat("A10. The cosine similarity between the spectra using only significant peaks is ",CS_rounded,".\n")
