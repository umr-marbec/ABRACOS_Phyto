###############################################################################################
# Ecology > Functional Ecology > Functionnal Diversity (FD) index calculation
#
# Calculation of FD according to the Oikos paper :
#
# Mouchet M., Guilhaumon F., Villéger S., Mason N., Tomasini J. A., Mouillot D.
# Towards a consensus for calculating dendrogram-based functional diversity indices.
# 2008. OIKOS (Volume 117 Issue 5, Pages 794 - 800)
#
# Code by François Guilhaumon (francois.guilhaumon@univ-montp2.fr)
# Xtree, Getlength.inner and Getlength functions by Owen Petchey (o.petchey@sheffield.ac.uk) 
# 
# Please reffer to upcited paper as appropriate
##############################################################################################


##############################################################################################
#  UPDATES
#
# 7 april 2008 : added $ultram to the result list
#
##############################################################################################

GFD <- function(mat,comp=NA,distances=NA,dissimMethod="spectral",methods=NA,verb=FALSE) {

################################################################
# Functional Diversity (FD) index calculation                  #
#                                                              #
# Required packages : cluster, clue, gtools                    #
#                                                              #
# Arguments :                                                  #
# mat : a sXt species/traits matrix (with species labels for   #
# communities FD calculations)                                 #
# comp : a composition data.frame  of the form                 #
# cbind(communities_labels,species_labels)                     #
# distances : distance(s) to use (NA means all)                #
# methods : clustering method(s) to use (NA means all)         #
# dissimMethod : how dissimilarities are calculated (          #
#        argument 'method' of the function cl_dissimilarity    #
# (package clue))                                              #
# verb : should verbosity be activated ?                       #
################################################################

require(gtools)
require(cluster)
require(clue)


########################################################
#    Xtree and GetLength function by Owen Petchey      #
########################################################
Xtree <- function(tree){

 species.names <- tree$labels
  
  merge <- tree$merge
  height <- tree$height
  
  s <- nrow(merge)+1
  nodeid <- matrix(integer(s*(s-1)), s, s-1)
  height1 <- matrix(double(s*(s-1)), s, s-1)
  branch <- matrix(double(s*(s-1)), s, s-1)
  branchid <- matrix(character(s*(s-1)), s, s-1)
  
  for(i in (s-1):1)
    for(j in 1:2){
      goto <- matrix(logical(i*2), i, 2)
      goto[i,j] <- 1;
      for(k in i:1)
        for(l in 1:2){
          speciescount <- 0
          if(goto[k,l] == 1){
            if(merge[k,l] > 0){
              goto[merge[k,l],1] <- 1
              goto[merge[k,l],2] <- 1
            }
            if(merge[k,l] < 0){
              nodeid[-merge[k,l],i] <- i
              height1[-merge[k,l],i] <- height[i]
            }
          }
        }
    }  

  for(i in 1:s){
    last <- 0
    for(j in 1:(s-1)){
      if(nodeid[i,j] > 0){
        branch[i,j] <- height1[i,j] - last
        last <- height1[i,j]
      }
    }
  }

  for(n in 1:(s-1))
    for(i in 1:s){
      if(nodeid[i,n] == 0)
        branchid[i,n] <- paste("0", sep="-")
      if(nodeid[i,n] > 0){
        if(n == 1)
          branchid[i,n] <- paste(i, "0", nodeid[i,n], sep="-")
        if(n > 1){                                
          j <- 1
          while(j < n && nodeid[i,n-j] == 0)
            j <- j+1
          if(j == n)
            j <- j - 1
          if(nodeid[i,n-j] == 0)
            branchid[i,n] <- paste(i, nodeid[i,n-j], nodeid[i,n], sep="-")
          else
           branchid[i,n] <- paste(nodeid[i,n-j], nodeid[i,n], sep="-")
        }
      }
    }

  unibranch <- unique(branchid[branchid != "0"])
  
  

  H1 <- matrix(integer(s*(2*s-2)), s, 2*s-2)
  h2.prime <- double(2*s-2)

  for(i in 1:s)
    for(j in 1:(s-1))
      for(k in 1:(2*s-2))  
        if(unibranch[k] == branchid[i,j])
          h2.prime[k] <- branch[i,j]

  names(h2.prime) <- unibranch

  for(i in 1:s)
    for(j in 1:(s-1))
      for(k in 1:(2*s-2))
        if(branchid[i,j] == unibranch[k])
           H1[i,k] <- 1

  dimnames(H1) <- list(species.names,unibranch)  
  
  xtree <- list("h2.prime"=h2.prime, "H1"=H1)
  class(xtree) <- "xtree"
  xtree
  }

Getlength.inner <- function(xtree){
	if(!is.matrix(xtree[[2]])) result <- 0
	if(is.matrix(xtree[[2]])) result = sum(xtree[[1]][colSums(xtree[[2]]) != 0 & colSums(xtree[[2]]) < length(xtree[[2]][,1])])
	result
}

Getlength <- function(xtree, comp=NA){
	if(!is.data.frame(comp)) result <- Getlength.inner(xtree)
	if(is.data.frame(comp)){
		S <- tapply(comp[,2], comp[,1], function(x) length(x))
		FD <- tapply(comp[,2], comp[,1], function(x) Getlength.inner(list(xtree[[1]], xtree[[2]][!is.na(match(dimnames(xtree[[2]])[[1]], x)),])))
		#FD <- FD/Getlength.inner(xtree)
		result <- data.frame(S=S, FD=FD)
	}
result
}

#########################################
#End of functions coded by Owen Petchey #
#########################################

options(warn=-1)


#Vector of available distances
supported.dists = c("gower", "euclidean")
#Vector of available methods (as named in the clue package) named as in the paper
supported.methods = c(ward="ward",single="single",complete="complete",UPGMA="average",UPGMC="centroid",WPGMC="median",WPGMA="mcquitty")

#Arguments tests, distances and methods selection
if(!is.na(distances)) {
	if(prod(is.element(distances,supported.dists))!=1) stop("distances should be in 'gower' or 'euclidean'")
	Distances = distances
}else{
	Distances = supported.dists
}

if (verb) cat("#################################################\n# Functional Diversity (FD) index calculation\n")


if (verb) cat("#--\n# Distances are : ",Distances,"\n")

if(!is.na(methods)) {
	if(prod(is.element(methods,names(supported.methods)))!=1) stop("methods should be in 'ward','single','complete','UPGMA','UPGMC','WPGMC' or 'WPGMA'")
	Methods = supported.methods[is.element(names(supported.methods),methods)]
}else{
	Methods = supported.methods
}

if (verb) cat("# Methods are : ",names(Methods),"\n")

if(!is.na(comp)) {
	if(!is.data.frame(comp)) stop("'communities compositions' should be a data.frame ")
	comp=comp
}#end of if(!is.na(comp))

#number of methods
nbMethods = length(Methods)

#list of method combinations
methodList = lapply(seq(1:nbMethods),function(index) combinations(nbMethods,index,Methods))
names(methodList) = paste("Comb",seq(1:length(methodList)),sep="_")

#scaling the trait matrix id numeric
if(is.numeric(mat)) {mat=scale(mat)}

#Creating a list of distances matrix
dMatList = sapply(Distances,function(distance) daisy(mat,metric=distance,stand=F) , simplify=F )

#Creating a list containing all the trees for each distance
calcTreesForADistance = function(dist,...) {
	dmat = dMatList[[dist]]
	#growing trees
	res=sapply(Methods,function(met) hclust(dmat,met) ,  simplify = FALSE)
	#naming trees with arguments methods
	#names(res) = available.methods[match(Methods,methods.corresp)]
	res
}#end of calcTreesForADistance
if (verb) cat("#--\n# Simple trees calculation ... \n")
treeList = sapply(Distances,calcTreesForADistance,dMatList=dMatList ,simplify = FALSE)

if (verb) cat("# Simple trees calculated ! \n")

#Calculating dissimilarities for simple trees
calcDissimilarityForATree = function(tree,distance) cl_dissimilarity(cl_ultrametric(tree),dMatList[[distance]],method=dissimMethod)
calcDissimilaritiesForADistance = function(distance) lapply(treeList[[distance]],calcDissimilarityForATree,distance=distance)
dissimTrees = sapply(Distances,calcDissimilaritiesForADistance)

if (verb) cat("# Dissimilarities for simple trees calculated ! \n")

#Getting the best simple Tree
if(length(Methods)==1 && length(Distances)==1){
	bestTreeValue = dissimTrees[[1]]
	bestTreeMethod =  Methods
	bestTreeDistance = Distances
}else {
	if(length(Methods)==1) {
		aa = as.matrix(dissimTrees)
		bestTreeNum = which.min(aa)
		bestTreeValue = aa[bestTreeNum][[1]]
		bestTreeMethod =  names(Methods)
		bestTreeDistance = substr(rownames(aa)[bestTreeNum],1,nchar(rownames(aa)[bestTreeNum])-(nchar(names(Methods))+1))

	}else {
		bestTreeNum = which.min(dissimTrees)
		bestTreeValue = dissimTrees[bestTreeNum][[1]]
		bestTreeMethod =  rep(rownames(dissimTrees),length(Distances))[bestTreeNum]
		bestTreeDistance = t(replicate(length(Methods),colnames(dissimTrees)))[bestTreeNum] #rep(colnames(dissimTrees),length(Methods))
	}#end of if/else

}#end of if/else

#Calculating and processing the consensuses if needed
calConsensusForARow = function(row,distance){
	#extracting the methods to use in the consensus
	methods.in.consensus = names(supported.methods[is.element(supported.methods,row)])

	#constructing an ensemble of trees with the methods combination
	hens = cl_ensemble(list=lapply(paste("treeList$",distance,"$",methods.in.consensus,sep=""),function(x)eval(parse(text=x))))

	#naming the ensemble
	names(hens) = methods.in.consensus

	#constructing and returning the corresponding consensus
	cl_consensus(hens, method="majority", weights=1, control=list(p=1))

}#end of calConsensusForARow

if (length(Methods) > 1) {

	if (verb) cat("#--\n# Consensuses calculation ...\n")

	consMethodList = methodList[-1]
	calConsensusForADistance = function(distance) lapply(consMethodList,function(mat) apply(mat,1,calConsensusForARow,distance=distance))

	consensuses = sapply(Distances,calConsensusForADistance,simplify=F)

	if (verb) cat("# Consensuses calculated ! \n")

	#Calculating dissimilarities for consensuses using an ugly loop ...
	dissimCons = list()

	for (g in 1:length(Distances)) {

		dissimCons[[g]] = list()

		for (i in 1:length(consMethodList)) {

			dissimCons[[g]][[i]] = list()

			for (j in 1:dim(consMethodList[[i]])[1]) {

			dissimCons[[g]][[i]][[j]]  = calcDissimilarityForATree(consensuses[[g]][[i]][[j]],Distances[g])
	
			}#end of for j

		}#end of for i

	}#end of for g

	names(dissimCons) = Distances


	#Getting the best consensus

	minCons=Inf
	bestConsDistance = character()
	bestConsMethods = character()
	bestCons = NULL
	posConsMin=integer()
	for (i in 1:length(dissimCons)) {

		for (j in 1:length(dissimCons[[i]])) {

			for (k in 1:length(dissimCons[[i]][[j]])) {

				x = dissimCons[[i]][[j]][[k]]
				if (!is.na(x)) if(x < minCons) {minCons=x ; bestConsDistance = Distances[i] ; bestConsMethods = consMethodList[[j]][k,] ; bestCons = consensuses[[i]][[j]][[k]]}

			}#end of for k

		}#end of for j

	}#end of for i

if (verb) cat("# Dissimilarities for consensuses calculated ! \n")

}#end of if length(Methods) > 1



res1 = list()
cons=FALSE

if (length(Methods) > 1) {
	#chossing the best : simple tree or consensus
	if(minCons < bestTreeValue) {
		cons=TRUE
		res1 = list(dissim=minCons,distance=bestConsDistance,methods= bestConsMethods,tree=as.hclust(bestCons),ultram=cl_ultrametric(as.hclust(bestCons)),data=mat,comp=comp,cons=cons)

	}else {

		res1 = list(dissim=bestTreeValue,distance=bestTreeDistance,methods= bestTreeMethod,tree=treeList[[bestTreeDistance]][[bestTreeMethod]],ultram=cl_ultrametric(treeList[[bestTreeDistance]][[bestTreeMethod]]),data=mat,comp=comp,cons=cons)

	}#end of if minCons < bestTreeValue
} else {

	#the best is the best simple tree
	res1 = list(dissim=bestTreeValue,distance=bestTreeDistance,methods= bestTreeMethod,tree=treeList[[bestTreeDistance]][[names(bestTreeMethod)]],ultram=cl_ultrametric(treeList[[bestTreeDistance]][[names(bestTreeMethod)]]),data=mat,comp=comp,cons=cons)


}#end of if/else length(Methods) > 1

#Xtree for the best
bestXtree = Xtree(as.hclust(res1$tree))

# FD is the sum of the product of i.prime and h2.prime
FD <- Getlength(bestXtree,comp=comp)

res2 = c(res1,FD=FD,valid=TRUE)

options(warn=0)



if (dissimMethod=="cophenetic"){

	res2$limit <- .36

	if (res1$dissim > .36) { res2$valid=FALSE ; cat("#-----------------------------------\n# Warning : Dissimilarity between initial distance and ultrametric distance matrices is too high ! \n# Be carrefull in FD interpretation !\n")
	}
}

if (dissimMethod=="spectral"){

	sigma2 <- var(dMatList[[res1$distance]]) + var(cl_ultrametric(res1$tree))
	threshold <- 2*sqrt(sigma2)*sqrt(dim(mat)[1])

	res2$limit <- threshold
	
	if(res1$dissim > threshold) {res2$valid=FALSE ; cat("#-----------------------------------\n# Warning : Dissimilarity between initial distance and ultrametric distance matrices is too high ! \n# Be carrefull in FD interpretation !\n")

	}
}


if(verb) cat("#--\n# Best tree as a dissimilarity value of : ",res1$dissim,",\n")
if(verb) if(cons) cat("# this resulting from a consensus,\n")
if (verb) cat("# it is constructed with the ",res1$distance," distance\n")
if (verb) cat("# and the following clustering method(s) : ",unlist(res1$methods),"\n")
if(verb) if(!is.data.frame(comp)) {cat("#--\n# FD : ",FD,"\n")}else {cat("#--\n# FD : ",res2$FD.FD,"\n") }
if(verb) cat("#########################################################\n")


invisible(res2)




}#end of GFD function

