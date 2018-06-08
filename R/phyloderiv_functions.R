#' Compute the Robinson-Foulds distance between pairs of trees at corresponding indices.
#'
#' @param est a list of estimated trees
#' @param true a list of true trees
#' @return a vector of RF distances between corresponding true and est trees.
compute.treedists.rf <- function(est, true) {
  if(length(est) != length(true)) {
    message("Families of trees are not the same size!")
  }
  dist <- rep(0,length(est))
  #For each tree in the list (both lists are the same size)
  for(i in 1:(length(est))) {
    currTree <- est[[i]] #ith tree from the estimated tree list
    nextTree <- true[[i]] #ith tree from the true tree list
    dist[i] <- RF.dist(currTree, nextTree)
  }
  return(dist)
}

#' Compute the path interval distance between pairs of trees at corresponding indices.
#'
#' @param est a list of estimated trees
#' @param true a list of true trees
#' @return a vector of path interval distances between corresponding true and est trees.
compute.treedists.intvl <- function(est, true) {
  if(length(est) != length(true)) {
    message("Families of trees are not the same size!")
  }

  #Function from treedist function
  coph <- function(x, path=FALSE){
    if (is.null(attr(x, "order")) || attr(x, "order") == "cladewise")
      x <- reorder(x, "postorder")
    nTips = as.integer(length(x$tip.label))
    parents = as.integer(x$edge[,1])
    kids = as.integer(x$edge[,2])
    lp= as.integer(length(parents))
    nNode = as.integer(x$Nnode)
    m = as.integer(max(x$edge))
    el = double(m)
    if(path) el <- rep(1.0, m)
    else el[kids] = x$edge.length
    dm <- .C("C_cophenetic", kids, parents, as.double(el), lp, m, nTips, nNode, double(nTips*(nTips-1L)/2L))[[8]]
    attr(dm, "Size") <- nTips
    attr(dm, "Labels") <- x$tip.label
    attr(dm, "Diag") <- FALSE
    attr(dm, "Upper") <- FALSE
    class(dm) <- "dist"
    dm
  }

  #Modified function from treedist function
  path.interval.dist <- function(tree1, tree2, check.labels=TRUE) {
    tree1 = unroot(tree1)
    tree2 = unroot(tree2)

    if (check.labels) {
      ind <- match(tree1$tip.label, tree2$tip.label)
      if (any(is.na(ind)) | length(tree1$tip.label) !=
          length(tree2$tip.label))
        stop("trees have different labels")
      tree2$tip.label <- tree2$tip.label[ind]
      ind2 <- match(1:length(ind), tree2$edge[, 2])
      tree2$edge[ind2, 2] <- order(ind)
    }

    tree1 = reorder(tree1, "postorder")
    tree2 = reorder(tree2, "postorder")

    if(!is.binary.tree(tree1) | !is.binary.tree(tree2))message("Trees are not binary!")

    l = length(tree1$tip.label)

    tree1$edge.length = rep(1, nrow(tree1$edge))
    tree2$edge.length = rep(1, nrow(tree2$edge))
    dt1 = coph(tree1)
    dt2 = coph(tree2)
    path.interval.distance = max(abs((dt1 - dt2)))
    return(path.interval.distance)
  }

  dists <- rep(0,length(est))
  #for each tree in the list (both lists are the same size)
  for(i in 1:length(est)) {
    tree1 = est[[i]] #ith tree from the estimated tree list
    tree2 = true[[i]] #ith tree from the true tree list

    path.interval.distance = path.interval.dist(tree1, tree2)
    dists[i] = path.interval.distance
  }
  return(dists)
}


#' Compute the phylogenetic derivative of a family of trees using Robinson-Foulds distance.
#'
#' @param trees a list of phylogenetic trees.
#' @return a vector with the phylogenetic derivative
compute.phyderiv.rf <- function(trees) {
  dist <- rep(0,length(trees)-1)
  for(i in 1:(length(trees)-1)) {
    currTree <- trees[[i]]
    nextTree <- trees[[(i+1)]]
    dist[i] <- RF.dist(currTree, nextTree)
  }
  return(dist)
}


#' Compute the phylogenetic derivative of a family of trees using the path interval distance.
#'
#' @param trees a list of phylogenetic trees.
#' @return a vector with the phylogenetic derivative
compute.phyderiv.intvl <- function(trees) {

  #Function from treedist function
  coph <- function(x, path=FALSE){
    if (is.null(attr(x, "order")) || attr(x, "order") == "cladewise")
      x <- reorder(x, "postorder")
    nTips = as.integer(length(x$tip.label))
    parents = as.integer(x$edge[,1])
    kids = as.integer(x$edge[,2])
    lp= as.integer(length(parents))
    nNode = as.integer(x$Nnode)
    m = as.integer(max(x$edge))
    el = double(m)
    if(path) el <- rep(1.0, m)
    else el[kids] = x$edge.length
    dm <- .C("C_cophenetic", kids, parents, as.double(el), lp, m, nTips, nNode, double(nTips*(nTips-1L)/2L))[[8]]
    attr(dm, "Size") <- nTips
    attr(dm, "Labels") <- x$tip.label
    attr(dm, "Diag") <- FALSE
    attr(dm, "Upper") <- FALSE
    class(dm) <- "dist"
    dm
  }

  #Modified function from treedist function
  path.interval.dist <- function(tree1, tree2, check.labels=TRUE) {
    tree1 = unroot(tree1)
    tree2 = unroot(tree2)

    if (check.labels) {
      ind <- match(tree1$tip.label, tree2$tip.label)
      if (any(is.na(ind)) | length(tree1$tip.label) !=
          length(tree2$tip.label))
        stop("trees have different labels")
      tree2$tip.label <- tree2$tip.label[ind]
      ind2 <- match(1:length(ind), tree2$edge[, 2])
      tree2$edge[ind2, 2] <- order(ind)
    }

    tree1 = reorder(tree1, "postorder")
    tree2 = reorder(tree2, "postorder")

    if(!is.binary.tree(tree1) | !is.binary.tree(tree2))message("Trees are not binary!")

    l = length(tree1$tip.label)

    tree1$edge.length = rep(1, nrow(tree1$edge))
    tree2$edge.length = rep(1, nrow(tree2$edge))
    dt1 = coph(tree1)
    dt2 = coph(tree2)
    path.interval.distance = max(abs((dt1 - dt2)))
    return(path.interval.distance)
  }

  dists <- rep(0,(length(trees)-1))
  for(i in 1:(length(trees)-1)) {
    tree1 = trees[[i]] #first tree
    tree2 = trees[[(i+1)]] #next tree
    path.interval.distance <- path.interval.dist(tree1, tree2)
    dists[i] = path.interval.distance
  }
  return(dists)
}

#'Read in the Blossoc trees
#'
#'@param treefile the path and file name of the file with the qb trees (type string)
#'@param indexfile the path and file name of the file with the index (type string)
#'@return a list of blossoc trees
readBlossoc <- function(treefile, indexfile) {
qb.trees <- list() # list of QB trees (a list of trees for each run)

#first run
for(i in 1:1) {
  temp.qb.list <- list()
  for(j in 1:length(ms.treepositions[[i]])) {
    qbtree <- read.tree(
      #text=scan(paste0("~/Desktop/Phylogenetics/Updated/NewData/data_Runs600-649/",'qbtrees_run',(i+599)),skip=(j-1)*3+1,
                #nlines=1,what='character',sep='\n')[1])
    text=scan(paste0(treefile,(i+599)),skip=(j-1)*3+1,
              nlines=1,what='character',sep='\n')[1])
    tips <- scan(paste0(indexfile,(i+599)), nlines = (length(qbtree$tip.label)),
              what = 'list', multi.line = FALSE)
    #tips <- scan(paste0("~/Desktop/Phylogenetics/Updated/NewData/data_Runs600-649/indextoorder_run",(i+599)), nlines = (length(qbtree$tip.label)),
                 #what = 'list', multi.line = FALSE)
    tips <- as.numeric(tips)
    for (k in 1:(length(qbtree$tip.label))) { # change tip labels in qtrees to match mstrees
      qbtree$tip.label[k] <- tips[2*(as.numeric(qbtree$tip.label[k]) + 1)] ## convert qbtree's tip lables
    }
    temp.qb.list[[j]] <- qbtree
  }
  qb.trees[[i]] <- temp.qb.list
}

  return(qb.trees)
}


#'Read in the Rent trees
#'
#'@param treefile the path and file name of the file with the rent trees (type string)
#'@return a list of blossoc trees
readRent <- function(treefile) {
  rent.trees <- list() #new list to store all of the rent trees (with regular recombination)

  for(i in 1:1) {
    temp.rent.list <- list()
    for(j in 1:length(ms.treepositions[[i]])) {
      #renttree<-read.tree(text=paste0(scan(paste0("~/Desktop/Phylogenetics/Updated/NewData/renttrees649/",'renttrees_run',(i+599)),
                                           #skip=(j-1),nlines=1,what='character',sep='\n')[1],';'))
      renttree<-read.tree(text=paste0(scan(paste0(treefile,(i+599)),
                                           skip=(j-1),nlines=1,what='character',sep='\n')[1],';'))
      renttree=collapse.singles(renttree)
      temp.rent.list[[j]] <- renttree
    }
    rent.trees[[i]] <- temp.rent.list
  }

  return(rent.trees)
}





