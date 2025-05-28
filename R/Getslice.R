###### Slicing function  ##########

library("treebalance")
library("geiger")
library("paleotree")
library("phytools")
library("ape")
library('treesliceR')


Getslice = function(realtr1, slice, t){
  
  #tree underslice (red)
  realtrundsl=timeSliceTree(realtr1, slice, plot = F, drop.extinct = F)
  realtrundsld=drop.extinct(realtrundsl)
  realtrundsld$root.time=NULL
  
  #tree overslice (blue)
  realtrovsl=timeSliceTree(realtr1, slice-t, plot = F, drop.extinct = F)
  realtrovsl$root.time=NULL
  
  #subtrees from overslice
  subtrees=treeSlice(realtrovsl, realtr1$root.time-slice, trivial = T)
  realtrovsld=drop.extinct(realtrovsl)
  
  ########################
  #this part is new
  #order subtrees same way as the labels in underslice by checking if the underslice label matches some of the subtree labels
  #one could do it smarter with some grep function I guess but this works
  subtreesord=list()
  for(j in 1:length(realtrundsld$tip.label)){
    isinsubtree=rep(F, length(subtrees))
    for (k in 1:length(subtrees)){
      isinsubtree[k]=realtrundsld$tip.label[j] %in% subtrees[[k]]$tip.label
    }
    subtreesord[[j]]=subtrees[[which(isinsubtree)]]
  }
  
  ##############################
  #from here on, I just changed subtrees->subtreesord
  
  #calculate tips
  tipcounts <- c()
  for (i in 1:length(subtreesord)){
    c = 0
    for ( j in 1:length(subtreesord[[i]]$tip.label)){
      ### count the number of subtree if tips are alive
      if (subtreesord[[i]]$tip.label[j] %in% realtrovsld$tip.label){
        c = c+1
      }
    }
    tipcounts <- append(tipcounts, c) 
  }
  
  
  return(list(realtrundsld, tipcounts, realtrovsl))
  
}
 





