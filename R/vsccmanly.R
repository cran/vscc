vsccmanly <-function(x, G=2:9, numstart=100,
                     selection="backward",forcereduction=FALSE,
                     initstart="k-means", seedval=2354){
  
  if(initstart!="k-means" & initstart!="hierarchical") {
    stop("No valid initialization method has been provided.")
  }
  if(selection!="forward" & selection!="backward" & selection!="none"){
    stop("No valid selection method has been provided.")
  }
  
  origx <- x 
  origG <- G
  x <- x
  p <- ncol(x)
  initruns=list()
  Gmod=list()
  initbic=vector(length=length(G))
  Mmod=list()
  initbic2=vector(length=length(G))
  listfcn=function(dim1){
    nameslist=c("la","tau","Mu","S","gamma","id","ll","bic","iter","flag" )
    list1 <- as.list(rep(NA, length(nameslist)))
    names(list1) <- nameslist
    return(list1)
  }
  listfcn2=function(dim1){
    nameslist=c("la","Mu","S","id","iter","flag" )
    list1 <- as.list(rep(NA, length(nameslist)))
    names(list1) <- nameslist
    return(list1)
  }
  
  if(selection=="backward"){
    
    for(j in 1:length(G)){
      Gtest=G[j]
      
      if(initstart=="k-means"){
        set.seed(seedval)
        id.start <- kmeans(x, Gtest, nstart=numstart)$cluster
      }
      else if(initstart=="hierarchical"){
        H <- hclust(dist(x), method = "ward.D")
        id.start<-cutree(H, Gtest)
      }
      
      #Following initialization steps aid in convergence
      la <- matrix(0.1, Gtest, p)
      C <- tryCatch({Manly.Kmeans(x, id = id.start, la = la)
      }, error=function(e) listfcn2(6))
      id.CK<-C$id
      
      suppressWarnings({
      Gmod[[j]] <- tryCatch({Manly.EM(x,id=id.CK,la=la)
      },error=function(e) listfcn(10))
      })
      initbic[j]<-Gmod[[j]]$bic
      
    }
    G=origG[which.min(initbic)]
    Gfinmod=Gmod[[which.min(initbic)]]
    initrun <-tryCatch({ Manly.select(x, model = Gfinmod, method = "backward") #transformation parameter selection
    },error=function(e) listfcn(10))
    x$guessid=initrun$id
  }
  else if (selection=="forward"){
    for(j in 1:length(G)){
      Gtest=G[j]
      if(initstart=="k-means"){
        set.seed(seedval)
        id.start <- kmeans(x, Gtest, nstart=numstart)$cluster
      }
      else if(initstart=="hierarchical"){
        H <- hclust(dist(x), method = "ward.D")
        id.start<-cutree(H, Gtest)
      }
      
      #Following initialization steps aid in convergence
      id.CK=id.start
      suppressWarnings({
      Gmod[[j]] <- tryCatch({Manly.EM(x,id=id.CK,la=matrix(0, Gtest, p)) 
      }, error=function(e) listfcn(10))
      })
      initbic[j]<-Gmod[[j]]$bic
    }
    G=origG[which.min(initbic)]
    Gfinmod=Gmod[[which.min(initbic)]]
    
    initrun <- tryCatch({
      Manly.select(x, model = Gfinmod, method = "forward")
    }, error=function(e) listfcn(10))
    x$guessid=initrun$id
  }
  
  else if (selection=="none"){ #full Manly
    for(j in 1:length(G)){
      Gtest=G[j]
      
      if(initstart=="k-means"){
        set.seed(seedval)
        id.start <- kmeans(x, Gtest, nstart=numstart)$cluster
      }
      else if(initstart=="hierarchical"){
        H <- hclust(dist(x), method = "ward.D")
        id.start<-cutree(H, Gtest)
      }
      #Following initialization steps aid in convergence
      la <- matrix(0.1, Gtest, p)
      C <- tryCatch({Manly.Kmeans(x, id = id.start, la = la)
      }, error=function(e) listfcn2(6))
      id.CK<-C$id
      
      Gmod[[j]] <- tryCatch({Manly.EM(x,id=id.CK,la=la) 
      },error=function(e) listfcn(10))
      initbic[j]<-Gmod[[j]]$bic 
      
    }
    G=origG[which.min(initbic)]
    initrun=Gmod[[which.min(initbic)]] #Gmod is used as there is no transformation parameter selection for the full model
    x$guessid=initrun$id
  }
  
  #Transforming data
  Xnew=as.vector(NA)
  transdata=matrix(NA,nrow=nrow(x),ncol=p)
  x$guessid=initrun$id
  for(j in 1:p){
    Xnew=as.vector(NA)
    for(i in 1:G){
      groupD=x[x$guessid==i,j]
      varlambda=initrun$la[i,j]
      if(varlambda==0){
        Xtrans=groupD
      }
      else{
        Xtrans=(exp(varlambda*groupD)-1)/varlambda
      }
      Xnew=c(Xnew,Xtrans)
    }
    Xnew=Xnew[-1]
    transdata[,j]<-Xnew
  }
  transdata<-as.data.frame(transdata)
  transdata<-scale(transdata)
  colnames(transdata)<-colnames(origx)
  #End of transformation
  
  initial <- initrun$id
  initunc=sum(1-apply(initrun$gamma,1,max)) #uncertainty
  
  G <- length(unique(initial))
  n <- nrow(x)
  zmat <- matrix(0,n,G)
  for(i in 1:G){
    for(j in 1:n){
      if(initial[j]==i){
        zmat[j, i]<-1 
      }
    }
  }
  
  ng <- colSums(zmat)
  meang <- matrix(0,G,p) #Mean for group and variable
  for(g in 1:G){
    meang[g,] <- colSums(zmat[,g]*transdata)/ng[g]
  }
  
  wvar=matrix(0,1,p) #Within-group variance for each variable.
  for(j in 1:p){
    for(g in 1:G){
      for(i in 1:n){
        wvar[,j] = wvar[,j]+zmat[i,g]*(transdata[i,j] - meang[g,j])^2
      }
    }
    wvar[,j]=wvar[,j]/n
  }
  colnames(wvar)<-colnames(transdata)
  
  #Lines 189-220 come from vscc function in the vscc package (Andrews and McNicholas, 2013).
  #The only difference is that selection occurs on the transformed variables and the original
  #variables are stored in 'select' as the Manly models are fit to the original data.
  sorted <- t(as.matrix(sort(wvar[1,])))
  select <- list()
  transelect <- list()
  varnames <- list()
  modrun <- list()
  numvars <- NA
  for(i in 1:5){ #Automatically select the variable that minimizes wmat into each set
    select[[i]] <-  matrix(data=origx[,colnames(sorted)[1]]) 
    transelect[[i]] <- matrix(data=transdata[,colnames(sorted)[1]])
    varnames[[i]] <- colnames(sorted)[1]
  }
  counts <- rep(2,5)
  for(k in 2:p){
    curname <- colnames(sorted)[k]
    for(i in 1:5){
      curcor <- cor(cbind(transdata[,curname],transelect[[i]]))
      if(all(abs(curcor[upper.tri(curcor)])<=(1-sorted[1,k]^i))){ #Selection criteria
        select[[i]] <- cbind(select[[i]],origx[,curname])
        transelect[[i]] <- cbind(transelect[[i]],transdata[,curname])
        varnames[[i]][counts[i]] <- curname
        counts[i] <- counts[i]+1
      }
    }
  }
  for(i in 1:5){
    colnames(select[[i]]) <- varnames[[i]]
  }
  
  moduncs <- Inf
  runif <- rep(TRUE,5)
  
  for(i in 1:4){
    for(j in (i+1):5){
      if(runif[j]==TRUE && identical(varnames[[i]],varnames[[j]])){
        runif[j]=FALSE # If subset contains same variables as another subset assign runif FALSE
      }
    }
  }
  
  for(i in 1:5){
    if(runif[i]){
      G=origG
      
      if(selection=="backward"){
        
        for(j in 1:length(G)){
          Gtest=G[j]
          
          if(initstart=="k-means"){
            set.seed(seedval)
            id.start <- kmeans(select[[i]], Gtest,nstart=numstart)$cluster
          }
          else if(initstart=="hierarchical"){
            H <- hclust(dist(select[[i]]), method = "ward.D")
            id.start<-cutree(H, Gtest)
          }
          #Following initialization steps aid in convergence
          p=ncol(select[[i]])
          la <- matrix(0.1, Gtest, p)
          C <- tryCatch({Manly.Kmeans(select[[i]], id = id.start, la = la)
          }, error=function(e) listfcn2(6))
          id.CK<-C$id
          
          suppressWarnings({
          Mmod[[j]] <- tryCatch({Manly.EM(select[[i]], id = id.CK,la=la)
          },error=function(e) listfcn(10))
          })
          initbic2[j]<-Mmod[[j]]$bic
          
        }
        Gnew=origG[which.min(initbic2)]
        Mfinmod=Mmod[[which.min(initbic2)]]
        
        modrun[[i]] <- tryCatch({Manly.select(select[[i]], model = Mfinmod, method = "backward")
        }, error=function(e) listfcn(10))
      }
      else if (selection=="forward"){
        for(j in 1:length(G)){
          Gtest=G[j]
          if(initstart=="k-means"){
            set.seed(seedval)
            id.start <- kmeans(select[[i]], Gtest,nstart=numstart)$cluster
          }
          else if(initstart=="hierarchical"){
            H <- hclust(dist(select[[i]]), method = "ward.D")
            id.start<-cutree(H, Gtest)
          }
          
          p=ncol(select[[i]])
          
          id.CK=id.start
          #Following initialization steps aid in convergence
          suppressWarnings({
          Mmod[[j]] <- tryCatch({Manly.EM(select[[i]], id = id.CK,la = matrix(0, Gtest, p))
          }, error=function(e) listfcn(10))
          })
          id.G <- Gmod$id
          initbic2[j]<-Mmod[[j]]$bic
        }
        Gnew=origG[which.min(initbic2)]
        Mfinmod=Mmod[[which.min(initbic2)]]
        
        modrun[[i]] <- tryCatch({
          inittest <- Manly.select(select[[i]], model = Mfinmod, method = "forward")
        }, error=function(e) listfcn(10))
        
        
      }
      else if (selection=="none"){
        for(j in 1:length(G)){
          Gtest=G[j]
          if(initstart=="k-means"){
            set.seed(seedval)
            id.start <- kmeans(select[[i]], Gtest,nstart=numstart)$cluster
          }
          else if(initstart=="hierarchical"){
            H <- hclust(dist(select[[i]]), method = "ward.D")
            id.start<-cutree(H, Gtest)
          }
          
          #Following initialization steps aid in convergence
          p=ncol(select[[i]])
          la <- matrix(0.1, Gtest, p)
          C <- tryCatch({Manly.Kmeans(select[[i]], id = id.start, la = la)
          }, error=function(e) listfcn2(6))
          id.CK<-C$id
          
          
          Mmod[[j]] <- tryCatch({Manly.EM(select[[i]], id = id.CK,la = la)
          }, error=function(e) listfcn(10))
          id.G <- Gmod$id
          initbic2[j]<-Mmod[[j]]$bic
        }
        Gnew=origG[which.min(initbic2)]
        modrun[[i]]=Mmod[[which.min(initbic2)]]
        
        
      }
      
      moduncs[i]=tryCatch({sum(1-apply(modrun[[i]]$gamma,1,max))
      }, error=function(e) NA)
      
    }
    else{
      modrun[[i]] <- "Same as simpler relation"
      moduncs[i] <- Inf
    }
  }
  
  #change bic to bigger is better
  suppressWarnings({
    for(i in 1:5){
      modrun[[i]]$bic=tryCatch({-1*modrun[[i]]$bic}, error=function(e) NA)
    }
  })
  
  #The following code, used to output results, comes from the vscc function in the vscc package (Andrews and McNicholas, 2013) 
  #to ensure consistency between functions.
  store <- list()
  store[["selected"]] <- select
  
  store[["initialrun"]] <- initrun
  
  if(forcereduction){
    store[["bestmodel"]] <- modrun[[which.min(moduncs)]]
    store[["chosenrelation"]] <- which.min(moduncs)
    store[["variables"]]<-names(as.data.frame(select[[which.min(moduncs)]]))
    store[["topselected"]] <- select[[which.min(moduncs)]]
    store[["uncertainty"]] <- min(moduncs)
  }
  else{
    if(min(moduncs)<initunc){
      store[["bestmodel"]] <- modrun[[which.min(moduncs)]]
      store[["chosenrelation"]] <- which.min(moduncs)
      store[["topselected"]] <- select[[which.min(moduncs)]]
      store[["uncertainty"]] <- min(moduncs)
      store[["variables"]]<-names(as.data.frame(select[[which.min(moduncs)]]))
      
    }
    else{
      store[["bestmodel"]] <- initrun
      store[["chosenrelation"]] <- "Full dataset"
      store[["topselected"]] <- origx
      store[["uncertainty"]] <- initunc
    }
  }
  
  store[["allmodelfit"]] <- modrun
  
  store[["wss"]] <- sorted
  class(store) <- "vsccmanly"
  store
}
