AMOEBA<-function(outc,neig,power=1,cpu=1){
  ### Checking
  # neig and outc
  if (!inherits(neig, "nb")) 
    stop(gettextf('%s must be a %s',sQuote('neig'),dQuote('list')))
  if (!is.vector(outc) | !is.numeric(outc))
    stop(gettextf('%s must be a %s',sQuote('outc'),dQuote('numeric vector')))
  if (length(neig)!=length(outc))
    stop(gettextf('%s and %s must have the same lenght',sQuote('neig'),dQuote('outc')))
  # power and cpu
  if (cpu<=0 | !is.numeric(cpu))
    stop(gettextf('%s must be a %s',sQuote('cpu'),dQuote('positive integer')))
  if (power<=0 | !is.numeric(power))
    stop(gettextf('%s must be a %s',sQuote('power'),dQuote('positive integer')))  
  is.wholenumber<-function(x,tol=.Machine$double.eps^0.5)
    abs(x - round(x)) < tol
  if (!is.wholenumber(cpu))
    stop(gettextf('%s must be a %s',sQuote('cpu'),dQuote('positive integer')))
  if (!is.wholenumber(power))
    stop(gettextf('%s must be a %s',sQuote('power'),dQuote('positive integer')))  
  super<-0
  while (length(outc) >= 3^(super+1))
    super<-super+1
  if (power>super)
    stop(gettextf(paste('%s must be equal or less than',super,', and a %s'),sQuote('power'),dQuote('positive integer')))
  
  ### Auxilary functions
  Classi<-function(x,y,cpu){
    ### Auxilary functions
    # Local Getis-Ord Statistic
    Getis.Ord<-function(num){
      # 1 neighbor
      if (length(num)==1){
        if (is.na(num))
          return(NA)
        else if ((num-total.mean)==0)
          return(0)
        else
          return((num-total.mean)/total.sd)
      }
      # 2 or more neighbors
      else{
        if ((sum(num)-(length(num)*total.mean))==0)
          return(0)
        else
          return((sum(num)-(length(num)*total.mean))/(total.sd*sqrt(((total*length(num))-(length(num)^2))/(total-1))))
      }
    }
    
    # Looking ecotope
    Clustering<-function(seed){
      ecotope<-seed             # Vector with the neighbors of the cluster
      out<-0                    # Vector with the eliminated neighbors
      stat0<-Getis.Ord(y[seed]) # Initial Getis-Ord local statistic
      if (is.na(stat0)){
        return(list(0,NA,ecotope))
      }
      else{
        stat1<-stat0
        # While we can improve
        while (abs(stat0)<=abs(stat1)){
          a.ecotope<-unique(unlist(x[ecotope]))                 # New neighbors
          a.ecotope<-a.ecotope[!is.element(a.ecotope,ecotope)]  # Eliminate those already
          a.ecotope<-a.ecotope[!is.element(a.ecotope,out)]      # Eliminate those that are no longer
          b.ecotope<-y[a.ecotope]                               # New values of outcome
          out<-c(out,a.ecotope[is.na(b.ecotope)])               # Removing missing values
          a.ecotope<-a.ecotope[!is.na(b.ecotope)]
          b.ecotope<-b.ecotope[!is.na(b.ecotope)]
          if (length(a.ecotope)>0 & length(b.ecotope)>0){
            # Order previous values
            if (stat0>0)
              c.ecotope<-order(b.ecotope,decreasing=T)  
            else
              c.ecotope<-order(b.ecotope)  
            d.ecotope<-a.ecotope[c.ecotope]                     # Orderly neighborhood
            while (length(d.ecotope)>0){
              stat2<-Getis.Ord(y[c(ecotope,d.ecotope[1])])      # New Getis-Ord local statistic
              # Checking if we improve
              if (abs(stat1)<abs(stat2) & ((stat1<0 & stat2<0) | (stat1>0 & stat2>0))){
                ecotope<-c(ecotope,d.ecotope[1])
                d.ecotope<-d.ecotope[-1]
                stat0<-stat1<-stat2
              }
              else{
                out<-c(out,d.ecotope)
                d.ecotope<-vector()
              }
            }
          }
          # No more improved
          else{
            if (stat1==0)
              sal<-list(2,stat1,ecotope)
            else if(stat1>0)
              sal<-list(3,stat1,ecotope)
            else
              sal<-list(1,stat1,ecotope)
            stat0<-1
            stat1<-0
          }
        }
      }
      return(sal)
    }
    
    ### Application
    # Total values
    total<-length(!is.na(y))
    total.mean<-mean(y,na.rm=T)
    total.sd<-sd(y,na.rm=T)
    
    # Ecotope for every neighbor
    sfInit(parallel=T,cpus=cpu,type='SOCK',socketHosts=rep('localhost',cpu))
    sfExport('x','y','Getis.Ord','total','total.mean','total.sd')
    sfClusterSetupRNG()
    sal<-sfLapply(1:length(y),Clustering)
    sfStop()
    
    # Overlapping
    res<-sapply(sal,'[[',1)
    index<-sapply(sal,'[[',2)
    index.list<-sapply(sal,'[[',3)
    rm(sal,total,total.mean,total.sd,Getis.Ord,Clustering)
    
    # High values
    while(any(res==3)){
      pos<-unlist(index.list[res==3 & max(index,na.rm=T)==index][1])
      res[pos]<-6
      index[pos]<-0
      pos2<-unlist(lapply(index.list[res==3], function(x) if (any(x%in%pos)) x[1]))
      res[pos2]<-5
      index[pos2]<-0
    }
    
    # Low values
    while(any(res==1)){
      pos<-unlist(index.list[res==1 & min(index,na.rm=T)==index][1])
      res[pos]<-4
      index[pos]<-0
      pos2<-unlist(lapply(index.list[res==1], function(x) if (any(x%in%pos)) x[1]))
      res[pos2]<-5
      index[pos2]<-0
    }
    
    # Middle values
    res<-ifelse(res==2,5,res)
    
    # Not values
    res<-ifelse(res==0,NA,res)
    
    # Return (4=Low,5=Middle,6=High)
    return(res)
  }
  
  ### First application
  sal<-Classi(neig,outc,cpu)
  
  ### Plus application
  if (power>1){
    cont<-2
    while (cont<=power){
      cate<-1
      movement<-unique(sal)[order(unique(sal))]
      movement<-movement[!is.na(movement)]
      sal2<-rep(NA,length(sal))
      while (length(movement)>0){
        aux1<-sort(sal==movement[1] & !is.na(sal),index.return=T,decreasing=T)
        aux2<-aux1$ix[1:sum(sal==movement[1],na.rm=T)]
        if (sum(aux1$x)>1){
          neig.bis<-subset(neig,(1:length(neig)%in%aux2))
          outc.bis<-outc[aux2]
          sal3<-Classi(neig.bis,outc.bis,cpu)
          sal2[aux2]<-ifelse(sal3==4,cate,ifelse(sal3==5,cate+1,cate+2))
        }
        else
          sal2[aux2]<-cate
        movement<-movement[-1]
        cate<-cate+3
      }
      sal<-sal2
      cont<-cont+1
    }
  }
  
  ### Export
  sal.aux1<-unique(as.numeric(names(table(sal))))
  sal.aux2<-sort(sal.aux1,index.return=T)
  sal.aux3<-apply(as.array(sal),1,function(x) if(is.na(x)) NA else sal.aux2$ix[sal.aux2$x==x])
  return(sal.aux3)
}