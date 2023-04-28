#setwd("C:/Users/Hu/Google Drive/HEA_bcc/data/ML/prediction")

#sys="Ti-V-Ta-Mo-W"

des_surf <- function (comp=comp){
#  file.1=paste("./comp_grid/comp_grid_",sys,".csv",sep = '')
#  file.2=paste("./descriptors/descriptors_surf_",sys,".csv",sep = '')
#  data.1=read.csv(file=file.1,header=T)
  comp_n=comp/sum(comp)
  a=c(1,1,1,1,1,1,1,1,1,1)
  data.1=data.frame(rbind(a,comp_n))
  
  file.3="elemental_features.csv"
  data.3=read.csv(file=file.3,header=T)
  para_ele=as.matrix(data.3[,2:4])
  
  
  load(file="features_surf.Rdata")
  
  feature.surf.select=feature.surf.all
  n.ob=length(data.1[,1])
  n.feature=length(feature.surf.select[1,1,])
  n.feature.ele=length(para_ele[1,])
  n.feature.all=n.feature+n.feature.ele
  
  parameter.surf=array(0,c(n.ob,10,n.feature))
  
  comp_ele=as.matrix(data.1)
  for (i in seq(n.feature)){
    parameter.surf[,,i]=comp_ele%*%as.matrix(feature.surf.select[,,i])
  }
  
  order=c(1)
  size=2*n.feature.all
  #+length(order)*n.feature*length(order)*(length(para_ele[1,]))
  #size=2*n.feature.all+2*2*n.feature*n.feature.ele
  predictors=array(0,c(n.ob,size))
  for (i in seq(n.ob)){
    #cat(i,'\r')
    #avg1=NULL
    avg=NULL
    nonzvec = which(comp_ele[i,]!=0)
    avg1=NULL
    for (k in seq(n.feature)){
      para.temp=parameter.surf[i,nonzvec,k]
      #avg1=NULL
      avgtemp=para.temp%*%comp_ele[i,nonzvec]
      #dvetemp=((para.temp-avgtemp[1])^2%*%comp_ele[i,nonzvec])^(1/2)
      dvetemp=((1/(1-sum(comp_ele[i,nonzvec]^2)))*((para.temp-avgtemp[1])^2%*%comp_ele[i,nonzvec]))^(1/2)
      avg1=c(avg1,avgtemp,dvetemp)
      #cat(para.temp,'\n')
      #cat(avgtemp,'\n')
      #cat(dvetemp,'\n')
      
    }
    #avg=c(avg,avg1)
    
    
    avg2=NULL
    for (n in seq(n.feature.ele)){
      avgtemp=comp_ele[i,nonzvec]%*%para_ele[nonzvec,n]
      #dvetemp=(((para_ele[nonzvec,n]-avgtemp[1])^2%*%comp_ele[i,nonzvec]))^(1/2)
      dvetemp=((1/(1-sum(comp_ele[i,nonzvec]^2)))*((para_ele[nonzvec,n]-avgtemp[1])^2%*%comp_ele[i,nonzvec]))^(1/2)
      avg2=c(avg2,avgtemp,dvetemp)
    }
    
    #  cat (avg2, '\n')
    #x=length(avg)
    #y=length(avg2)
    #interact=NULL
    #for (x in seq(length(avg2))){
    #  for (y in seq(length(avg1))){
    #    temp3=avg2[x]*avg1[y]
    #    interact=c(interact,temp3)
    #  }
    #}
    avg=c(avg1,avg2)
    
    predictors[i,]=avg
  }
  
  
  all=data.frame(predictors)
  all[is.na(all)]=0
  return(all)
}




