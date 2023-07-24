##AFS package+ Analyse factorielle g�n�rale

library("ggplot2")
library("ggpmisc")
library("ggrepel")

dist_geod<- function(x,y){
  return (2*acos(sum(sqrt(x*y))))
}

Hellingerizer<-function(x){

  return (sign(x)*sqrt(abs(x)))
}
turn_to_freq<-function(x){
  return(x/sum(abs(x)))
}


setClass("matrice_d",representation(mat="matrix",icol="numeric",ili="numeric"))

depistage<-function(f){
  J<-length(f[1,])
  I<-length(f[,1])
  F2<-f
  ind_li<-c()
  ind_col<-c()
  for(j in 1:J){
    if(sum(f[,j])==0){F2<-f[,-j]
    
    ind_col<-c(ind_col,j)}}
  for (i in 1:I){
    if(sum(f[i,])==0){F2<-f[-i,]
    ind_li<-c(ind_li,j)}}
  if(is.null(ind_col)){ind_col<-c(0)}
  if(is.null(ind_li)){ind_li<-c(0)}
  r=new("matrice_d",mat=F2,ili=ind_li,icol=ind_col)
  return(r)
}


Nuage_X<-function(t,mI=poids_lignes(t)){
  Z<-depistage(t)
  X<-Z@mat
  if (Z@ili != c(0)){
    linull<-Z@ili
    n<-length(linull)
    for (x in 1:n){ i<-linull[x]
    mI<-mI[-i]}
  }
  X<-diag(1/sqrt(abs(mI)),ncol=length(mI),nrow=length(mI))%*%X
  return(X)
}

Nuage_Y<-function(t,mJ=poids_colonne(t)){
  Z<-depistage(t)
  Y<-Z@mat
  if (Z@icol != c(0)){
    conull<-Z@icol
    n<-length(conull)
    for (x in 1:n){j<-conull[x]
    mJ<-mJ[-j]}}
  Y<-Y%*%diag(1/sqrt(abs(mJ)),ncol=length(mJ),nrow=length(mJ))

  return(Y)
}
poids_colonne<-function(f){
  J<-length(f[1,])
  mJ<-rep(0,J)
  for (j in 1:J){
    mJ[j]<-sum(f[,j])
                    }
  return(mJ)
}
poids_lignes<-function(f){
  J<-length(f[,1])
  mI<-rep(0,J)
  for (j in 1:J){
    mI[j]<-sum(f[j,])}
  return(mI)}


        
         
I_alpha<-function (p,q,alpha){
  return ((1/(alpha-1))*log2((sum(((p^alpha)/(q^alpha))*p))))

}

dist_Hellinger<-function(p,q){
  return (sum((sqrt(p)-sqrt(q))^2))
}

setClass("AF",representation(t="matrix",lbd="numeric",U_vects="matrix",W_vects="matrix",compos_F="matrix",compos_G="matrix",
                                 cos2="matrix",cos2_du="matrix",CTR="matrix",CTR_du="matrix",INR="matrix",INR_du="matrix",selection_des_axes="data.frame"))


produit_des_marges<-function(f){
  
  phi<-matrix(poids_lignes(f),ncol=1)%*%aperm(matrix(poids_colonne(f),ncol = 1))

return (phi)
}


AF_corres<-function(t){
  t<-turn_to_freq(depistage(t)@mat)
  phi<-turn_to_freq(depistage(produit_des_marges(t))@mat)
  Tableau<-(t-phi)/sqrt(phi)
  mI=poids_lignes(t)
  mJ=poids_colonne(t)
  A<-AF(Tableau,mI,mJ)
  return (A)
}


AF_plot<-function(B,pcx,pcy,choix="both",RN=row.names(B@t),CN=colnames(B@t),Title="Analyse factorielle"){
  if (pcx>5 ||	pcy>5){return("L'analyse factorielle ne s'effectue pas au dela du 5�me axe ")}
  else {
    if (is.null(RN)){ RN<-as.character(1:length(B@t[,1])) }
    
    if (is.null(CN)){ CN<-as.character(1:length(B@t[1,])) }
       
  
  comp_F<-as.data.frame(cbind(B@compos_F))
  comp_G<-as.data.frame(cbind(B@compos_G))
  charx<-paste("V",pcx,sep="")
  chary<-paste("V",pcy,sep="")
  Ix<-(trunc((B@lbd[pcx])*(1000/sum(B@t^2)),digits=3))/10
  Iy<-(trunc((B@lbd[pcy])*(1000/sum(B@t^2)),digits=3))/10
  chx<-paste(Ix,"%",sep="")
  chy<-paste(Iy,"%",sep="")
  chix<-paste(chx,"de l'inertie",sep=" ")
  chiy<-paste(chy,"de l'inertie",sep=" ")
  lax<-paste("Axe",pcx,sep=" ")
  lay<-paste("Axe",pcy,sep=" ")
  labx<-paste(lax,chix,sep="  ")
  laby<-paste(lay,chiy,sep="  ")
 
   if (choix=="ind"){
  
  ggplot() +
    geom_label_repel(data=comp_F,hjust=1,vjust=1,aes_string(x=charx, y=chary),colour="white",size=2.3,fill="darkgreen" ,label=as.vector(RN) )+
      geom_point(data=comp_F,aes_string(x=charx, y=chary),colour="darkgreen" )+
       ggtitle(Title)+xlab(labx)+ylab(laby)+ 
    geom_quadrant_lines(linetype = "solid")+
      scale_x_continuous(limits = symmetric_limits(c(B@compos_F[,pcx]))) +
      scale_y_continuous(limits = symmetric_limits(c(B@compos_F[,pcy])))+
     xlab(labx)+ylab(laby)+
    theme_minimal()
  
  }
  else {if(choix=="var"){
    
    ggplot() +
      geom_label_repel(data=comp_G,hjust=1,vjust=1,aes_string(x=charx,y=chary),colour="white",fill="darkred",size=2.3,label=as.vector(CN))+
      geom_point(data=comp_G,aes_string(x=charx, y=chary),colour="darkred" )+
      ggtitle(Title)+xlab(labx)+ylab(laby)+ 
      geom_quadrant_lines(linetype = "solid")+
      scale_x_continuous(limits = symmetric_limits(c(B@compos_G[,pcx]))) +
      scale_y_continuous(limits = symmetric_limits(c(B@compos_G[,pcy])))+
       
      xlab(labx)+ylab(laby)+
      theme_minimal()
      
  }
    else {
      
      ggplot() +
        geom_label_repel(data=comp_F,hjust=0,vjust=0,aes_string(x=charx, y=chary),colour="white",size=2.3 ,label=as.vector(RN),fill="darkgreen")+
        geom_label_repel(data=comp_G,hjust=0,vjust=0,aes_string(x=charx,y=chary),colour="white",size=2.3,label=as.vector(CN),fill="darkred")+
        geom_point(data=comp_F,aes_string(x=charx, y=chary),colour="darkgreen" )+
      geom_point(data=comp_G,aes_string(x=charx, y=chary),colour="darkred" )+
        ggtitle(Title)+xlab(labx)+ylab(laby)+
        geom_quadrant_lines(linetype = "solid")+
        scale_x_continuous(limits = symmetric_limits(c(B@compos_F[,pcx],B@compos_G[,pcx]))) +
        scale_y_continuous(limits = symmetric_limits(c(B@compos_F[,pcy],B@compos_G[,pcy])))+
        theme_minimal()
    }}
  }}





AF<-function(t,mI=rep(1,length(depistage(t)@mat[,1])),mJ=rep(1,length(depistage(t)@mat[1,])),nbvrp=5){
  

  ##Choix de l'analyse 
  
  X<-Nuage_X(t,mI)
  Y<-Nuage_Y(t,mJ)
  I<-length(t[,1])
  J<-length(t[1,])
  ##diagonalisation de l'inertie du nuage de plus petite dimension:
  
  V<-aperm(t)%*%t
  lambda<-eigen(V)$values
  Dlambda<-diag(1/sqrt(abs(lambda)))
  U<-eigen(V)$vectors 
  
  ##w vecteur propre de TT' associé à lambda 
  w<-(t%*%U)%*%Dlambda
  composante_G<-matrix(0,J,nbvrp)
  composante_F<-matrix(0,I,nbvrp)
  
  for (k in 1:nbvrp){
    for (i in 1:I){
      composante_F[i,k]<-X[i,]%*%U[,k]
    }
    for (j in 1:J){
      composante_G[j,k]<-Y[,j]%*%w[,k]
    }
  }
  INR<-matrix(0,I,1)
  row.names(INR)
  INR_du<-matrix(0,J,1)
  
  CTR<- matrix(0,I,nbvrp)
  CTR_du<- matrix(0,J,nbvrp)
  cos_2<-CTR
  cos_2_du<-CTR_du
  ##INDICATEURS: 
  for  (i in 1:I){
    INR[i]<-sum(t[i,]^2)/sum(t^2) 
    for (k in 1:nbvrp){
      CTR[i,k]<-(mI[i]*(composante_F[i,k])^2)/lambda[k]  ##INRibution d'un individu i à l'inertie de l'axe k
      cos_2[i,k]<-(lambda[k]*w[i,k]^2)/sum(t[i,]^2) ##part visible de l'intertie de l'individu i sur l'axe k
    }}
  for  (j in 1:J){
    INR_du[j]<-sum(t[,j]^2)/sum(t^2) ##INRibution d'un axe à l'inertie totale
    for (k in 1:nbvrp){
      CTR_du[j,k]<-U[j,k]^2    ##INRibution de la variable j à l'inertie de l'axe k
      cos_2_du[j,k]<-(lambda[k]*U[j,k]^2)/sum(t[,j]^2) ##part visible de l'interie de la variable j sur l'axe k
    }}
  row.names(cos_2)<-row.names(t)
  row.names(cos_2_du)<-colnames(t)
  
  row.names(CTR)<-row.names(t)
  row.names(CTR_du)<-colnames(t)
  colnames(INR)<-"Contribution des lignes a l'inertie"
  colnames(INR_du)<-"Contribution des colonnes a l'inertie"
  row.names(INR)<-row.names(t)
  row.names(INR_du)<-colnames(t)
  ## minimum des deux dimensons en valeur par défaut dans la sélection des axes
  
  
  selection_des_axes<-data.frame("Rang"=(1:nbvrp),"Proportion d'inertie"=lambda[1:nbvrp]/sum(t^2),"Proportion d'inertie cumulee"=cumsum(lambda[1:nbvrp]/sum(t^2)))
  

  
  
  
  
  AF=new("AF",t=t,lbd=lambda[1:nbvrp],U_vects=U,W_vects=w,compos_F=composante_F,compos_G=composante_G,cos2=cos_2,cos2_du=cos_2_du,CTR=CTR,CTR_du=CTR_du,INR=INR,INR_du=INR_du,selection_des_axes=selection_des_axes)
  
  
  
  return(AF)
  
}  

##fonction qui trace un cercle
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(xc = xx, yc = yy))
}

AFS_corres<-function(t){
  t<-depistage(t)@mat 
  t<-turn_to_freq(t)
  
  mI<-poids_lignes(t)
  mJ<-poids_colonne(t)
  phi<-produit_des_marges(t)
  
  Tableau_final<-Hellingerizer(t)-Hellingerizer(phi)
  row.names(Tableau_final)<-row.names(t)
  colnames(Tableau_final)<-colnames(t)
  
  A<-AF(Tableau_final,mI,mJ)

return(A)
  
}


AFS_comp<-function(t,phi,mI=rep(1,length(t[,1])),mJ=rep(1,length(t[1,]))){

  t<-turn_to_freq(t)
  phi<-turn_to_freq(phi)
  Tableau_final<-Hellingerizer(t)-Hellingerizer(phi)
  Tableau_final<-depistage(Tableau_final)@mat
  row.names(Tableau_final)<-row.names(t)
  colnames(Tableau_final)<-colnames(t)
  
  A<-AF(Tableau_final,mI,mJ)
  
  return(A)
  
}
AFS_echanges<-function(t){
  t<-turn_to_freq(depistage(t)@mat)
  Tableau_final<-t-diag(diag(t),nrow=length(t[1,]),ncol=length(t[1,]))
  Mi<-poids_lignes(Tableau_final)
  Mj<-poids_colonne(Tableau_final)
  row.names(Tableau_final)<-row.names(t)
  colnames(Tableau_final)<-colnames(t)
  
  Res<-AF(Tableau_final,Mi,Mj)
  
  return(Res)
  
}

AFS_transition<-function(t){
  
  t1<-depistage(t)@mat
  t2<-turn_to_freq(diag(diag(t),nrow=length(t[1,]),ncol=length(t[1,])))
  row.names(t1)<-row.names(t)
  Res<-AFS_comp(t=t1,phi=t2)
  return(Res)
  
  
}


QLT<-function(B,vec,choix="ind",RN=row.names(B@t),CN=colnames(B@t))
 
   {if (choix=="ind"){
  QLT<-matrix(0,length(B@t[,1]),ncol = 1)
  for (i in 1:length(B@t[,1])){
    QLT[i,]<-sum(B@cos2[i,vec])
    row.names(QLT)<-RN
    char1<-paste("Qualite de representation des lignes dans le plan(",vec[1],sep="")
    char2<-paste(char1,vec[2],sep=",")
    char3<-paste(char2,")",sep="")
    colnames(QLT)<-char3
    }
    }

  else{QLT<-matrix(0,length(B@t[1,]),ncol=1)
  for (j in 1:length(B@t[1,])){
      QLT[j,]<-sum(B@cos2_du[j,vec])
      row.names(QLT)<-CN
      char1<-paste("Qualite de representation des colonnes dans le plan(",vec[1],sep="")
      char2<-paste(char1,vec[2],sep=",")
      char3<-paste(char2,")",sep="")
      colnames(QLT)<-char3
      
  }
  }
  return(QLT)
}

missSVD<-function(t){
  NA_places<-which(is.na(t))
  
  
  for (i in 1:200){
    
SVD_obj<- AF(t,nbvrp =min(length(t[,1]),length(t[1,])))
  df<-SVD_obj@U_vects%*%diag(SVD_obj@lbd)%*%SVD_obj@W_vects
  t<-df
  }
  return(t)


}



