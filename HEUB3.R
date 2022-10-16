# 'creates test statistic value and critical value
#' @export
#' @param x data matrix
#' @param alpha numeric variable
HEUB3=function(x,alpha)
{

  p=3
  k=3
  sequence_of_n=c(nrow(x[[1]]),nrow(x[[2]]),nrow(x[[3]]))
  y=vector(mode = "list", length = k)
  k_hat_0=vector(mode = "list", length = k)
  k_hat_hat_0=vector(mode = "list", length = k)
  mu_hat_0=vector(mode = "list", length = k)
  mu_hat_hat_0=vector(mode = "list", length = k)
  g=vector(mode = "list", length = k)
  z=vector(mode = "list", length = k)
  r_bar=array(0,k)
  r_bar1=array(0,k)
  q_bar=array(0,k)
  q_bar1=array(0,k)
  l=array(0,k)
  l1=array(0,k)

  n11=vector(mode = "list", length = k)
  n1=vector(mode = "list", length = k)
  n1_hat=vector(mode = "list", length = k)
  gamma1=array(0,k)

  Q1=vector(mode = "list", length = k)
  A=vector(mode = "list", length = k)
  # D1=diag(1,nrow = (p),ncol = (p))nbt
  # O=diag(0,nrow = (p),ncol = (p))
  # l1=list(D1=D1,D2=-D1,D3=O)
  # m=array(c("D1","D3","D2","D1","D3","D2"),c((k-1),k))
  # C=blockmatrix(value=m,list=l1)
  # C1=as.matrix(C)
  mu_0=matrix(c(1,1,-1)/sqrt((1)^2+2),nrow = p,ncol = 1)

  H=1000 # no of inner bootstrap simulation
  inner_loop_statistic_value_array=array(0,H)
  for(i in 1:k)
  {

    y[[i]]=c(sum(x[[i]][,1]),sum(x[[i]][,2]),sum(x[[i]][,3]))

    r_bar[i]=sqrt(sum(y[[i]]**2))/sequence_of_n[i]

    r_bar1[i]=sqrt(sum(y[[i]]**2))

    mu_hat_0[[i]]=y[[i]]/sqrt(sum(y[[i]]**2))
    k_hat_0[[i]]=matrix(c(mu_hat_0[[i]]),nrow = p,ncol = 1)
    l[i]=(r_bar[i]*(p-r_bar[i]^2))/(1-r_bar[i]^2)

  }
  #H1=rbind(x[[1]],x[[2]],x[[3]])
  X=matrix(c(unlist(mu_hat_0)),nrow = k*(p),ncol = 1,byrow = FALSE) #X^ value
  Y=matrix(c(0.5773503,0.5773503,-0.5773503,0.5773503,0.5773503,-0.5773503,0.5773503,0.5773503,-0.5773503),nrow = k*(p),ncol = 1,byrow = FALSE)
  #X1=(X-Y)
  #N=sum(sequence_of_n)    #calculations of Q^
  Q=c(sqrt(sequence_of_n[1])*(X[1:3,]-Y[1:3,]),sqrt(sequence_of_n[2])*(X[4:6,]-Y[4:6,]),sqrt(sequence_of_n[3])*(X[7:9,]-Y[7:9,]))
  #Q2=matrix(c(Q),nrow=2,ncol=2,byrow = T)
  #Q=sqrt(N)*(X1)
  ###################################################
  statistic_value<-sum(abs(Q)**2)

  M1= matrix(0,nrow = sequence_of_n[1],ncol = p)
  M2= matrix(0,nrow = sequence_of_n[2],ncol = p)
  M3= matrix(0,nrow = sequence_of_n[3],ncol = p)



  M=list(M1,M2,M3)

  for(i in 1:k)
  {
    n11[[i]]=(t(mu_0)%*%k_hat_0[[i]])
    gamma1[i]=acos(c(n11[[i]]))

    n1[[i]]=k_hat_0[[i]]-mu_0*c(n11[[i]])
    n1_hat[[i]]=n1[[i]]/sqrt(sum(n1[[i]]**2))

    A[[i]]=mu_0%*%t(n1_hat[[i]])-n1_hat[[i]]%*%t(mu_0)
    Q1[[i]]=diag(1,p)+sin(gamma1[i])*A[[i]]+(cos(gamma1[i])-1)*(mu_0%*%t(mu_0)+n1_hat[[i]]%*%t(n1_hat[[i]]))
    for(j in 1:sequence_of_n[i])
    {
      M[[i]][j,]=t(Q1[[i]]%*%matrix(c(x[[i]][j,]),nrow = p,ncol = 1))
    }
  }

  for(h in 1:H)
  {

    for(i in 1:k)
    {

      #g[[i]]=rmovMF(sequence_of_n[i],l[i]*b1)
      g[[i]]=M[[i]][sample(nrow(M[[i]]),replace = TRUE),]



      z[[i]]=c(sum(g[[i]][,1]),sum(g[[i]][,2]),sum(g[[i]][,3]))

      q_bar[i]=sqrt(sum(z[[i]]**2))/sequence_of_n[i]
      q_bar1[i]=sqrt(sum(z[[i]]**2))

      mu_hat_hat_0[[i]]=z[[i]]/q_bar1[i]
      k_hat_hat_0[[i]]=matrix(c(mu_hat_hat_0[[i]]),nrow = p,ncol = 1)

      #l1[i]=(q_bar[i]*(p-q_bar[i]^2))/(1-q_bar[i]^2)
      #N2[i]=t(k_hat_hat_0[[i]]-mu_0)%*%ginv(V/(l1[i]*a.p.k(p,l1[i])))%*%(k_hat_hat_0[[i]]-mu_0)

    }
    X2=matrix(c(unlist(mu_hat_hat_0)),nrow = k*(p),ncol = 1,byrow = FALSE)
    #X3=(X2-Y)
    N=sum(sequence_of_n)
    #Q2=sqrt(N)*X3
    Q2=c(sqrt(sequence_of_n[1])*(X2[1:3,]-Y[1:3,]),sqrt(sequence_of_n[2])*(X2[4:6,]-Y[4:6,]),sqrt(sequence_of_n[3])*(X2[7:9,]-Y[7:9,]))
    inner_loop_statistic_value_array[h]=sum(abs(Q2)**2)

    #inner_loop_statistic_value_array[h]=2*(l1[1]*q_bar1[1]*N2[1]+l1[2]*q_bar1[2]*N2[2]+l1[3]*q_bar1[3]*N2[3])
    #Q1=1-m2%*%mu_0
    #inner_loop_statistic_value_array[h]=1-m2%*%mu_0
    #inner_loop_statistic_value_array[h]=2*b3-2*sqrt(sum(b2**2))*(m2%*%mu_0)
  }
  critical_value=quantile(inner_loop_statistic_value_array,1-alpha)
  return(list(statistic_value,critical_value))

}
