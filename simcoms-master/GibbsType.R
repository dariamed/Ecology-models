

rm(list=ls())
library(pracma)
library(Brobdingnag)
library(copula)
#####Prior distribution on the number of groups.
kval<-10
theta<-1 
sigma<-0.5

############### V for PY########################
v_py<- function(kval, sigma,theta,npoints){
  c_v<-1:(kval-1)
  v_nk<- (theta +sigma*c_v)
  Vup<- prod(v_nk)
  n_vec<- 1:(npoints-1)
  Vlow<- prod(theta +n_vec)
  V_t_nk<-Vup/Vlow
 return(V_t_nk)
}

###############Generalized coefficient########################

gen_fac_coef<-function(kval,sigma,npoints){
  sum<-0
  kfac<-factorial(kval)
  for(i in (0:kval)){
    n_vec<- 0:(npoints-1)
    sn<- prod(-i*sigma +n_vec)
    ckn<- choose(kval,i)
    sum<- sum + ((-1)^i)*ckn*sn 
  }
  sumf<- sum/kfac
  return(sumf)
}


############### V for NGG########################

v_ng<- function(beta, sigma, kval, npoints){
  sum<-0
  coef_low<-gamma(npoints)
  coef_high<-exp(beta)* sigma^(kval-1)
  coef<- coef_high/coef_low
  for(i in (0:(npoints-1))){
    gn<- incgam(beta,kval - i/sigma)
    ckn<- choose(npoints-1,i)
    sum<- sum + ((-1)^i)*(beta^(i/sigma))*ckn*gn
  }
  sumf<- sum/coef
  return(sumf) 
}

########### density for  NGG #############################


prob_ng<- function(kg, npoints, sigma, beta){
  pb_v<- v_ng(beta, sigma, kg, npoints)
  pb_gen<- gen_fac_coef(kg,sigma, npoints)
  prob<- (pb_v*pb_gen)/(sigma^kg)
  return(prob)
}

########### density for PY#############################


prob_py<- function(kg, npoints, sigma, theta){
  pb_v<- v_py(kg,sigma, theta,npoints)
  pb_gen<- gen_fac_coef(kg,sigma, npoints)
  prob<- (pb_v*pb_gen)/(sigma^kg)
  return(prob)
}

########### density for Dirichlet#############################
prob_dir<- function(k, npoints, theta){
  n_vec<- 0:(npoints-1)
  theta_n<- prod(theta +n_vec)
  prob<- ((theta^k) *(abs(Stirling1(npoints,k))))/theta_n
  return(prob)
}


prob_dir_large_dim<- function(k, npoints, theta){
  n_vec<-as.brob( 0:(npoints-1))
  theta_n<- prod(theta +n_vec)
  stir<- as.brob(abs(Stirling1(npoints,k)))
  powerk<- as.brob((theta^k))
  prob_brob<- powerk*(stir/theta_n)
  prob<- as.numeric(prob_brob)
  return(prob)
}


n_vec<-as.brob( 0:(npoints-1))
theta_n<- prod(theta +n_vec)
stir<- as.brob(abs(Stirling1(npoints,k)))
powerk<- as.brob((theta^k))
prob_brob<- powerk*(stir/theta_n)
prob<- as.numeric(prob_brob)



####################################################################################################################
#############################Plotting###############################################################################
k_vec<-seq(1,50,by=49/9)
sigma_vec<-seq(0.2,0.8, by=0.6/9)
z<- outer(k_vec,sigma_vec,Vectorize(prob_py),npoints=50, theta=1)

# p<- plot_ly(showscale = TRUE) %>%
#   add_surface(x=k_vec, y=sigma_vec,z =z, cmin = min(z), cmax = max(z),colorbar=list(title='PY'), colorscale = list(c(0,1),c("rgb(255,112,184)","rgb(128,0,64)")),opacity = 0.98) %>%
#   layout(title="Prior distribution", scene = list(xaxis= list(title="K"),yaxis= list(title="sigma"),zaxis= list(title="N",range = c(min(z),max(z)))))
# p
# 
# k_vec<- seq(1,50, by =1)
# y<- c()
# for(l in (1: length(k_vec))){
#   y<- c(y, prob_py(k_vec[l],npoints=50,theta=1,sigma=0.15))
# }
# 
# plot(k_vec, y)
# k_vec<-seq(1,50,by=49/9)
# sigma_vec<-seq(0.2,0.8, by=0.6/9)
# z<- outer(k_vec,sigma_vec,Vectorize(prob_py),npoints=50, theta=1)
# 

#######################################PY##################################
sigma_df <- rep(seq(0.2,0.8, by=0.6/9), each=50)
k_df<- rep(1:50, 10)
df<- as.data.frame(matrix(NA, nrow=500, ncol=1))

df$k<- k_df
df$sigma<- sigma_df
for(l in (1: nrow(df))){
  df$val[l]<- prob_py(df$k[l],npoints=50,theta=1,sigma=df$sigma[l])
}

df$sigma<- as.factor(df$sigma)
p<- plot_ly(df, x =df$k , y = df$sigma, z = df$val, split = df$sigma, type = "scatter3d", mode = "lines") %>% 
  layout(title="Prior density", scene = list(xaxis= list(title="K"),yaxis= list(title="sigma"), zaxis= list(title="Probability")))                                                                                          


p

#####################################################################################################################

sigma_df <- rep(seq(0.2,0.8, by=0.6/9), each=50)
k_df<- rep(1:50, 10)
df<- as.data.frame(matrix(NA, nrow=500, ncol=1))

df$k<- k_df
df$sigma<- sigma_df
for(l in (1: nrow(df))){
  df$val[l]<- prob_ng(df$k[l],npoints=50,beta=1,sigma=df$sigma[l])
}

df$sigma<- as.factor(df$sigma)
p<- plot_ly(df, x =df$k , y = df$sigma, z = df$val, split = df$sigma, type = "scatter3d", mode = "lines") %>% 
  layout(title="Prior density", scene = list(xaxis= list(title="K"),yaxis= list(title="sigma"), zaxis= list(title="Probability")))                                                                                          


p

##########Dirichlet
#prob_dir<- function(k, npoints, theta)


theta_df <- rep(seq(1,100, by=99/9), each=150)
k_df<- rep(1:150, 10)
df<- as.data.frame(matrix(NA, nrow=1500, ncol=1))

df$k<- k_df
df$theta<- theta_df
for(l in (1: nrow(df))){
  df$val[l]<- prob_dir_large_dim(df$k[l],npoints=150,theta=df$theta[l])
}

df$theta<- as.factor(df$theta)
p150<- plot_ly(df, x =df$k , y = df$theta, z = df$val, split = df$theta, type = "scatter3d", mode = "lines") %>% 
  layout(title="Prior density for Dirichlet prior for n=50", scene = list(xaxis= list(title="K"),yaxis= list(title="alpha"), zaxis= list(title="Probability")))                                                                                          


p150



#prob_dir(10,npoints=150,theta=60)


###################
# 
# StirlingFirst<- function(n, k){
#   if((k == 0)&(n != 0)){
#     return(0)
#   }
#   else{
#     if(k==n){return(1)}
#     else{
#       return((n-1)*StirlingFirst(n-1, k) + StirlingFirst(n-1, k-1))
#       }
#   }
# }


k<-60
npoints<-150
n_vec<-as.brob( 0:(npoints-1))
theta_n<- prod(theta +n_vec)
stir<- as.brob(abs(Stirling1(npoints,k)))
powerk<- as.brob((theta^k))
prob_brob<- powerk*(stir/theta_n)
prob<- as.numeric(prob_brob)


x <- as.brob(1:10)
y <- 1e10
x+y
as.numeric((x+y)-1e10)
x^(1/y)


