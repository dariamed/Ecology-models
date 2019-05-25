
A<-diag(5)
M<- (1/10^8)*matrix(1,5,5)
A_s<-solveRcpp(A+M)
A_1<- solve(A)
g<- sum(diag(M%*%A_1))
A_ninv2<- A_s + (1/(1+g))*A_1%*%M%*%A_1
A_s<-solveRcpp(A_ninv2)


















data<-sim_data$FacCompSparseSp20

data <- list(
  Y = subset(data, select = -env),
  X = cbind(1, scale(poly(data$env, 2))),
  covx = cov(cbind(1, scale(poly(data$env, 2)))),
  K = 3,
  J = ncol(data) - 1,
  n = nrow(data),
  I = diag(ncol(data) - 1),
  df = ncol(data)
)

dat<- cbind(data$Y,data$X)
dat<-dat[,-6]
dat_tr<-melt(dat,id=c("1","2"))
colnames(dat_tr)<-c("env","env2","sp","presence")

x1 <- runif(100)
x2 <- runif(100)
y  <- sample.int(3 , 100 , replace = T)

ggplot( dat_tr[which(dat_tr$presence==1),] )+
  geom_text( aes(env ,env2 ,label=sp, colour = factor(sp)))



