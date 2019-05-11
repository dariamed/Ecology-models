

A<-diag(5)
M<- (1/10^8)*matrix(1,5,5)

A_s<-solveRcpp(A+M)
A_1<- solve(A)
g<- sum(diag(M%*%A_1))
A_ninv2<- A_s + (1/(1+g))*A_1%*%M%*%A_1




A_s<-solveRcpp(A_s)
