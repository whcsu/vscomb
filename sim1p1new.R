source("/u/home/w/wanghong/sisnew/screenfun.R")
  set.seed(123)
#########################################
## Example for uncensored dataset
#### Following Meier, Geer, and Bühlmann (2009), we generate the data from the following additive model:
 sim_dat1<-function(N, p,rho=0.5){
  active=c(1:4)
  sig=diag(p) ## covariance matrix
  sig <- rho^abs(row(sig)-col(sig));
  g1<-function(x) (x)
  g2<-function(x) ((2*x-1)^2)
  g3<-function(x) (sin(2*pi*x)/(2-sin(2*pi*x)))
  g4<-function(x)(0.1*sin(2*pi*x)+0.2*cos(2*pi*x)+0.3*sin(2*pi*x)^2+0.4*cos(2*pi*x)^3+0.5*sin(2*pi*x)^3)
  X=mvrnorm(N,mu=rep(0,p),Sigma=sig)
  Y=0.5*g1(X[,1])+3*g2(X[,2])+4*g3(X[,3])+16*g4(X[,4])+sqrt(1.74)*rnorm(N,0,1) 
  beta <- rep(0,p);
  beta[1:4]=c(0.5,3,4,16)
  return(list(x=X,y=Y,active=active,beta=beta))
 
}
#Case 2 we use t(4) distributiont to generate predictors
sim_dat2<-function(N, p){
  active=c(1:4)
  X<-matrix(rt(n*p,4),n,p)
  Y<- -3*sin(2*X[,1])+(X[,2]*X[,2]-25.0/12)-1.5*X[,3]+(exp(-X[,4])-2.0/5*sinh(5/2))+rt(n,1)
   return(list(x=X,y=Y,active=active))
}

sim_dat3<-function(N, p){
  active=c(1:4)
 #Case 3 uniform(-2.5,2.5) distribution to generate predictors
 X<-matrix(runif(n*p,-2.5,2.5),n,p)
 Y<- -3*sin(2*X[,1])+(X[,2]*X[,2]-25.0/12)-1.5*X[,3]+(exp(-X[,4])-2.0/5*sinh(5/2))+rt(n,1)
   return(list(x=X,y=Y,active=active))
}
## https://getd.libs.uga.edu/pdfs/qiu_debin_201605_phd.pdf
## page 64

sim_dat4<-function(N, p,rho=0.6){
  active=c(1:15)
 #Case 4: The grouped variable selection
p0 <- 5 # number of predictors in each group
J <- floor(p/p0) # group size
group <- rep(1:J,each = 5) # group indices
##autoregressive correlation
Sigma <- rho^abs(matrix(1:(p0*J),p0*J,p0*J) - t(matrix(1:(p0*J),p0*J,p0*J)))
X <- mvrnorm(N,seq(0,5,length.out = p0*J),Sigma)
betaTrue <- runif(15,-2,5)
mu <- X%*%matrix(c(betaTrue,rep(0,p0*J-15)),ncol = 1)
# normal distribution
Y <- mu + rnorm(N)

   return(list(x=X,y=Y,active=active))
}


reptime=100
nv=c(200)
p=1000
myrhopv=c(0.5,0.6,0.7,0.8,0.9)


for (k in 1:length(nv)){
    n=nv[k]
	for (j in 1:length(myrhopv)){
	    myrho=myrhopv[j]
	  	SIS=NULL
		  RRCS=NULL
	  	HOLP=NULL
		  RF=NULL
		  R_S=NULL
		  R_R=NULL
		   R_H=NULL
		for (i in 1:reptime){
        
		#dat=sim_dat4(N=n,p=p) 
		dat=sim_dat4(N=n,p=p,rho=myrho)
		x=dat$x
		y=dat$y

		active=dat$active

		num_select=floor(dim(x)[1]/log(dim(x)[1]))
		#if(num_select%%2!=0) num_select=num_select+1

		ss1=screening(x,y,num.select=num_select,method="sis")
		ss2=screening(x,y,num.select=num_select,method="rrcs")
		ss3=screening(x,y,num.select=num_select,method="holp")
		    rfs=Boruta(x,y); 
			ssall=order(attStats(rfs)$medianImp,decreasing = T)
		ss4=ssall[1:num_select]
		#ss4_half=ssall[1:(num_select/2)]

		  comb1=orderedunion(ss1$screen[1:(num_select)],ss4)
		  comb2=orderedunion(ss2$screen[1:(num_select)],ss4)
		  comb3=orderedunion(ss3$screen[1:(num_select)],ss4)
		  
		  SIS=cbind(SIS,ss1$screen)
		  RRCS=cbind(RRCS,ss2$screen)
		  HOLP=cbind(HOLP,ss3$screen)
		  RF=cbind(RF,ss4)
		  R_S=cbind(R_S,comb1)
		  R_R=cbind(R_R,comb2)
		  R_H=cbind(R_H,comb3)  

		}


		A1=Sel.rate(n,active,SIS) 
		A2=Sel.rate(n,active,RRCS) 
		A3=Sel.rate(n,active,HOLP) 
		A4=Sel.rate(n,active,RF) 
		A5=Sel.rate(n,active,R_S) 
		A6=Sel.rate(n,active,R_R)
		A7=Sel.rate(n,active,R_H)

		B1=PA.rate(active,SIS) 
		B2=PA.rate(active,RRCS) 
		B3=PA.rate(active,HOLP) 
		B4=PA.rate(active,RF) 
		B5=PA.rate(active,R_S) 
		B6=PA.rate(active,R_R)
		B7=PA.rate(active,R_H)
		
		C1=PC.rate(active,SIS) 
		C2=PC.rate(active,RRCS) 
		C3=PC.rate(active,HOLP) 
		C4=PC.rate(active,RF) 
		C5=PC.rate(active,R_S) 
		C6=PC.rate(active,R_R)
		C7=PC.rate(active,R_H)


		#boxplot(A1,A2,A3,A4,A5,A6,A7)

		resulta=data.frame(A1,A2,A3,A4,A5,A6,A7)
		colnames(resulta)=c("SIS","RRCS","HOLP","RF","R_S","R_R","R_H")
		csvfn1=paste0("n=",n,"p=",p,"rho=",myrho,"1-single.csv")
		write.csv(resulta,csvfn1)
		resultb=c(B1,B2,B3,B4,B5,B6,B7)
		names(resultb)=c("SIS","RRCS","HOLP","RF","R_S","R_R","R_H")
		csvfn2=paste0("n=",n,"p=",p,"rho=",myrho,"1-all.csv")
		write.csv(resultb,csvfn2)
		
		resultc=data.frame(C1,C2,C3,C4,C5,C6,C7)
		colnames(resultc)=c("SIS","RRCS","HOLP","RF","R_S","R_R","R_H")
		csvfn3=paste0("n=",n,"p=",p,"rho=",myrho,"1-PC.csv")
		write.csv(resultc,csvfn3)
		
		
		csvfn=paste0("n=",n,"p=",p,"rho=",myrho,"1.RData")
		save.image(csvfn)

	}#pho
}#n
# A1=PA(active,SIS) 
# A2=PA(active,RRCS) 
# A3=PA(active,HOLP) 
# A4=PA(active,RF) 
# A5=PA(active,R_S) 
# A6=PA(active,R_R)
# A7=PA(active,R_H)
