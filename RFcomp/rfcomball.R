source("/u/home/w/wanghong/sisnew/screenfun.R")
  set.seed(123)
  library(Pomona)
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
p=5000
myrhopv=c(0.6)


for (k in 1:length(nv)){
    n=nv[k]
	for (j in 1:length(myrhopv)){
	    myrho=myrhopv[j]
	  	  DC=NULL
        RF1=NULL
		    RF2=NULL
	  	  RF3=NULL
        R_1=NULL
		    R_2=NULL
		    R_3=NULL
		for (i in 1:reptime){
        
		dat=sim_dat3(N=n,p=p) 
		#dat=sim_dat1(N=n,p=p,rho=myrho)
		x=dat$x
		y=dat$y

		active=dat$active

		num_select=floor(dim(x)[1]/log(dim(x)[1]))
		#if(num_select%%2!=0) num_select=num_select+1

       ss4=screenIID(x,y, method="DC-SIS")
       
	   rfs=Boruta(x,y); 
	   ssall=order(attStats(rfs)$medianImp,decreasing = T)
       rf1=ssall[1:num_select]
	   
	  
      #perm
      varsel=var.sel.perm(x=x,y=y,no.perm = 10)
       # DELETE X FROM VAR
      rfvar1=as.numeric(gsub("[a-zA-Z ]", "", varsel$var))
      rf2=rfvar1[1:num_select]

      #vita
       varsel=var.sel.vita(x=x,y=y)
       # DELETE X FROM VAR
       rfvar1=as.numeric(gsub("[a-zA-Z ]", "", varsel$var))
       rf3=rfvar1[1:num_select]

       #Boruta-DCSIS
       comb1=orderedunion(ss4$rank[1:(num_select)],rf1)
	   	
	   #Perm-DCSIS
	   comb2=orderedunion(ss4$rank[1:(num_select)],rf2)
	   
	   #Vita-DCSIS
	   comb3=orderedunion(ss4$rank[1:(num_select)],rf3)
	   
	
   
	   
	  DC=	cbind(DC,ss4$rank[1:(num_select)])
		RF1=cbind(RF1,rf1)
		RF2=cbind(RF2,rf2)
		RF3=cbind(RF3,rf3)
    R_1=cbind(R_1,comb1)
		R_2=cbind(R_2,comb2)
		R_3=cbind(R_3,comb3)  
	 

	   
		  
		  

		}


	  A0=Sel.rate(n,active,DC) 
		A1=Sel.rate(n,active,RF1) 
		A2=Sel.rate(n,active,RF2) 
		A3=Sel.rate(n,active,RF3)
		A4=Sel.rate(n,active,R_1)
    A5=Sel.rate(n,active,R_2)
		A6=Sel.rate(n,active,R_3)
		
 	  B0=PA.rate(active,DC)     
		B1=PA.rate(active,RF1) 
		B2=PA.rate(active,RF2) 
		B3=PA.rate(active,RF3)
		B4=PA.rate(active,R_1)
    B5=PA.rate(active,R_2)
		B6=PA.rate(active,R_3)
	
	  C0=PC.rate(active,DC) 
		C1=PC.rate(active,RF1) 
		C2=PC.rate(active,RF2) 
		C3=PC.rate(active,RF3)
		C4=PC.rate(active,R_1)
		C5=PC.rate(active,R_2)
		C6=PC.rate(active,R_3)


		#boxplot(A1,A2,A3,A4,A5,A6,A7)

		resulta=data.frame(A0,A1,A2,A3,A4,A5,A6)
		colnames(resulta)=c("DC","RF1","RF2","RF3","R_1","R_2","R_3")
		csvfn1=paste0("RFn=",n,"p=",p,"rho=",myrho,"m3-single.csv")
		write.csv(resulta,csvfn1)
		resultb=c(B0,B1,B2,B3,B4,B5,B6)
		names(resultb)=c("DC","RF1","RF2","RF3","R_1","R_2","R_3")
		csvfn2=paste0("RFn=",n,"p=",p,"rho=",myrho,"m3-all.csv")
		write.csv(resultb,csvfn2)
		
		resultc=data.frame(C0,C1,C2,C3,C4,C5,C6)
		colnames(resultc)=c("DC","RF1","RF2","RF3","R_1","R_2","R_3")
		csvfn3=paste0("RFn=",n,"p=",p,"rho=",myrho,"m3-PC.csv")
		write.csv(resultc,csvfn3)
		
		
		#csvfn=paste0("RFn=",n,"p=",p,"rho=",myrho,"m3.RData")
		#save.image(csvfn)

	}#pho
}#n
# A1=PA(active,SIS) 
# A2=PA(active,RRCS) 
# A3=PA(active,HOLP) 
# A4=PA(active,RF) 
# A5=PA(active,R_S) 
# A6=PA(active,R_R)
# A7=PA(active,R_H)
