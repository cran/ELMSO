

#' Main ELMSO Function
#'
#' This function allows you to allocate budget to a set of websites based on the cost curve of the websites
#' and a matrix of pageviews for those sites.
#' @param z An n by p matrix of pageviews
#' @param CPM A p-dimensional vector of the average CPM values at each website. This is used to calculate the cost curve from a shifted logistic function. You may instead enter values for a p-dimensional "a" vector to define your own shifted logistic cost curve.
#' @param a A p-dimensional vector of values controlling the steepness of the shifted logistic cost curve. You may instead enter values for a p-dimensional vector of average CPM values to have the curve calculated for you.
#' @param tau A p-dimensional vector of total pageviews (in thousands) for each website. Defaults to the total pageviews in the matrix for each website (i.e., assumes z matrix represents all website pageviews) divided by 1000.
#' @param step A value to control the step size of the lambda grid (distance between budget points). Default is 0.05.
#' @param size A value to control the number of lambda values tried (number of budget points). Default is 100.
#' @param tol A value to control the convergence tolerance of the coordinate descent procedure. Default is 10^-3.
#' @param iters A value to control the number of iterations until algorithm should exit if convergence tolerance is not reached. Default is 200.
#' @return bid: A matrix of bid values by website at each budget
#' @return spend: a matrix of total spend by website at each budget
#' @return budget:  a vector of budget values
#' @return lambda: a vector of lambda values
#' @return a: a vector of a values (used to calculate shifted logistic curves and reach in reach.ELMSO function)
#' @references Courtney Paulson, Lan Luo, and Gareth M. James (2018) Efficient Large-Scale Internet Media Selection Optimization for Online Display Advertising. Journal of Marketing Research: August 2018, Vol. 55, No. 4, pp. 489-506.
#' @export
#' @examples
#' z=matrix(round(abs(rnorm(5000,0,0.7))),1000,5)
#' CPM.avg=c(3,4,5,6,7)
#' tau.values=rep(1000,5) #Note tau here is in thousands of pageviews
#'
#' allocation=ELMSO(z=z,CPM=CPM.avg,tau=tau.values)
#' allocation$bid
#' allocation$spend
#' allocation$budget
#' allocation$lambda
#' allocation$a

ELMSO <-
  function(z,CPM=NULL,a=NULL,tau=NULL,step=0.05,size=100,tol=10^-3,iters=200){


    if(is.null(CPM)&is.null(a)){stop("Must provide either average CPM values or logistic shaping parameter a")}

    B=0
    n=nrow(z)
    p=ncol(z)

    if(is.null(a)){a=log(3)/CPM}

    if(is.null(tau)){tau=round(colSums(z)/1000)}

    c.til=rep(0,p)
    c=rep(0,p)
    w.til=rep(0,p)
    w=rep(0,p)

    bid.c=w.mat=w.sum=lam.grid=NULL

    h=0


    for(j in 1:p){

      L1.temp = -2*a[j]*z[,j]*exp(a[j]*c.til[j]-z[,j]*(exp(a[j]*c.til[j])-1)/(exp(a[j]*c.til[j])+1))/(exp(a[j]*c.til[j])+1)^2
      L2.temp = -((a[j]*(exp(2*a[j]*c.til[j])+2*z[,j]*exp(a[j]*c.til[j])-1))/(exp(a[j]*c.til[j])+1)^2)*L1.temp

      H.temp=sum(c.til[j]*L2.temp-L1.temp)

      if(h<abs(H.temp)){h=abs(H.temp)}

    }

    lam=(2*h)^0.1
    lam.min=((2*h)^0.1-step*size)


    while(lam>lam.min){
      converge=1
      iter=0
      while(converge>tol & iter<iters){

        iter=iter+1

        for(j in 1:p){

          L1= -2*a[j]*z[,j]*exp(a[j]*c[j]-z[,j]*(exp(a[j]*c[j])-1)/(exp(a[j]*c[j])+1))/(exp(a[j]*c[j])+1)^2
          L2= -((a[j]*(exp(2*a[j]*c[j])+2*z[,j]*exp(a[j]*c[j])-1))/(exp(a[j]*c[j])+1)^2)*L1.temp

          H=sum(c.til[j]*L2-L1)
          sum.L=sum(L2)


          if(abs(H)>10^lam){
            c[j]=c.til[j]-(sum(L1)+10^lam*sign(H))/sum.L
            w[j]=c[j]*(exp(c[j]*a[j])-1)/(exp(c[j]*a[j])+1)*tau[j]
          }

          else {
            c[j]=w[j]=0
          }

        }

        if(sum(w)==0){
          converge=0
        }

        else{
          converge=sum(abs(w.til-w))/sum(w)
        }

        w.til=w
        c.til=c
        B=sum(w.til)


      }

      if(iter>iters){

        print("Max iterations reached. Converge value: ",converge)

      }

      bid.c=cbind(bid.c,c)
      w.mat=cbind(w.mat,w)
      w.sum=cbind(w.sum,sum(w))
      lam.grid=c(lam.grid,10^lam)
      lam=lam-step

    }

    list(bid=bid.c,spend=w.mat,budget=w.sum,lambda=lam.grid,a=a)

  }

#' Fixed ELMSO Function (fixed advertising costs, no cost curve)
#'
#' This function allows you to allocate budget to a set of websites when cost is fixed at each website based on
#' a matrix of pageviews for those sites.
#' @param z An n by p matrix of pageviews
#' @param CPM A p-dimensional vector of the (fixed) CPM values at each website
#' @param tau A p-dimensional vector of total pageviews (in thousands) for each website. Defaults to the total pageviews in the matrix for each website (i.e., assumes z matrix represents all website pageviews) divided by 1000.
#' @param step A value to control the step size of the lambda grid (distance between budget points). Default is 0.05.
#' @param size A value to control the number of lambda values tried (number of budget points). Default is 100.
#' @param tol A value to control the convergence tolerance of the coordinate descent procedure. Default is 10^-3.
#' @param iters A value to control the number of iterations until algorithm should exit if convergence tolerance is not reached. Default is 200.
#' @return spend: a matrix of total spend by website at each budget
#' @return budget:  a vector of budget values
#' @return lambda: a vector of lambda values
#' @references Courtney Paulson, Lan Luo, and Gareth M. James (2018) Efficient Large-Scale Internet Media Selection Optimization for Online Display Advertising. Journal of Marketing Research: August 2018, Vol. 55, No. 4, pp. 489-506.
#' @export
#' @examples
#' z=matrix(round(abs(rnorm(5000,0,0.7))),1000,5)
#' CPM.fixed=c(3,4,5,6,7)
#' tau.values=rep(100,5) #Note tau here is in thousands of pageviews
#'
#' allocation=ELMSO.fixed(z=z,CPM=CPM.fixed,tau=tau.values)
#' allocation$spend
#' allocation$budget
#' allocation$lambda

ELMSO.fixed  <-
  function(z,CPM,tau=NULL,step=0.05,size=100,tol=10^-3,iters=200){

    B=0
    n=nrow(z)
    p=ncol(z)

    if(is.null(tau)){tau=round(colSums(z)/1000)}

    w=rep(0,p)
    w.til=rep(0,p)

    w.mat=w.sum=lam.grid=NULL

    h=0

    for(j in 1:p){

      L1.temp = -1/(CPM[j]*tau[j])*z[,j]*exp(-1/(CPM[j]*tau[j])*z[,j]*w.til[j])
      L2.temp = -1/(CPM[j]*tau[j])*z[,j]*L1.temp

      H.temp=sum(w.til[j]*L2.temp-L1.temp)

      if(h<abs(H.temp)){h=abs(H.temp)}


    }

    lam=(2*h)^0.1
    lam.min=((2*h)^0.1-step*size)


    while(lam>lam.min){
      converge=1
      iter=0
      while(converge>tol & iter<iters){

        iter=iter+1

        for(j in 1:p){

          L1= -1/(CPM[j]*tau[j])*z[,j]*exp(-1/(CPM[j]*tau[j])*z[,j]*w.til[j])
          L2= -1/(CPM[j]*tau[j])*z[,j]*L1.temp

          H=sum(w.til[j]*L2-L1)
          sum.L=sum(L2)

          if(abs(H)>10^lam){
            w[j]=w.til[j]-(sum(L1)+10^lam*sign(H))/sum.L
          }

          else {
            w[j]=0
          }

        }

        if(sum(w)==0){
          converge=0
        }

        else{
          converge=sum(abs(w.til-w))/sum(w)
        }

        w.til=w
        B=sum(w)


      }

      if(iter>iters){

        print("Max iterations reached. Converge value: ",converge)

      }

      w.mat=cbind(w.mat,w)
      w.sum=cbind(w.sum,sum(w))
      lam.grid=c(lam.grid,10^lam)
      lam=lam-step

    }

    list(spend=w.mat,budget=w.sum,lambda=lam.grid)

  }

#' Calculating Reach from Main ELMSO Function
#'
#' This function allows you to calculate reach achieved at a given budget value from the ELMSO output.
#' @param bid A p-dimensional vector of the bidded CPM at each website for a particular budget value
#' @param a A p-dimensional vector of steepness values for the cost curves associated with each website
#' @param z An n by p matrix of pageviews
#' @return A value between 0 and 1 specifying the reach achieved with the given budget allocation.
#' @references Courtney Paulson, Lan Luo, and Gareth M. James (2018) Efficient Large-Scale Internet Media Selection Optimization for Online Display Advertising. Journal of Marketing Research: August 2018, Vol. 55, No. 4, pp. 489-506.
#' @export
#' @examples
#' z=matrix(round(abs(rnorm(5000,0,0.7))),1000,5)
#' CPM.avg=c(3,4,5,6,7)
#' tau.values=rep(100,5) #Note tau here is in thousands of pageviews
#'
#' allocation=ELMSO(z=z,CPM=CPM.avg,tau=tau.values)
#' reach.ELMSO(allocation$bid[,101],allocation$a,z)

reach.ELMSO<-function(bid,a,z) {

  p=ncol(z)
  gam=z

  for(j in 1:p){
    gam[,j]=(exp(a[j]*bid[j])-1)/(exp(a[j]*bid[j])+1)*z[,j]
  }

  1-1/nrow(z)*sum(exp(-rowSums(gam)))

}

#' Calculating Reach from Fixed ELMSO Function
#'
#' This function allows you to calculate reach achieved at a given budget value from the fixed ELMSO output.
#' @param CPM A p-dimensional vector of the fixed CPM at each website for a particular budget value
#' @param w A p-dimensional vector of amount spent at each website
#' @param z An n by p matrix of pageviews
#' @param tau A p-dimensional vector of total pageviews (in thousands) for each website. Defaults to the total pageviews in the matrix for each website (i.e., assumes z matrix represents all website pageviews) divided by 1000.
#' @return A value between 0 and 1 specifying the reach achieved with the given budget allocation.
#' @references Courtney Paulson, Lan Luo, and Gareth M. James (2018) Efficient Large-Scale Internet Media Selection Optimization for Online Display Advertising. Journal of Marketing Research: August 2018, Vol. 55, No. 4, pp. 489-506.
#' @export
#' @examples
#' z=matrix(round(abs(rnorm(5000,0,0.7))),1000,5)
#' CPM.fixed=c(3,4,5,6,7)
#' tau.values=rep(100,5) #Note tau here is in thousands of pageviews
#'
#' allocation=ELMSO.fixed(z=z,CPM=CPM.fixed,tau=tau.values)
#' reach.ELMSO.fixed(CPM.fixed,allocation$spend[,101],z,tau.values)
#'
#'

reach.ELMSO.fixed<-function(CPM,w,z,tau=NULL) {

  if(is.null(tau)){tau=round(colSums(z)/1000)}

  p=ncol(z)
  gam=z

  for(j in 1:p){
    gam[,j]=1/(CPM[j]*tau[j])*z[,j]*w[j]
  }

  1-1/nrow(z)*sum(exp(-rowSums(gam)))

}
