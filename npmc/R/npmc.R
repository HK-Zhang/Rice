"npmc" <-
function(dataset, control=NULL, df=2, alpha=0.05)
{
  mvtnorm <- require(mvtnorm, quietly=TRUE);
  if (!mvtnorm)
  {
    msg <- paste("npmc requires the mvtnorm-package to calculate",
                 "p-values for the test-statistics. This package is not",
                 "available on your system, so these values and the",
                 "confidence-limits will be missing values in the output!\n",
                 sep="\n");
    warning(msg);
  }
  
  if (any(df==0:2)) 
    dfm <- c(3,2,1)[df+1]
  else
  {
    warning("Wrong value for df\nUsing Satterthwaite t-approximation\n");
    dfm <- 1;
  }

  if (alpha<=0 || alpha>=1)
    stop("alpha must be a value between 0 and 1");

  
  name <- attr(dataset, "name");
  desc <- attr(dataset, "description");


  ##=== Function definitions ===================================================
  
  
  ##
  ## ssq:
  ## ----
  ## Calculates a vector's sum of squares
  ##
  ssq <- function(x) sum(x*x);

  ##
  ## force.ps: 
  ## ---------
  ## Forces a matrix to be positive semidefinite by replacing 
  ## all negative eigenvalues by zero.
  ## 
  force.ps <- function(M.in)
  {
    eig <- eigen(M.in, symmetric=TRUE);
    spec <- eig$values;
    if (adjusted <- any(spec<0))
    {
      spec[spec<0] <- 0;
      M <- eig$vectors %*% diag(spec) %*% t(eig$vectors);
      ginv <- diag(1/sqrt(diag(M)));
      M.out <- ginv %*% M %*% ginv;
      ##if ((msg <- all.equal(M.in,M.out))!=TRUE) attr(M.out, "msg") <- msg;
    }
    else
    {
      M.out <- M.in;
    }
    attr(M.out,"adjusted") <- adjusted; 
    return (M.out);
  }

  
  ##
  ## z.dist:
  ## -------
  ## Calculates the p-values for the teststatistics using the mvtnorm-package.
  ## The 'sides' parameter determines whether the p-values for the one-
  ## or two-sided test are calculated.
  ## The statistic is supposed to follow a multivariate t-statistic with
  ## correlation-matrix 'corr' and 'df' degrees of freedom. If df=0, the
  ## multivariate normal-distribution is used.
  ## We use the mvtnorm-package by Frank Bretz (www.bioinf.uni-hannover.de) 
  ## to calculate the corresponding p-values. These algorithms originally 
  ## calculate the value P(X<=stat) for the mentioned multivariate distributions, 
  ## i.e. the 1-sided p-value. In order to gain the 2-sided p-value as well, 
  ## we used the algorithms on the absolute value of the teststatistic in 
  ## combination with the inflated correlation-matrix
  ##   kronecker(matrix(c(1,-1,-1,1),ncol=2), corr)
  ##
  z.dist <- function(stat, corr, df=0, sides=2)
  {
    if (!mvtnorm) return (NA);
    
    if (sides==2)
    {
      corr <- kronecker(matrix(c(1,-1,-1,1),ncol=2), corr);
      stat <- abs(stat);
    }
    n <- ncol(corr);
    sapply(stat, function(arg) 
         {
           pmvt(
               lower=rep(-Inf, n), 
               upper=rep(arg, n), 
               df=as.integer(df), 
               corr=corr, 
               delta=rep(0,n)
               )[1];
         });     
  } 


  ##
  ## z.quantile:
  ## -----------
  ## Calculates the corresponding quantiles of z.dist p-values
  ## (used for the confidence-intervals)
  ##
  z.quantile <- function(p=0.95, start=0, corr, df=0, sides=2)
  {
    if (!mvtnorm) return (NA);

    if (z.dist(start,corr=corr,df=df,sides=sides) < p)
    {
      lower <- start;
      upper <- lower+1;
      while(z.dist(upper,corr=corr,df=df,sides=sides) < p)
        upper <- upper+1;
    }
    else
    {
      upper <- start;
      lower <- upper-1;
      while(z.dist(lower,corr=corr,df=df,sides=sides) > p)
        lower <- lower-1;
    }
    ur <- uniroot(f=function(arg) p-z.dist(arg,corr=corr,df=df,sides=sides),
                  upper=upper, lower=lower
                  );
    ur$root;
  }
  

  ##=== Calculations ===========================================================
      
  ## sort the dataset by factor
  dataset$class <- factor(dataset$class);
  datord <- order(dataset$class);
  attrs <- attributes(dataset);
  dataset <- data.frame(lapply(dataset, "[", datord));
  attributes(dataset) <- attrs;
  
  ## general data characteristics
  #attach(dataset);
  fl <- levels(dataset$class);             # factor-levels
  a <- nlevels(dataset$class);             # number of factor-levels
  samples <- base::split(dataset$var, dataset$class);    # split the data in separate sample-vectors
  n <- sapply(samples, length);    # sample sizes
  #detach(dataset);

  if (is.null(control))
  {
    ## create indexing vectors for the all-pairs situation
    tmp <- expand.grid(1:a, 1:a);
    ind <- tmp[[1]] > tmp[[2]];
    vi <- tmp[[2]][ind];
    vj <- tmp[[1]][ind];
  }
  else
  {
    ## create indexing vectors for the many-to-one situation
    if (!any(fl==control))
    {
      msg <- paste("Wrong control-group specification\n",
                   "The data does not contain a group with factor-level ",
                   control,
                   sep="");
      stop(msg, FALSE);
    }
    cg <- which(fl==control);
    vi <- which((1:a)!=cg);
    vj <- rep(cg, a-1);
  }

  ## number of comparisons ( a*(a-1)/2 for all-pairs, (a-1) for many-to-one )
  nc <- length(vi);              
  
  ## labels describing the compared groups 
  cmpid <- paste(vi, "-", vj, sep="");
  
  ## pairwise pooled sample-sizes
  gn <- n[vi]+n[vj];
  
  ## internal rankings
  intRanks <- lapply(samples, rank);
  
  ## pairwise rankings
  pairRanks <- lapply(1:nc, function(arg) 
                    {
                      rank(c(samples[[vi[arg]]], samples[[vj[arg]]]));  
                    });
 
  ## estimators for the relative effects
  pd <- sapply(1:nc, function(arg)
             {
               i <- vi[arg]; 
               j <- vj[arg];
               (sum(pairRanks[[arg]][(n[i]+1):gn[arg]])/n[j]-(n[j]+1)/2)/n[i];  
             });
    
  ## Calculations for the BF-test ###################################
  ##
  dij <- dji <- list(0);

  sqij <- sapply(1:nc, function(arg) 
               {
                 i <- vi[arg]; 
                 j <- vj[arg];
                 pr <- pairRanks[[arg]][(n[i]+1):gn[arg]];
                 dij[[arg]] <<- pr - sum(pr)/n[j] - intRanks[[j]] + (n[j]+1)/2;
                 ssq(dij[[arg]])/(n[i]*n[i]*(n[j]-1));
               });
  
  sqji <- sapply(1:nc, function(arg)
               {
                 i <- vi[arg];  
                 j <- vj[arg];
                 pr <- pairRanks[[arg]][1:n[i]];
                 dji[[arg]] <<- pr - sum(pr)/n[i] - intRanks[[i]] + (n[i]+1)/2;
                 ssq(dji[[arg]])/(n[j]*n[j]*(n[i]-1));
               });

  ## diagonal elements of the covariance-matrix
  vd.bf <- gn*(sqij/n[vj] + sqji/n[vi]);

  ## mark and correct zero variances for further calculations
  singular.bf <- (vd.bf==0);
  vd.bf[singular.bf] <- 0.00000001;
  
  ## standard-deviation
  std.bf <- sqrt(vd.bf/gn);

  ## teststatistic
  t.bf <- (pd-0.5)*sqrt(gn/vd.bf);
  
  ## Satterthwaite approxiamtion for the degrees of freedom
  df.sw <- (n[vi]*sqij + n[vj]*sqji)^2 / 
    ((n[vi]*sqij)^2/(n[vj]-1) + (n[vj]*sqji)^2/(n[vi]-1));
  df.sw[is.nan(df.sw)] <- Inf;

  ## choose degrees of freedom 
  df <- if (dfm<3) max(1, if (dfm==2) min(gn-2) else min(df.sw)) else 0;


  ## Calculations for the Steel-test ################################
  ##
  ## the Steel-type correlation factors
  lambda <- sqrt(n[vi]/(gn+1));
      
  ## diagonal elements of the covariance-matrix
  vd.st <- sapply(1:nc, function(arg) ssq(pairRanks[[arg]]-(gn[arg]+1)/2)) / 
    (n[vi]*n[vj]*(gn-1));

  ## mark and correct zero variances for further calculations
  singular.st <- (vd.st==0);
  vd.st[singular.st] <- 0.00000001;
  
  ## standard-deviation
  std.st <- sqrt(vd.st/gn);

  ## teststatistic
  t.st <- (pd-0.5)*sqrt(gn/vd.st);
  

  ## Calculate the correlation-matrices (for both, BF and Steel) ####
  ##    
  rho.bf <- rho.st <- diag(nc);
  for (x in 1:(nc-1))
  {
    for (y in (x+1):nc)
    {
      i <- vi[x]; j <- vj[x];
      v <- vi[y]; w <- vj[y];
      p <- c(i==v, j==w, i==w, j==v);
      if (sum(p)==1) 
      {      
        cl <- list(
                   function()  (t(dji[[x]]) %*% dji[[y]]) / (n[j]*n[w]*n[i]*(n[i]-1)),
                   function()  (t(dij[[x]]) %*% dij[[y]]) / (n[i]*n[v]*n[j]*(n[j]-1)),
                   function() -(t(dji[[x]]) %*% dij[[y]]) / (n[i]*n[w]*n[i]*(n[i]-1)),
                   function() -(t(dij[[x]]) %*% dji[[y]]) / (n[j]*n[v]*n[j]*(n[j]-1))
                   );
        case <- (1:4)[p];
        rho.bf[x,y] <- rho.bf[y,x] <- 
          sqrt(gn[x]*gn[y]) / sqrt(vd.bf[x]*vd.bf[y]) * cl[[case]]()
        ;
        rho.st[x,y] <- rho.st[y,x] <- 
        {if (case>2) -1 else 1}*lambda[x]*lambda[y]
        ;
      }
    }
  }
  rho.bf <- force.ps(rho.bf);
  rho.st <- force.ps(rho.st);
      

  ## Calculate the p-values     (BF and Steel) ######################
  ##
  p1s.bf <- 1 - z.dist(t.bf, corr=rho.bf, df=df, sides=1);
  p2s.bf <- 1 - z.dist(t.bf, corr=rho.bf, df=df, sides=2);
   
  p1s.st <- 1 - z.dist(t.st, corr=rho.st, sides=1);
  p2s.st <- 1 - z.dist(t.st, corr=rho.st, sides=2);

  
  ## Calculate the confidence-limits (BF and Steel) #################
  ##
  z.bf <- z.quantile(1-alpha, corr=rho.bf, df=df, sides=2);
  lcl.bf <- pd - std.bf*z.bf;
  ucl.bf <- pd + std.bf*z.bf;

  z.st <- z.quantile(1-alpha, corr=rho.st, sides=2);
  lcl.st <- pd - std.st*z.st;
  ucl.st <- pd + std.st*z.st;

  
  ##=== Output ==================================================================
      
  ## Create the result-datastructures ###############################
  ##    
  dataStructure <- data.frame("group index"=1:a, 
                              "class level"=fl, 
                              "nobs"=n
                              );
  
  test.bf <- data.frame("cmp"=cmpid, 
                        "gn"=gn, 
                        "effect"=pd,
                        "lower.cl"=lcl.bf,
                        "upper.cl"=ucl.bf,
                        "variance"=vd.bf, 
                        "std"=std.bf, 
                        "statistic"=t.bf, 
                        "p-value 1s"=p1s.bf, 
                        "p-value 2s"=p2s.bf, 
                        "zero"=singular.bf
                        ); 

  test.st <- data.frame("cmp"=cmpid, 
                        "gn"=gn, 
                        "effect"=pd, 
                        "lower.cl"=lcl.st,
                        "upper.cl"=ucl.st,
                        "variance"=vd.st, 
                        "std"=std.st, 
                        "statistic"=t.st, 
                        "p-value 1s"=p1s.st, 
                        "p-value 2s"=p2s.st, 
                        "zero"=singular.st
                        ); 

  result <- list("data"=dataset,
                 "info"=dataStructure, 
                 "corr"=list("BF"=rho.bf, "Steel"=rho.st),
                 "test"=list("BF"=test.bf, "Steel"=test.st),
                 "control"=control,
                 "df.method"=dfm,
                 "df"=df,
                 "alpha"=alpha
                 );
  
  class(result) <- "npmc";
  
  return (result);
  
}
