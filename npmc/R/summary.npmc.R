"summary.npmc" <-
function(object, type="both", info=TRUE, short=TRUE, corr=FALSE, ...)
{
  x <- object;
  if (info)
  {
    name <- attr(data, "name");
    desc <- attr(data, "desc");
    df <- x$df;
    df.method <- x$df.method;
    alpha <- x$alpha;
  
    apm <- c(paste("Satterthwaite t-approximation (df=",df,")",sep=""),
             paste("simple t-approximation (df=",df,")",sep=""),
             "standard normal approximation"
             );
    msg <- c(paste("npmc executed", if (!is.null(name)) paste("on", name)),
             if (is.null(desc)) "" else c("","Description:",desc,""),
             "NOTE:",
             paste("-Used", apm[df.method]),
             paste("-Calculated simultaneous (1-", alpha, ") confidence intervals",sep=""),
             "-The one-sided tests 'a-b' reject if group 'a' tends to",
             " smaller values than group 'b'"
             );
    report(msg, style=2, char="/");
    report();
  }

  if (short)
  {
    bf <- st <- c("cmp","effect","lower.cl","upper.cl","p.value.1s","p.value.2s");
  }
  else
  {
    bf  <- names(x$test$BF);
    st <- names(x$test$Steel);
  }

  
  content <- list();
  if (info)
    content <- c(content, list("Data-structure"=x$info));
  if (corr && type!="Steel")
    content <- c(content, list("Behrens-Fisher type correlation-matrix"=x$corr$BF));
  if (corr && type!="BF")
    content <- c(content, list("Steel type correlation-matrix"=x$corr$Steel));
  if (type!="Steel")
    content <- c(content, list("Results of the multiple Behrens-Fisher-Test"=x$test$BF[bf]));
  if (type!="BF")
    content <- c(content, list("Results of the multiple Steel-Test"=x$test$Steel[st]));
  
  ##h <- (list("Data-structure"=x$info, 
  ##           "Behrens-Fisher type correlation-matrix"=x$corr$BF, 
  ##           "Steel type correlation-matrix"=x$corr$Steel,
  ##           "Results of the multiple Behrens-Fisher-Test"=x$test$BF[bf], 
  ##           "Results of the multiple Steel-Test"=x$test$Steel[st]
  ##           ));

  print(content);
}
