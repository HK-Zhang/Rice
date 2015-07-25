"report" <-
function(msg=NULL, style=0, char="-")
{
  if (is.null(msg)) msg <- "";
  
  if (is.vector(msg))
    msg <- unlist(msg)
  else
    stop("msg must be of type vector");
  
  char <- substr(char, 1, 1);

  underlined <- function (arg)
  {
    c(arg, paste(rep(char, max(nchar(msg))), collapse=""));
  }

  border <- function(arg) 
  {
    n <- length(msg);
    ml <- max(nchar(msg));
    space <- paste(rep(" ", ml), collapse="");
    line <- paste(rep(char, ml+4), collapse="");
    msg <- paste(msg, substr(rep(space, n), rep(1, n), ml-nchar(msg)), sep=""); 
    c(line, paste(char, msg, char), line);          
  }

  sfun <- list(underlined = underlined,
               border = border
               );
  
  if (is.numeric(style) && length(style)==1 && any(style==1:length(sfun)))
    msg <- sfun[[style]](msg)
  else if (is.character(style) && length(style)==1 && !is.null(sfun[[style]]))
    msg <- sfun[[style]](msg)
  
  m <- matrix(msg, ncol=1);
  colnames(m) <- "";
  rownames(m) <- rep("", length(msg));
  print.noquote(m); 
}
