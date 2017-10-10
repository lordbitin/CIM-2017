helper.function <- function()
{
  return(1)
}

###
# 
###
helper.extractOIDFromCompleteOIDString <- function (oid)
{
  return(mongo.oid.from.string(str_sub(oid,10,-2)))
}

##
# 
##
helper.getUserName <- function (clientIP,createdAt)
{
  return(
    paste(str_sub(createdAt,-2),str_sub(clientIP,-3),sep = "_")
  )
}

helper.dataFrameDiff <- function(x.1,x.2,...)
{
  x.1p <- do.call("paste", x.1)
  x.2p <- do.call("paste", x.2)
  x.1[! x.1p %in% x.2p, ]
}

helper.splitIntoUnequalChunks  <- function(vector,chunkSizes) {
  
  if (length(vector)!=sum(chunkSizes))
    stop("Vector length must be equal to the sum of chunk sizes")
  
  split(vector,rep(1:length(chunkSizes),chunkSizes))
}

#It is not always possible to chunk in equal lengths!!!!!
helper.splitIntoEqualChunks <- function(vector,numberOfChunks) {
  split(vector, (seq(length(vector))-1) %/% (length(vector)/numberOfChunks))
}

helper.slice<-function(x,n) {
  N<-length(x);
  lapply(seq(1,N,n),function(i) x[i:min(i+n-1,N)])
}

helper.goodUnlist <- function (x) {
  v <- as.data.frame(t(unlist(x)),stringsAsFactors = FALSE)
  vc <- t(apply(v,1,as.character))
  vn <- t(apply(v,1,as.numeric))
  nas <- which(is.na(vn))
  v[-nas] <- vn[-nas]
  v[nas] <- vc[nas]
  v
}

helper.combine <- function(..., prefix = "", sep = "_") {
  paste0(prefix, levels(interaction(..., sep = sep)))
}

