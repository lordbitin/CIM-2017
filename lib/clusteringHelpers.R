clusteringHelpers.medoids = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}