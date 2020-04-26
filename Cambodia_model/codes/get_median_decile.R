decile_index = floor(quantile(seq(1:nsim),c(seq(0.1,1,by=0.1))))
names = c()
for(i in 1:10)
{
  names = c(names,paste0("x",i))
}

get_median_decile = function(matrices,decile_index)
{
  x = apply(matrices, 2, function(x) max(x, na.rm = TRUE))
  y = cbind.data.frame(seq(1:length(x)),x)
  y = y[order(y[,2]),]
  z = y[decile_index,1]
  v=cbind.data.frame(seq(1:500),matrices[,z])
  colnames(v) = c("time",names)
  v =melt(v, id="time")
  return(v)
  
}

return_index_median_decile = function(matrices,decile_index)
{
  x = apply(matrices, 2, function(x) max(x, na.rm = TRUE))
  y = cbind.data.frame(seq(1:length(x)),x)
  y = y[order(y[,2]),]
  z = y[decile_index,1]
 
  return(z)
  
}