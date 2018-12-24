estimatetype2 = function(score1, cutoff){
  if(is.na(cutoff)){
    stop('The cutoff is NA')
  }else{
    return(mean(score1 <= cutoff))
  }
}
