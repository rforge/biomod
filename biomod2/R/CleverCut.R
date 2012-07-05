CleverCut <- function(x){
  switch(EXPR=x,
         '1' = return(c(1,1)),
         '2' = return(c(1,2)),
         '3' = return(c(2,2)),
         '4' = return(c(2,2)),
         '5' = return(c(2,3)),
         '6' = return(c(2,3)),
         '7' = return(c(3,3)),
         '8' = return(c(3,3)),
         return(c(3,3)))
}