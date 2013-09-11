## arbutus:::edgeBM

## This fxn rescales a edge according to a BM  process

## Takes 3  arguments:
## br.start -- branching time of start of branch
## br.end -- branching time of end of branch
## sigsq -- rate of variance increase

## Note that br.end must be greater than or equal to  br.start

edgeBM <- function(br.start, br.end, sigsq){

  if (br.end < br.start)
    stop("start of branch must proceed end of branch")
  
  bl <- (br.end - br.start) * sigsq
  
  bl
  
}

