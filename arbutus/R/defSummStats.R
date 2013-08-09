## arbutus:::defSummStats

## function to build list of default summary statistics
## used internally in traitStats

## takes no arguments

defSummStats <- function(){
	
	list("REML.sigsq"=sigsqReml, "KS.D"=ksPic, "Var.pic"=varPic, "m.pic.var"=slopePicBl, "m.pic.asr"=slopePicAsr, "m.pic.nh"=slopePicNh)
	
}