## arbutus:::plot.picBl

## fxn for plotting results for summary stat 4: the slope of the contrasts versus their variances


## NOTE: not sure if I want to include these plotting fxns in the package or not


plot.picBl <- function(unit.tree, data=NULL, colours=c("dodgerblue4", "darkblue")){
	
	## make sure the unit.tree is of class 'unit.tree'
	## if not create a unit.tree using supplied data
	if (!inherits(unit.tree, "unit.tree"))
		unit.tree <- as.unit.tree(unit.tree, data)

	## sd of the pics
	sd.pic <- sqrt(unit.tree$pics[,"variance"])
	
	## absolute value of the pics
	abs.pic <- abs(unit.tree$pics[,"contrasts"])
	
	pp <- cbind.data.frame(abs.pic, sd.pic)
	colnames(pp) <- c("contrasts", "variance") 

	p <- ggplot(pp, aes(variance, contrasts), environment=environment()) + stat_smooth(method="lm", fill=colours[1], colour=colours[2], size=1, alpha=0.1) + geom_point() + xlab("Var[Contrast]") + ylab("Contrast") + ggtitle("Variation relative to branch length")
	
	print(p)
	
}