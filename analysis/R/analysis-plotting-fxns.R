## Generate plots for analysis
require(ggplot2)

cd <- getwd()


## read in sla results
sla.fam <- read.csv(file.path(cd, "output", "results-ml-sla-family.csv"),
                header=TRUE, as.is=TRUE)
sla.ord <- read.csv(file.path(cd, "output", "results-ml-sla-order.csv"),
                header=TRUE, as.is=TRUE)
sla.time <- read.csv(file.path(cd, "output", "results-ml-sla-timeslice.csv"),
                header=TRUE, as.is=TRUE)

sla <- rbind(sla.fam, sla.ord, sla.time)

## read in seedMass results
sm.fam <- read.csv(file.path(cd, "output", "results-ml-seedMass-family.csv"),
                   header=TRUE, as.is=TRUE)
sm.ord <- read.csv(file.path(cd, "output", "results-ml-seedMass-order.csv"),
                   header=TRUE, as.is=TRUE)
sm.time <- read.csv(file.path(cd, "output", "results-ml-seedMass-timeslice.csv"),
                    header=TRUE, as.is=TRUE)

sm <- rbind(sm.fam, sm.ord, sm.time)





## little fxn to rename the rank for the time slice
clean.ml.results <- function(x){

    ## convert NAs to names so not discarded
    tmp <- which(is.na(x[,"rank"]))
    x[tmp, "rank"] <- rep("timeslice", length(tmp))

    tmp2 <- which(is.na(x[, "taxa"]))
    x[tmp2, "taxa"] <- rep("random", length(tmp2))

    ## remove all datapoints where mv.modelad
    x <- x[!is.na(x[,"mv.modelad"]), ]

    ## log mv.modelad
    x[,"mv.modelad"] <- log(x[,"mv.modelad"])

    ## make rank a factor
    x[,"rank"] <- as.factor(x[,"rank"])

    as.data.frame(x)
    
}



modelad.age.plot <- function(data){

    trait <- unique(data[,"trait"])
    
    .e <- environment()

    p <- ggplot(data, aes(log(age), mv.modelad), environment=.e)
    p <- p + geom_point(aes(colour=rank), size=3, alpha=0.9)
    p <- p + theme_bw()
    p <- p + xlab("Age of clade")
    p <- p + ylab("Log mahalanobis distance")
    p <- p + ggtitle(paste("Model adequacy versus clade age -", trait, sep=" "))
    print(p)
}




modelad.size.plot <- function(data){

    trait <- unique(data[,"trait"])
    
    .e <- environment()

    p <- ggplot(data, aes(log(size), mv.modelad), environment=.e)
    p <- p + geom_point(aes(colour=rank), size=3, alpha=0.9)
    p <- p + theme_bw()
    p <- p + xlab("Log number of taxa")
    p <- p + ylab("Log mahalanobis distance")
    p <- p + ggtitle(paste("Model adequacy versus clade size -", trait, sep=" "))
    print(p)
}



## Make plots

## Clean data
dat.sla <- clean.ml.results(sla)
dat.sm <- clean.ml.results(sm)

## model adequacy versus age
pdf(file.path(cd, "output", "results-ml-sla-adequacy-age.pdf"))
modelad.age.plot(dat.sla)
dev.off()


pdf(file.path(cd, "output", "results-ml-seedMass-adequacy-age.pdf"))
modelad.age.plot(dat.sm)
dev.off()


## model adequacy versus size
pdf(file.path(cd, "output", "results-ml-sla-adequacy-taxa.pdf"))
modelad.size.plot(dat.sla)
dev.off()


pdf(file.path(cd, "output", "results-ml-seedMass-adequacy-taxa.pdf"))
modelad.size.plot(dat.sm)
dev.off()


## combine all the datasets together
ml.res <- rbind(dat.sla, dat.sm)

modelad.size.plot.alldata <- function(data){
    
    .e <- environment()

    p <- ggplot(data, aes(log(size), mv.modelad), environment=.e)
    p <- p + geom_point(aes(colour=trait, shape=rank), size=3, alpha=0.6)
    p <- p + scale_colour_brewer(palette="Set1")
    p <- p + theme_bw()
    p <- p + xlab("Log number of taxa")
    p <- p + ylab("Log mahalanobis distance")
    p <- p + ggtitle("Model adequacy versus clade size")
    print(p)
}

pdf("output/results-ml-alldata-taxa.pdf")
modelad.size.plot.alldata(ml.res)
dev.off()













## ternacular plot


model.support.ternacular.plot <- function(res){
## Load Packages

library (cwhmisc)
library (grid)
library (scales)

# Select Data Source

# Select csv file
data = res

# Define apices
#t = top, l = left, r = right
data$l = data$aic.w.ou		
data$t = data$aic.w.bm
data$r = data$aic.w.eb


# Plot customisation

# Customise: Colours
colour.1 = "white"		# Panel background colour
colour.2 = "white"		# Plot area background colour
colour.3 = "black"		# Outer text colour
colour.4 = "black"		# Plot outline colour
colour.5 = "lightgrey"		# Gridline colour
colour.6 = "white"		# Legend box outline colour

# Customise: Font sizes
fs.r = 1.168
fs.1 = 6
fs.2 = (fs.1*2) * (fs.r^4)
fs.3 = fs.1 * (fs.r^5)

# Labels & Title
title = "Relative model support (AIC weights)"
apices = c("OU", "EB", "BM") 	# Left, Right, Top

# Other
gridlines = F		# Turn gridlines on/off with T/F






##################################################################################

# Convert data to %
data$total = data$l+data$r+data$t
data$t = (data$t/data$total)*100
data$l = (data$l/data$total)*100
data$r = (data$r/data$total)*100

# Convert abc into xy
data$x = (100-((data$t/2)+((100-data$t)*(data$l/(data$l+data$r)))))
data$y = data$t

#if t = 100 then x should be 50
n = nrow (data)
for (i in 1:n) {
    temp = data[i,]
    if (temp$y == 100) {temp$x = 50}
    data[i,] = temp
}

# Coordinates for triangle
tri = data.frame (c(0, 100, 50, 0), c(0, 0, 100, 0))
colnames(tri) = c("x", "y")

# Coordinates for grid
x1 = seq (10, 90, 10)
x2 = seq (5, 45, 5)
x3 = seq (55, 95, 5)
x4 = seq (45, 5, -5)
y1 = 0
y2 = seq (10, 90, 10)
y3 = seq (90, 10, -10)
grid = data.frame (x1, x2, x3, x4, y1, y2, y3)


# Open new quartz window
dev.new (height = 10*0.866025, width = 10)
if (gridlines == T) {gridalpha = 1} else {gridalpha = 0}




p <- ggplot (data = data, aes (x = x, y = y)) +

# Shapes, lines, text
geom_polygon (data = tri, fill = colour.2) +

geom_segment (data = grid,
aes (x = x1, y = y1, xend = x2, yend = y2),
linetype = "dashed", size = 0.5, colour = colour.5, alpha = gridalpha) +
geom_segment (data = grid,
aes (x = x1, y = y1, xend = x3, yend = y3),
linetype = "dashed", size = 0.5, colour = colour.5, alpha = gridalpha) +
geom_segment (data = grid,
aes (x = x4, y = y3, xend = x3, yend = y3),
linetype = "dashed", size = 0.5, colour = colour.5, alpha = gridalpha) +
annotate ("text", label = apices, x = c(-5, 105, 50), y = c(-5, -5, 105), size = fs.1, colour = colour.3) +
labs (title  = title) +
geom_path (data = tri, size = 0.5, colour = colour.4) +

# Appearance controls
coord_fixed (ratio = 0.866025) +
scale_x_continuous (limits = c(-20, 120), expand = c(0,0)) +
scale_y_continuous (limits = c(-20, 120), expand = c(0,0)) +
theme (

# Panel and plot attributes
plot.title = element_text (size = fs.2, colour = colour.3), 	# Plot title
plot.margin = unit (c(3, 4, 3, 3), "lines"), 					# Plot margins
plot.background = element_rect (colour = F, fill = colour.1),
panel.border = element_rect (colour = F, fill = F, size = 1), 	# Axis colours
panel.grid.major = element_blank (), 							# Remove major grid
panel.grid.minor = element_blank (), 							# Remove minor grid
panel.background = element_rect (fill = colour.1),

# Legend attributes
legend.background = element_rect (colour = colour.6, size = 0.3, fill = "white"),
legend.justification = c(1, 1),
legend.position = c(0.9, 0.9), 				# Put the legend INSIDE the plot area
legend.key = element_blank (), # switch off the rectangle around symbols in the legend
legend.box.just = "bottom",
legend.box = "horizontal",
legend.title = element_text (size = fs.3, colour = colour.3), # switch off the legend title
legend.text = element_text (size = fs.3, colour = colour.3), #sets the attributes of the legend text

# Axis attributes
axis.title.x = element_blank (),
axis.title.y = element_blank (),
axis.text.x = element_blank (),
axis.text.y = element_blank (),
axis.ticks = element_blank ()
) +

# Data layer
geom_point (aes(colour = trait, shape=rank), size = 3, alpha=0.5) +
scale_colour_brewer(palette="Set1")

print(p)
}

model.support.ternacular.plot(ml.res)
ggsave(filename="output/model-support-ml-ternacular.pdf")
dev.off()







## update names to more recent version
old2new.summ.stat.names <- function(x){
    cc <- colnames(x)
    cc[cc == "reml.sigsq"] <- "sigsq.est"
    cc[cc == "var.con"] <- "var.contrast"
    cc[cc == "slope.con.var"] <- "cor.contrast.var"
    cc[cc == "slope.con.asr"] <- "cor.contrast.asr"
    cc[cc == "slope.con.nh"] <- "cor.contrast.nh"
    cc[cc == "ks.dstat"] <- "ks.contrast"

    colnames(x) <- cc
    x
}




build.pvalue.table <- function(res, summ.stat.names){
    data <- data.frame()
    nm <- summ.stat.names
    for (i in 1:nrow(res)){
        for (j in 1:length(nm)){
            tmp <-cbind.data.frame(res[i, c("taxa", "rank", "trait", "size", "age")], nm[j], res[i, nm[j]])
            colnames(tmp) <- c(colnames(tmp[1:5]), "summ.stat", "p.value")
            data <- rbind(data, tmp)
        }
    }
    data
}
        




## histogram of p-values
pval.histogram <- function(data){

    data <- data 
    .e <- environment()
    p <- ggplot(data, aes(p.value, fill=trait), environment=.e)
    p <- p + geom_bar(binwidth=0.025)
    p <- p + scale_fill_brewer(palette="Set1")
    p <- p + theme_bw()
    p <- p + facet_wrap(~summ.stat)
    p <- p + xlab("p-value")
    p <- p + theme(axis.text = element_text(size = 7))
    p <- p + ggtitle("Distribution of p-values for summary statistic")
    print(p)
}




## plot histogram
res <- old2new.summ.stat.names(ml.res)
nm <- c("sigsq.est", "var.contrast", "cor.contrast.var", "cor.contrast.asr", "cor.contrast.nh", "ks.contrast")
dd <- build.pvalue.table(res, nm)

pval.histogram(dd)
ggsave("output/summ-stat-pvalue-dist.pdf")
dev.off()

    


    
