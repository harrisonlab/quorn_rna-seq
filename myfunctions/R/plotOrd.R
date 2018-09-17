# ordination (PCA) plotter
plotOrd <- function (	
	obj, 
	colData,
	design,
	shapes,
	labels=F,
	cluster=NULL,
	continuous=F,
	colourScale=c(low="#0072B2", high="#D55E00"),
	cbPalette=T,
	pointSize=2,
	textSize=11,
	textFont="Helvetica",
	xlims=NULL,
	ylims=NULL,
	legend=T,
	title=NULL,
	xlabel="Dimension 1",
	ylabel="Dimension 2",
	dimx=1,
	dimy=2,
	...
) {
	suppressPackageStartupMessages(require(ggplot2))

	d <- data.frame(PC1 = obj[, dimx], PC2 = obj[, dimy])
	
	if (!missing(design)) {
		design.df <- as.data.frame(colData[, design,drop = FALSE])
		group<- if (length(design) > 1) {
			factor(apply(design.df, 1, paste, collapse = " : "))
		}
		else {
			colData[[design]]
		}
		d <- cbind(d,design=design,design.df)
		colnames(d)[grep("group", colnames(d))] <- design
	}
	
	if (!missing(shapes)) {
		shape <- if (length(shapes) > 1) {
			factor(apply(as.data.frame(colData[, shapes,drop = FALSE]), 1, paste, collapse = " : "))
		} else {
			as.factor(colData[[shapes]])
		}
		d <- cbind(d,shapes = shape)
	}
	
	if(!is.null(cluster)) {
		km <- kmeans(obj,...)
		d$Cluster<-km$cluster
	}	

	cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

	g <- ggplot(data=d,aes(x=PC1, y=PC2))
	g <- g + coord_fixed(ratio = 1, xlim = xlims, ylim = ylims, expand = TRUE)

	g <- g + theme_classic_thin(textSize,textFont) # theme_classic_thin is in the metabarcoding functions
	# g <- g + theme_classic(textSize,textFont)

	if(!legend) {
		g <- g + theme_classic(textSize,textFont) %+replace% theme(legend.position="none")
	} else {
		g <- g + theme_classic(textSize,textFont)
	}

	if(!missing(design)) {
		g <- g + aes(colour=group)
	}
	
	if (!missing(shapes)) {
		g <- g + aes(shape=shapes)
		g <- g + scale_shape_discrete(name=shapes)
	}
	
	g <- g + geom_point(size=pointSize)
	
	if(labels) {
		g <- g + aes(label=row.names(obj))
		g <- g + geom_text(size=(pointSize+1), vjust=2, hjust=0.5)
	}
	
	if(continuous) {
		g <- g + scale_color_gradient(low=colourScale[1], high=colourScale[2],name=design)
	} else {
		if(cbPalette) {
			g<-g+scale_colour_manual(values=cbbPalette)			
		} else {
			g<-g+scale_colour_manual(colours=rainbow())
		}
	}
	
	if(!is.null(cluster)) {
		g<-g+stat_ellipse(geom="polygon", level=cluster, alpha=0.2)
	}
	
	g <- g + xlab(xlabel)
	g <- g + ylab(ylabel)
	
	return(g)
}
