# ordination (PCA) plotter
plotOrd <- function (	
	obj, 
	colData,
	design = "condition",
	shapes,
	labels=F,
	cluster=NULL,
	continuous=F, 
	xlims=NULL,
	ylims=NULL,
	legend=T,
	xlabel="Dimension 1",
	ylabel="Dimension 2",
	dimx=1,
	dimy=2,
	pointSize=2,
	cbPalette=F,
	...
) {

	suppressPackageStartupMessages(require(ggplot2))
   
	if (!all(design %in% names(colData))) {
		stop("the argument 'design' should specify columns of colData")
	}
	design.df <- as.data.frame(colData[, design,drop = FALSE])
	group <- if (length(design) > 1) {
		factor(apply(design.df, 1, paste, collapse = " : "))
	}
	else {
		if(continuous) {
			colData[[design]]
		} else {	
			as.factor(colData[[design]])
		}
	}
	
	if (!missing(shapes)) {
		shape <- as.factor(colData[[shapes]])
		d <- data.frame(PC1 = obj[, dimx], PC2 = obj[, dimy], group = group,design.df,shape = shape)
		colnames(d)[grep("shape", colnames(d))] <- shapes
	} else {
		d <- data.frame(PC1 = obj[, dimx], PC2 = obj[, dimy], group = group,design.df)
	}

	colnames(d)[grep("group", colnames(d))] <- design

	cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

	g <- ggplot()
	g <- g + coord_fixed(ratio = 1, xlim = xlims, ylim = ylims, expand = TRUE)
	g <- g + theme_bw()
	g <- g + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	g <- g + theme(axis.line.x = element_line(size=0.5,colour = "black"),axis.line.y = element_line(size=0.5,colour = "black"),axis.text = element_text(colour = "black"))
	if(!legend) {
		g <- g+ theme(legend.position="none")
	}	
	if (!missing(shapes)) {
		g <- g + geom_point(data=d, mapping=aes(x=PC1, y=PC2, colour=group, shape=shape),size=pointSize)
		g <- g + scale_shape_discrete(name=shapes)
	} else {
		g <- g + geom_point(data=d, mapping=aes(x=PC1, y=PC2, colour=group),size=pointSize)
	}
	if(labels) {
		g <- g + geom_text(data=d, mapping=aes(x=PC1, y=PC2, label=row.names(obj),colour=group), size=(pointSize+1), vjust=2, hjust=0.5)
	}


	if(continuous) {
		#g <- g + scale_color_gradientn(colours = rainbow(5))
		g <- g + scale_color_gradient(low="#edf8b1", high="#2c7fb8",name=design)
	} else {
		if(cbPalette) {
			g<-g+scale_colour_manual(values=cbbPalette)			
		} else {	
			g <- g + scale_colour_discrete(name=design)
		}

		if(!is.null(cluster)) {
			km <- kmeans(obj,...)
			d$Cluster<-km$cluster
			g<-g+stat_ellipse(data=d,mapping=aes(x=PC1,y=PC2,fill=factor(Cluster)), geom="polygon", level=cluster, alpha=0.2)
		}
	
	}
	g <- g + xlab(xlabel)
	g <- g + ylab(ylabel)
	return(g)
}

