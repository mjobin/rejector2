# Functions for processing rejector output
# Load by putting this in your R working directory and
# typing source("rejector.r")

# Written by Matthew Jobin, Department of Anthropology, Stanford University


rejplotsens <- function(iname, xcols, ycols)
{

qr<- read.delim(iname, header=TRUE)

name <- unlist(strsplit(iname, split = "\\."))
#oname <- paste(name[1], ".tiff", sep ="")
#pdf(oname, title=oname)


	for(x in xcols){
		for(y in ycols){
		oname <- paste(name[1], names(qr)[x], "-" , names(qr)[y], ".tiff", sep ="")
		tiff(oname)
		
			 plot(qr[,x], qr[,y], main = paste("From file", iname, "\nSensitivity of", names(qr)[y], " to changes in ", names(qr)[x]), xlab = names(qr)[x], ylab = names(qr)[y])
		dev.off()
		#end y loop
		}
	#end x loop
	}

}

rejfiltertolerance <- function(iname, test, tolerance, target)
{

# This function reads in a rejector output file, then makes an output file based on closeness within tolerance of the given 
# target. this duplictaes the function of some of rejecdtor itself, but is useful to re-set the tolerance
# if all the data is saved by using the command line -p switch.
# NOTE: This can be done to only one statistic or test per use of the function.
# The given target should be read by running the input file with the -rs switch.

qr<- read.delim(iname, header=TRUE)

name <- unlist(strsplit(iname, split = "\\."))
oname <- paste(name[1], "-",names(qr)[test], "-",tolerance,".txt", sep ="")

upper <- target + (target*tolerance)
cat("upper " , upper, "\n")
lower <- target - (target*tolerance)
cat("lower " , lower, "\n")

qs <- subset(qr, qr[test]>lower & qr[test]<upper)
#qs <- subset(qr, qr[test]<upper)

cat("\nwhee\n")
qs
#cat(unlist(qs), file=oname)
#print(qs, file=oname)

write.table(qs, file=oname, sep="\t", append=FALSE)

}

rejoned <- function(iname, tests, outcol, rlimit, priortype="noprior", v1=0, v2=1, truemark = NULL, density = "FALSE")
{

# This function reads in rejector output, then takes each of the T/F columns specified in
# tests and takes only those rows marked T. It then makes histograms of those columns in the
# subset specified in outcol
# usage:
#   use rejoned(
# 	iname, - input name of the rejector outfile (must be processed by gluedat.pl first)
#   tests, - vector of column numbers of the T/F output columns you would like to examine
#   outcol, - vector the column numbers of population parameter values you would like to process
#   rlimit, - the total number of iterations in the run of rejector
#   priortype - the type of prior distribution  to superimpose over each histogram. Leaving blank will print no prior distribution
#   v1 - The minimum of a uniform distn, the mean of a gaussian distribution
#   v2 - The maximum of a uniform distribution, or the standard deviation of a gaussian distribution
#   truemark - A line drawn at the true parameter value, if you know what it is (usually for testing)
#   density - TRUE if you prefer a density plot to a histogram
#	)
#

qr<- read.delim(iname, header=TRUE)

name <- unlist(strsplit(iname, split = "\\."))
oname <- paste(name[1], ".pdf", sep ="")






# Calculate limits of of the plot if a prior distribution is specified 
if(priortype == "gaussian" | priortype == "normal") {
	#v1 then becomes the mean, and v2 the s.d.
	xmin = 0
	xmax = 2000
	#xmin <- (v1-(3*v2));
	#xmax <- (v1+(3*v2));
}

else if (priortype == "uniform") {
	xmin = v1
	xmax = v2

}


pdf(oname, title=oname)




for(h in tests){
		if(length (tests) == 1) {
			if(tests == 0){
				qs <- subset(qr, select=c(outcol))
				
				testname = ""
				
			}
			
	else {
	 	testname = c( "for ", names(qr[h]))
		qs <- subset(qr, qr[,h] == TRUE, select=c(outcol))
	
	}
			
			
		}
	else {
	 	testname = c( "for ", names(qr[h]))
		qs <- subset(qr, qr[,h] == TRUE, select=c(outcol))
	
	}
	
	for(i in 1:NCOL(qs)){
				if(NROW(qs[,i]) > 0){

					
					if(priortype == "noprior") pxlim = range(qs[,i])
					else pxlim = c(xmin, xmax)
					#pxlim = range(qs[,i])
		
					if(density == "TRUE"){
						plot(density(qs[,i],from=rmin, to=rmax), main = paste("Density plot of", names(qs)[i], "for", names(qr)[h], "\nwith", NROW(qs), "runs accepted from a run limit of", rlimit))
					}
					else{
					#hist(qs[,i], freq=FALSE, breaks=seq(0,10, by=0.05), xlim = pxlim, xlab = names(qs)[i])

#hist(qs[,i], freq=FALSE, breaks=seq(0,40000, by=1000), xlim = pxlim, main = paste("Histogram of", names(qs)[i], "for", names(qr)[h], "\nwith", NROW(qs), "runs accepted from a run limit of", rlimit), xlab = names(qs)[i])
					hist(qs[,i], freq=FALSE, xlim = pxlim, main = paste("Histogram of", names(qs)[i], "for", names(qr)[h], "\nwith", NROW(qs), "runs accepted from a run limit of", rlimit), xlab = names(qs)[i])
		 			
		 			}
		 			# breaks=seq(-1,1, by=0.000025) breaks=40
		 						 

				 if(is.null(truemark) == FALSE){
				 	truepoint = c(truemark, 1)
				

					 points(truemark, 0, type = "h",  col='blue', lwd=15);
				 
				 }
		 		
				rejprior1d(priortype, v1, v2)
		 		
		 		
		 		}
			}
		}
	


dev.off()
}


###############################################################

rejonedb <- function(stats, outcol, rlimit)
{
	flist <- list.files(path=".", pattern="-out.txt")
	

	

	for (i in 1:length(flist)) {
		for(j in stats){
			for(k in outcol){
			rejoned(flist[i], stats, outcol, rlimit)
			}
		}
	}


}

###############################################################

rejonednewaxes <- function(iname, tests, outcol, rlimit, priortype="noprior", v1=0, v2=1, truemark = NULL, density = "FALSE")
{

# This function reads in rejector output, then takes each of the T/F columns specified in
# tests and takes only those rows marked T. It then makes histograms of those columns in the
# subset specified in outcol
# usage:
#   use rejoned(
# 	iname, - input name of the rejector outfile (must be processed by gluedat.pl first)
#   tests, - vector of column numbers of the T/F output columns you would like to examine
#   outcol, - vector the column numbers of population parameter values you would like to process
#   rlimit, - the total number of iterations in the run of rejector
#   priortype - the type of prior distribution  to superimpose over each histogram. Leaving blank will print no prior distribution
#   v1 - The minimum of a uniform distn, the mean of a gaussian distribution
#   v2 - The maximum of a uniform distribution, or the standard deviation of a gaussian distribution
#   truemark - A line drawn at the true parameter value, if you know what it is (usually for testing)
#   density - TRUE if you prefer a density plot to a histogram
#	)
#

qr<- read.delim(iname, header=TRUE)

name <- unlist(strsplit(iname, split = "\\."))
oname <- paste(name[1], ".pdf", sep ="")






# Calculate limits of of the plot if a prior distribution is specified 
if(priortype == "gaussian" | priortype == "normal") {
	#v1 then becomes the mean, and v2 the s.d.
	xmin = 0
	xmax = 2000
	#xmin <- (v1-(3*v2));
	#xmax <- (v1+(3*v2));
}

else if (priortype == "uniform") {
	xmin = v1
	xmax = v2

}


pdf(oname, title=oname)




for(h in tests){
		if(length (tests) == 1) {
			if(tests == 0){
				qs <- subset(qr, select=c(outcol))
				
				testname = ""
				
			}
			
	else {
	 	testname = c( "for ", names(qr[h]))
		qs <- subset(qr, qr[,h] == TRUE, select=c(outcol))
	
	}
			
			
		}
	else {
	 	testname = c( "for ", names(qr[h]))
		qs <- subset(qr, qr[,h] == TRUE, select=c(outcol))
	
	}
	
	for(i in 1:NCOL(qs)){
				if(NROW(qs[,i]) > 0){

					
					if(priortype == "noprior") pxlim = range(qs[,i])
					else pxlim = c(xmin, xmax)
					#pxlim = range(qs[,i])
		
					if(density == "TRUE"){
						plot(density(qs[,i],from=rmin, to=rmax), main = paste("Density plot of", names(qs)[i], "for", names(qr)[h], "\nwith", NROW(qs), "runs accepted from a run limit of", rlimit))
					}
					else{
					#hist(qs[,i], freq=FALSE, breaks=seq(0,10, by=0.05), xlim = pxlim, xlab = names(qs)[i])

#hist(qs[,i], freq=FALSE, breaks=seq(0,40000, by=500), xlim = pxlim, main = paste(""), ylab = "Number of Accepted Simulations", xlab="Population Size (Number of Chromosomes)")
#					hist(qs[,i], freq=FALSE, xlim = pxlim, main = paste(""), ylab = "Number of Accepted Simulations", xlab="Bottleneck Severity (X-fold Reduction)")
					hist(qs[,i], freq=FALSE, xlim = pxlim, main = paste(""), ylab = "Number of Accepted Simulations", xlab="Bottleneck Time (# Generations Before Present)")
		 			
		 			}
		 			# breaks=seq(-1,1, by=0.000025) breaks=40
		 						 

				 if(is.null(truemark) == FALSE){
				 	truepoint = c(truemark, 1)
				

					 points(truemark, 0, type = "h",  col='blue', lwd=15);
				 
				 }
		 		
				rejprior1d(priortype, v1, v2)
		 		
		 		
		 		}
			}
		}
	


dev.off()
}




###############################################################
rejprior1d <- function(priortype, v1, v2)
# This function draws a specified prior distribution. Designed to be called by
# rejoned to superimpose priors over posterior distributions
# usage:
#   use rejprior(
#   priortype - the type of prior distribution  to superimpose over each histogram. Leaving blank will print no prior distribution
#   v1 - The minimum of a uniform distn, the mean of a gaussian distribution
#   v2 - The maximum of a uniform distribution, or the standard deviation of a gaussian distribution
#   )
#

{

			#Draw prior
		 		if(priortype == "gaussian" | priortype == "normal") {
				#v1 then becomes the mean, and v2 the s.d.
				xmin <- (v1-(1*v2));
				xmax <- (v1+(1*v2));
				curve(dnorm(x, v1, v2), add=TRUE, col='red', lwd=3);
			}
			else if (priortype == "gamma") {
				#v1 becomes scale, v2 shape parameter
				curve(dgamma(x, shape=v2, scale=v1), add=TRUE, col='red', lwd=3);
				}
			
				
			else if (priortype == "uniform") {
				curve(dunif(x, v1, v2), add=TRUE, col='red', lwd=3);
				}

}





###############################################################
rejprior2d <- function(priortype1, v11, v12, num1, priortype2, v21, v22, num2, col1min, col1max, col2min, col2max)
# This function draws a specified 2d density plot from two prior distributions. Designed to be called by
# rejtwod to superimpose priors over posterior distributions
# usage:
#   use rejprior(
#   priortype - the type of prior distribution  to superimpose over each histogram. Leaving blank will print no prior distribution
#   v1 - The minimum of a uniform distn, the mean of a gaussian distribution
#   v2 - The maximum of a uniform distribution, or the standard deviation of a gaussian distribution
#   )
#

{

	dist1 <- seq(col1min, col1max, length=1000)
	dist2 <- seq(col2min, col2max, length=1000)

		#Obtain prior distn using num1 number of entries for the first prior

			#Draw prior
		 		if(priortype1 == "gaussian" | priortype1 == "normal") {
				#v1 then becomes the mean, and v2 the s.d.
				#contour(dist1, dist2, outer(dnorm(dist1, v11, v12),dnorm(dist2,v21, v22)), col="grey82", drawlabels = FALSE)
				pdist1 = dnorm(dist1,v11,v12)
			}
			else if (priortype1 == "gamma") {
				#v1 becomes scale, v2 shape parameter
				
				#contour(dist1, dist2, outer(dgamma(dist1, shape = v12, scale = v11),dgamma(dist2, shape = v22, scale = v21)), col="grey82", drawlabels = FALSE)
				pdist1 = dgamma(dist1, shape=v12, scale=v11)
			}
			
				
			else if (priortype1 == "uniform") {
				pdist1 = dunif(dist1, v11, v12)
				#contour(dist1, dist2, outer(dunif(dist1, v11, v12),dunif(dist2,v21, v22)), col="grey82", drawlabels = FALSE)
			}



			#Draw prior
		 		if(priortype2 == "gaussian" | priortype2 == "normal") {
				#v1 then becomes the mean, and v2 the s.d.
					pdist2 = dnorm(dist2,v21,v22)
			}
			else if (priortype2 == "gamma") {
				#v1 becomes scale, v2 shape parameter
				

				pdist2 = dgamma(dist2, shape=v22, scale=v21)
			}
			
				
			else if (priortype2 == "uniform") {
			
						pdist2 = dunif(dist2, v21, v22)
			}

		

	contour(dist1, dist2, outer(pdist1, pdist2), col="grey82", drawlabels = FALSE)
}












###############################################################

rejtwod <- function(iname, tests, col1, col1min, col1max, col2,  col2min, col2max, rlimit, priortype1="noprior", v11=0, v12=1, priortype2="noprior", v21=0, v22=1, truemark1=NULL, truemark2=NULL) 
{

# This function reads in rejector output, then takes each of the T/F columns specified in
# tests and takes only those rows marked T. It then makes histograms of those columns in the
# subset specified in outcol
# usage:
#   use rejtwod(
# 	iname, - input name of the rejector outfile (must be processed by gluedat.pl first)
#   tests, - vector of column numbers of the T/F output columns you would like to examine
#   col1, -  column number of the first set of population parameter values you would like to process (e.g. DemeSize.1, HEvent0-Time)
#   col1min, col1max - graphing boundaries for col1
#   col2, -  column number of the second set of  population parameter values you would like to process
#   col2min, col2max - graphing boundaries for col2
#   rlimit, - the total number of iterations in the run of rejector
#   priortype1, priortype2 - the type of prior distributions  to superimpose over each histogram. Leaving blank will print no prior distribution
#   v11, v21 - The minimum of a uniform distn, the mean of a gaussian distribution, the scale of a gamma distribution
#   v12, v22 - The maximum of a uniform distribution, or the standard deviation of a gaussian distribution, the shape of a gamma distribution
#   truemark1, truemark2 - The two known parameters values (if you know them) for col1 and col2
#	)
#

library(MASS)

qr<- read.delim(iname, header=TRUE)

name <- unlist(strsplit(iname, split = "\\."))
oname <- paste(name[1], ".pdf", sep ="")



	pdf(oname, title=oname)
	for(h in tests){
		qs <- subset(qr, qr[,h] == TRUE, select=c(col1, col2))
		if(NROW(qs[,1]) > 1){
			if(NROW(qs[,2]) > 1){
				# Prior plot
				if(priortype1 != "noprior" && priortype2 != "noprior") {
					rejprior2d (priortype1, v11, v12, length(qs[,1]), priortype2, v21, v22, length(qs[,2]), col1min, col1max, col2min, col2max)
					#Posterior Plot
					contour(kde2d(qs[,1], qs[,2], lims=c(col1min, col1max, col2min, col2max)), add=TRUE, drawlabels = FALSE, main = paste("Figure XX: Density plot of", names(qs)[1], "and\n", names(qs)[2], "for", names(qr)[h], "\nwith", NROW(qs), "runs accepted from a run limit of", rlimit), xlab = names(qs)[1], ylab = names(qs)[2])
				
					}
					
					else {
					
						contour(kde2d(qs[,1], qs[,2], lims=c(col1min, col1max, col2min, col2max)), drawlabels = FALSE, main = paste("Figure XX: Density plot of", names(qs)[1], "and\n", names(qs)[2], "for", names(qr)[h], "\nwith", NROW(qs), "runs accepted from a run limit of", rlimit), xlab = names(qs)[1], ylab = names(qs)[2])
					
					}
				
				
				#Posterior Plot
				contour(kde2d(qs[,1], qs[,2], lims=c(col1min, col1max, col2min, col2max)), add=TRUE, drawlabels = FALSE, main = paste("Figure XX: Density plot of", names(qs)[1], "and\n", names(qs)[2], "for", names(qr)[h], "\nwith", NROW(qs), "runs accepted from a run limit of", rlimit), xlab = names(qs)[1], ylab = names(qs)[2])
				

				
				if(is.null(truemark1) == FALSE && is.null(truemark2) == FALSE){

				 	points(truemark1, truemark2,  col='black', lwd=5);

	
				 
				 }
				
				
				
			}	
		}
		
		



	}
	dev.off()


}






###############################################################


rejtwodb <- function(tests, outcol, rlimit, rmin, rmax) 
{
	flist <- list.files(path=".", pattern="-out.txt")

	for (i in 1:length(flist)) {
		cat(flist[i],"\n")
		rejdens(flist[i],tests, outcol, rlimit, rmin, rmax);
	
	}


}


###############################################################

rejtwodsuper <- function(iname, tests, col1, col1min, col1max, col2,  col2min, col2max, rlimit, truemark1=NULL, truemark2=NULL) 
{

# This function reads in rejector output, then takes each of the T/F columns specified in
# tests and takes only those rows marked T. It then makes histograms of those columns in the
# subset specified in outcol
# usage:
#   use rejtwod(
# 	iname, - input name of the rejector outfile (must be processed by gluedat.pl first)
#   tests, - vector of column numbers of the T/F output columns you would like to examine
#   col1, -  column number of the first set of population parameter values you would like to process (e.g. DemeSize.1, HEvent0-Time)
#   col1min, col1max - graphing boundaries for col1
#   col2, -  column number of the second set of  population parameter values you would like to process
#   col2min, col2max - graphing boundaries for col2
#   rlimit, - the total number of iterations in the run of rejector
#   truemark1, truemark2 - The two known parameters values (if you know them) for col1 and col2
#	)
#

library(MASS)

qr<- read.delim(iname, header=TRUE)

name <- unlist(strsplit(iname, split = "\\."))
oname <- paste(name[1], ".pdf", sep ="")



	pdf(oname, title=oname)
	for(h in tests){
		qs <- subset(qr, qr[,h] == TRUE, select=c(col1, col2))
		
		qt <- subset(qr, TRUE, select=c(col1, col2))
		
		contour(kde2d(qt[,1], qt[,2], lims=c(col1min, col1max, col2min, col2max)),  col="grey82", drawlabels = FALSE, main = paste("Figure XX: Density plot of", names(qs)[1], "and\n", names(qs)[2], "for", names(qr)[h], "\nwith", NROW(qs), "runs accepted from a run limit of", rlimit), xlab = names(qs)[1], ylab = names(qs)[2])
				
		
		
		if(NROW(qs[,1]) > 1){
			if(NROW(qs[,2]) > 1){
				
				
				
				#Posterior Plot
				contour(kde2d(qs[,1], qs[,2], lims=c(col1min, col1max, col2min, col2max)), add=TRUE, drawlabels = FALSE, main = paste("Figure XX: Density plot of", names(qs)[1], "and\n", names(qs)[2], "for", names(qr)[h], "\nwith", NROW(qs), "runs accepted from a run limit of", rlimit), xlab = names(qs)[1], ylab = names(qs)[2])
				

				
				if(is.null(truemark1) == FALSE && is.null(truemark2) == FALSE){

				 	points(truemark1, truemark2,  col='black', lwd=5);

	
				 
				 }
				
				
				
			}	
		}
		
		



	}
	dev.off()


}


###############################################################


rejtwodsuperb <- function(tests, col1, col1min, col1max, col2,  col2min, col2max, rlimit, truemark1=NULL, truemark2=NULL) 
{
	flist <- list.files(path=".", pattern="-out.txt")

	for (i in 1:length(flist)) {
		cat(flist[i], tests, col1, col1min, col1max, col2,  col2min, col2max, rlimit, truemark1=NULL, truemark2=NULL,"\n")
		rejtwodsuper(flist[i], tests, col1, col1min, col1max, col2,  col2min, col2max, rlimit, truemark1=NULL, truemark2=NULL) 
	
	}


}


###############################################################

rejtwodstatssuper <- function(iname,  col1,  col2, tol=0.05) 
{
	library(MASS)
	qr<- read.delim(iname, header=TRUE)
	
	
	cat("Rejector Stats for Output File: ", iname, "\n")
	#cat("Column: ", names(qr)[col1], " real value is: " , qr[1,col1], "\n")
	#cat("Column: ", names(qr)[col2], " real value is: " , qr[1,col2], "\n")
	
	
	c1 <- qr[1,col1]
	c1l <- c1 - (c1*tol)
	c1h <- c1 + (c1*tol)
	
	
	c2 <- qr[1,col2]
	c2l <- c2 - (c2*tol)
	c2h <- c2 + (c2*tol)
	
	cat("Column: ", names(qr)[col1], " real value is: " , c1, " | " , c1l, "<>" , c1 , "<>" , c1h, "\n")
	cat("Column: ", names(qr)[col2], " real value is: " , c2, " | " , c2l, "<>", c2 , "<>" , c2h, "\n")
	
	
	
	qs <- subset(qr, qr[,col1] != "-999" & qr[,col2] != "-999", select=c(col1, col2))
	
	
	cat("\n\nNumber non-missing: ", NROW(qs), "\n")
	
	min1 <- min(qs[,1])
	max1 <- max(qs[,1])
	min2 <- min(qs[,2])
	max2 <- max(qs[,2])
	
	
	cat(names(qs)[1], " Min: ", min1, " Max: ",max1, "\n")
	cat(names(qs)[2], " Min: ", min2, " Max: ",max2, "\n")
	#print(qs)
	
	qt <- subset(qs, qs[,1] <= c1h & qs[,1] >= c1l & qs[,2] <= c2h & qs[,2] >= c2l )
	
	#print (qt)
	
	cat("\n\nNumber accepted within tolerance of " , tol , ": ", NROW(qt), "\n")
	cat(names(qt)[1], " Min: ", min(qt[,1]), " Max: ",max(qt[,1]), "\n")
	cat(names(qt)[2], " Min: ", min(qt[,2]), " Max: ",max(qt[,2]), "\n")
	

	
	name <- unlist(strsplit(iname, split = "\\."))
	oname <- paste(name[1], ".pdf", sep ="")
	pdf(oname, title=oname)
	
	contour(kde2d(qs[,1], qs[,2], lims=c(min1, max1, min2, max2)), nlevels=20, col="grey82", drawlabels = FALSE, main = paste("Density plot of file ", iname , " columns ", names(qs)[1], "and\n", names(qs)[2]), xlab = names(qs)[1], ylab = names(qs)[2])
	
	cat("\nPrinting within tolerance graph.\n")
	
	contour(kde2d(qt[,1], qt[,2], lims=c(min1, max1, min2, max2)), add=TRUE, drawlabels = FALSE)
				
	points(c1, c2,  col='blue', lwd=5);
		
	dev.off()
}



###############################################################

rejthreed <- function(iname, tests, col1, col2,   col3,  rlimit) 
{

# This function reads in rejector output, then takes each of the T/F columns specified in
# tests and takes only those rows marked T. It then makes 3d scatterplots of the data.
# usage:
#   use rejthreed(
# 	iname, - input name of the rejector outfile (must be processed by gluedat.pl first)
#   tests, - vector of column numbers of the T/F output columns you would like to examine
#   col1, -  column number of the first set of population parameter values you would like to process (e.g. DemeSize.1, HEvent0-Time)
#   col2, -  column number of the second set of  population parameter values you would like to process
#   col3, -  column number of the third set of  population parameter values you would like to process
#   rlimit, - the total number of iterations in the run of rejector
#	)
#

library(misc3d)
library(rgl)

qr<- read.delim(iname, header=TRUE)

name <- unlist(strsplit(iname, split = "\\."))





	for(h in tests){
		qs <- subset(qr, qr[,h] == TRUE, select=c(col1, col2, col3))
		if(NROW(qs[,1]) > 1){
			if(NROW(qs[,2]) > 1){
			if(NROW(qs[,3]) > 1){
					

				open3d()
				
				#Posterior Plot
				#d <- kde3d(qs[,1], qs[,2],qs[,3], n = 40, lims=c(col1min, col1max, col2min, col2max, col3min, col3max))
				
				#contour3d(d$d, exp(-12), d$x, d$y, d$z, color = "blue", color2 = "gray", engine = "standard")
				#contour3d(d$d, exp(-12), d$x, d$y, d$z, color = "blue", color2 = "gray")
				
				plot3d(qs[,1], qs[,2], qs[,3], names(qs)[1], names(qs)[2], names(qs)[3])
				#open3d()
				#rgl.surface(qs[,1], qs[,2], qs[,3])
				
				
				}
			}	
		}
		
		



	}



}




###############################################################


rejthreedb <- function(tests, col1,  col2,  col3,  rlimit) 
{
	flist <- list.files(path=".", pattern="-out.txt")

	for (i in 1:length(flist)) {
		cat(flist[i], tests, col1, col2,  col3, rlimit,"\n")
		rejthreed(flist[i], tests, col1, col2,   col3,   rlimit) 
	
	}


}


###############################################################
rejsum <- function (iname,  tests, outcol, oname="NONE")
{
	# This function reads rejector output, and prints the mean and variance of the T/F columns
	# selected as in rejdens above to standard output.

qr<- read.delim(iname, header=TRUE)


	name <- unlist(strsplit(iname, split = "\\."))
	sname <- paste(name[1], "-summary.txt", sep ="")

if(oname != "NONE") sname <- oname




#Print Header of Summary File
if(oname == "NONE") {
	cat("Rejector Summary of Output File: ", iname, "\n", file=sname)
	cat("Total entiries in output file ", NROW(qr), "\n", file=sname)
	cat("Printed: ", date(), "\n\n", file=sname, append=TRUE)
	cat("COL", " ", "NAME", " ", "ACCEPTS", " ",  "MEAN", " ", "VAR", " ", "MEDIAN", " " , "2.5%Int" , " " , "97.5%Int" , "\n", file=sname, append=TRUE)

}




	for(h in tests){
		qs <- subset(qr, qr[,h] == TRUE, select=c(outcol))
		
		for(i in 1:NCOL(qs)){


			if(NROW(qs[,i]) > 1){
			
				# get de mode
				
				#modetbl <- table(qs[,])
				#themode <- as.numeric(names(modetbl)[modetbl==max(modetbl)])
				
				if(oname != "NONE") cat (iname, " ", file=sname, append=TRUE)

				cat(names(qs)[i], " ", names(qr)[h], " ", NROW(qs[,i]), " ", mean(qs[,i]), " ", var(qs[,i]), quantile(qs[,i], probs = c(0.5, 0.025, 0.975)),  "\n", file=sname, append=TRUE)
				

		 	}
		}
	}

cat("\n", file=sname, append=TRUE)
}

###############################################################

rejsumb <- function(oname, tests, outcol)
{
	flist <- list.files(path=".", pattern="-out.txt")
	
	cat(date(), "\n\n", file=oname)
		cat("FILE", " " , "COL", " ", "NAME", " ", "ACCEPTS", " ",  "MEAN", " ", "VAR", " ", "MEDIAN", " " , "2.5%Int" , " " , "97.5%Int" , "\n", file=oname, append=TRUE)


	for (i in 1:length(flist)) {
		cat(flist[i],"\n")
		#cat(flist[i],"\n", file=oname, append=TRUE)
		rejsum(flist[i],  tests, outcol, oname)

	
	}


}






###############################################################
rejnames <- function(iname)

{
	
	qr<- read.delim(iname, header=TRUE)
	
	names(qr)
		


}



###############################################################
rejde <- function(iname)

{
	
	qr<- read.delim(iname, header=TRUE)
	
	de(qr)
		


}


rejhilo <- function(loname, origname, hiname)
{
	lo<- read.delim(loname, header=TRUE)
	orig<- read.delim(origname, header=TRUE)
	hi <- read.delim(hiname, header=TRUE)
	
	name <- unlist(strsplit(origname, split = "\\."))
	oname <- paste(name[1], "hiloplots.pdf", sep ="")
	
	
	for(i in 1:NCOL(lo)){
		hist(lo[,i], freq=FALSE)
		hist(orig[,i], freq=FALSE)
		hist(hi[,i], freq=FALSE)

	
	

	}

}



###############################################################

rejprob <- function(iname, tests, col1, col1min, col1max, col2,  col2min, col2max, rlimit, priortype1="noprior", v11=0, v12=1, priortype2="noprior", v21=0, v22=1, truemark1=NULL, truemark2=NULL) 
{

# This function reads in rejector output, then takes each of the T/F columns specified in
# tests and takes only those rows marked T. It then makes histograms of those columns in the
# subset specified in outcol
# usage:
#   use rejtwod(
# 	iname, - input name of the rejector outfile (must be processed by gluedat.pl first)
#   tests, - vector of column numbers of the T/F output columns you would like to examine
#   col1, -  column number of the first set of population parameter values you would like to process (e.g. DemeSize.1, HEvent0-Time)
#   col1min, col1max - graphing boundaries for col1
#   col2, -  column number of the second set of  population parameter values you would like to process
#   col2min, col2max - graphing boundaries for col2
#   rlimit, - the total number of iterations in the run of rejector
#   priortype1, priortype2 - the type of prior distributions  to superimpose over each histogram. Leaving blank will print no prior distribution
#   v11, v21 - The minimum of a uniform distn, the mean of a gaussian distribution, the scale of a gamma distribution
#   v12, v22 - The maximum of a uniform distribution, or the standard deviation of a gaussian distribution, the shape of a gamma distribution
#   truemark1, truemark2 - The two known parameters values (if you know them) for col1 and col2
#	)
#

library(MASS)

qr<- read.delim(iname, header=TRUE)

name <- unlist(strsplit(iname, split = "\\."))
oname <- paste(name[1], ".pdf", sep ="")



	pdf(oname, title=oname)
	for(h in tests){
		qs <- subset(qr, select=c(col1, col2))
		if(NROW(qs[,1]) > 1){
			if(NROW(qs[,2]) > 1){

			#cat(qs[,1], " ", qs[,2])

					#Posterior Plot
					#contour(kde2d(qs[,1], qs[,2], lims=c(col1min, col1max, col2min, col2max)), drawlabels = FALSE)
				
					

				
				
			}	
		}
		
		



	}
	dev.off()


}

###############################################################

rejxstatsb <- function(oname, stats, outcol)
{
	flist <- list.files(path=".", pattern="-out.txt")
	
	cat(date(), "\n\n", file=oname, append=TRUE)
	#cat("**************************************\n", file=oname, append=TRUE)

	for (i in 1:length(flist)) {
		cat(flist[i],"\n")
		#cat(flist[i],"\n", file=oname, append=TRUE)
		
		qr<- read.delim(flist[i], header=TRUE)
		
		
		
		for(h in outcol){
			cat(mean(qr[,h]),"\n", file=oname, append=TRUE)
			qt <- c(qt, mean(qr[,h]))
			#cat(qt, "\n")
			qt
		#stats loop ends
		}
		

		#cat("**************************************\n", file=oname, append=TRUE)
	
	}


}

