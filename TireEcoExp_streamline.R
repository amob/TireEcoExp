###
###Script goals
###analyze tire global change ecology experiment data for effects on:
		## duckweed fitness, 
		## microbe ''fitness'', 
		## fitness correlations, 
		## and duckweed trait "greenness"
###NOTE to users: paths to inputs and outputs will not be the same on your local system. 

library(MCMCglmm)
library(MuMIn)

################################
################################
##DEFINING FUNCTIONS
###############################
###############################

bufferX <- function(x,p) {  #range, but with added buffer (for figure axes)
	r<- range(x,na.rm=T)
	add <- c(-1,1)*p*(r[2]-r[1])
	return(r+add)
	}
	
range01=function(x){  #convert a variable to a range between 0 and 1
newnums=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
return(newnums)
}

std.error <- function(x){sd(x)/sqrt(length(x))} 

#function for analyzing data in map and dat files
MapToWellsT <- function(dat,map,firstcol, sumcols, meancols){ #columns labeled "plate" "round" "roi", "image",   roi is specific to dat only.
	n <- nrow(map)
	outdata <- matrix(NA,nrow=n, ncol=(length(sumcols)+length(meancols)+1))
	for(i in 1:n){
		rois <- map[i,firstcol:ncol(map)] #what happens to NA values
		p <- map$plate[i]
		r <- map$round[i]
		welldat <- dat[dat$plate==p & dat$round==r & dat$roi%in%rois,]
		welldat.sums <- colSums(welldat[,sumcols])
		welldat.means <- colMeans(welldat[,meancols],na.rm=T)
		wellstats <- c(nrow(welldat),welldat.sums,welldat.means)
		outdata[i,] <- wellstats
		}
	mappeddata <- cbind(map[,1:(firstcol-1)],outdata)
	colnames(mappeddata) <- c(colnames(map)[1:(firstcol-1)],"particles",colnames(dat)[sumcols],colnames(dat)[meancols])
	return(mappeddata)
}

#This function converts OD to cells/ul to a desired number of cells in inocula, according to the formula at http://www.reric.org/wordpress/archives/169
# Build a function to calculate cells/ul for each sample
#we used this function to standardize inocula concentrations, as in O'Brien et al 2020 Microbial Ecology & O'Brien et al 2019 American Journal of Botany
ODtoCells <- function(x){
	log10cellspml <- 9.9 - 4/(1 + (x/.14)^0.8 ) #log base 10 cellspermL
	cellpul <-  (10^log10cellspml)/1000 #per ul
	return(cellpul)
} #note that the minimum estimate is 794 cells / ul... (relevant for expressing OD as cell count, if of interest, below. not in O'Brien et al.)


################################
################################
##reading in raw data (i.e. from image J), processing, reading in treatments, combining into analysis dataframes
###############################
###############################

#read in temperature measures and inferred points then impute weighted mean temperature for each well
ImputedT <- read.csv("~/Dropbox/TireEcoExp/raw-input-data/TireEcoExpPlateTemps_inferred_imputed_s.csv",header=T,stringsAsFactors=T)
#THIS is Raw input file for analysis 2
cornercols <- which(colnames(ImputedT)%in%c("a1.skirt","h1.skirt","a12.skirt","h12.skirt"))
rndpltavs <- lapply(1:3, function(r) lapply(1:3, function(p) colMeans(ImputedT[ImputedT$round==r & ImputedT$plate==p,cornercols])  ) )
rndpltavsmat <- matrix(unlist(rndpltavs),ncol=4,byrow=T)

#make FULL temp vector for trts file. i.e. own temp for each well based on a linear interpolation of the 4 corners
#plates are 12 columns, and 8 rows, trts file reports full columns first
platetemps <- function(cornertemps,rows=8,cols=12){  ##
	horizontals <- matrix(nrow=rows,ncol=cols)
	verticals <- matrix(nrow=rows,ncol=cols)
	agreement <- matrix(nrow=rows,ncol=cols)
	horizontals[1,] <- seq(from=cornertemps[1], to = cornertemps[3],length.out=cols)
	horizontals[8,] <- seq(from=cornertemps[2], to = cornertemps[4],length.out=cols)
	verticals[,1] <- seq(from=cornertemps[1], to = cornertemps[2],length.out=rows)
	verticals[,12] <- seq(from=cornertemps[3], to = cornertemps[4],length.out=rows)
	for(i in 2:7){
		horizontals[i,] <- seq(from=verticals[i,1],to=verticals[i,12],length.out=cols)
	}
	for(j in 2:11){
		verticals[,j] <- seq(from=horizontals[1,j],to=horizontals[8,j],length.out=rows)
	}
	agreement <- (horizontals + verticals)/2
return(agreement)
}
simplates <- list(lapply(1:3, function(z) platetemps(rndpltavsmat[z,])),lapply(4:6, function(z) platetemps(rndpltavsmat[z,])),lapply(7:9, function(z) platetemps(rndpltavsmat[z,]) ))


#
#read in start and end image analyses and roi to well map files
tirestartdat <- read.csv("~/Dropbox/TireEcoExp/raw-input-data/Tire_EcoExp_StartDat_s.csv",header=T,stringsAsFactors=F)#
#THIS is Raw data 3
tireenddat <- read.csv("~/Dropbox/TireEcoExp/raw-input-data/Tire_EcoExp_EndDat_s.csv",header=T,stringsAsFactors=F)#
#THIS is Raw data 4
tirestartmap <- read.csv("~/Dropbox/TireEcoExp/raw-input-data/StartFrondsandMap_s.csv",header=T,stringsAsFactors=F)
#THIS is Raw data 1
tireendmap <- read.csv("~/Dropbox/TireEcoExp/raw-input-data/EndFrondsandMap_s.csv",header=T,stringsAsFactors=F)
#THIS is Raw data 2

trts <- read.csv("~/Dropbox/TireEcoExp/raw-input-data/TireEcoExp_Design_s.csv")
#THIS is Raw input file for analysis 1
trts$numrow <- as.numeric(as.factor(trts$row))
#trts$co2 <- seq(from = 1000, to = 400, length.out=8)[trts$numrow]
trts$temp <- sapply(1:nrow(trts), function(z) simplates[[trts$round[z]]][[trts$plate[z]]][trts$numrow[z],trts$col[z]] )

#optical density files
r1od <- read.csv("~/Dropbox/TireEcoExp/raw-input-data/AO TireEcoRound1_sort_s.csv",header=T,stringsAsFactors=F)#
#THIS is Raw data 5
r2od <- read.csv("~/Dropbox/TireEcoExp/raw-input-data/AO TireEcoRound2_sort_s.csv",header=T,stringsAsFactors=F)#
#THIS is Raw data 6
r3od <- read.csv("~/Dropbox/TireEcoExp/raw-input-data/AO TireEcoRound3_sort_s.csv",header=T,stringsAsFactors=F)#
#THIS is Raw data 7


#should be true if sheets organized as expected.
 colnames(r1od)[c(1:6,11,12)] ==  colnames(r2od)[c(1:6,11,12)]
 colnames(r1od)[c(1:6,11,12)] ==  colnames(r3od)[c(1:6,15,16)]

oddat <- rbind(r1od[,c(1:6,11,12)],r2od[,c(1:6,11,12)],r3od[,c(1:6,15,16)])

#preprocess dat files
startplate <- rep(NA,length.out=nrow(tirestartdat))
startround <- rep(NA,length.out=nrow(tirestartdat))
startplate[grep("Plate1",tirestartdat$image)] <- 1
startplate[grep("Plate2",tirestartdat$image)] <- 2
startplate[grep("Plate3",tirestartdat$image)] <- 3
startround[grep("Round1",tirestartdat$image)] <- 1
startround[grep("Round2",tirestartdat$image)] <- 2
startround[grep("Round3",tirestartdat$image)] <- 3
endplate <- rep(NA,length.out=nrow(tireenddat))
endround <- rep(NA,length.out=nrow(tireenddat))
endplate[grep("Plate1",tireenddat$image)] <- 1
endplate[grep("Plate2",tireenddat$image)] <- 2
endplate[grep("Plate3",tireenddat$image)] <- 3
endround[grep("Round1",tireenddat$image)] <- 1
endround[grep("Round2",tireenddat$image)] <- 2
endround[grep("Round3",tireenddat$image)] <- 3
#
plateround <- c("Plate1Round1","Plate2Round1","Plate3Round1","Plate1Round2","Plate2Round2","Plate3Round2","Plate1Round3","Plate2Round3","Plate3Round3")
plateroundp <- c(1,2,3,1,2,3,1,2,3)
plateroundr <- c(1,1,1,2,2,2,3,3,3)
pasteplateround <- paste(plateroundp,plateroundr,sep=".")
endimages <- sapply(plateround, function(z) tireenddat$image[grep(z,tireenddat$image)[1]])
imagename <- sapply(1:nrow(tireendmap), function(z) endimages[pasteplateround == paste(tireendmap$plate,tireendmap$round,sep=".")[z] ])
tireendmap2 <- data.frame(cbind(imagename,tireendmap))
colnames(tireendmap2)[1] <- "image"
tireenddat2<-data.frame(cbind(endplate,endround,tireenddat))
startimages <- sapply(plateround, function(z) tirestartdat$image[grep(z,tirestartdat$image)[1]])
imagenameS <- sapply(1:nrow(tirestartmap), function(z) startimages[pasteplateround == paste(tirestartmap$plate,tirestartmap$round,sep=".")[z] ])
tirestartmap2 <- data.frame(cbind(imagenameS,tirestartmap))
colnames(tirestartmap2)[1] <- "image"
tirestartdat2<-data.frame(cbind(startplate,startround,tirestartdat))
#color data extraction
colors <- c("(blue)","(green)","(red)")
colorrows <- lapply(colors,function(z) grep(z,tireenddat$image))
colorrowsS <- lapply(colors,function(z) grep(z,tirestartdat$image))
#
tireenddat3 <- tireenddat2[-unlist(colorrows),]
tireenddat3$redraw <- tireenddat2[colorrows[[3]],]$mean
tireenddat3$greenraw <- tireenddat2[colorrows[[2]],]$mean
tireenddat3$blueraw <- tireenddat2[colorrows[[1]],]$mean
##CHECK THAT THEY LINE UP, should be 1/4 nrow(tireenddat2)
sum(tireenddat3$area == tireenddat2[colorrows[[3]],]$area)
sum(tireenddat3$area == tireenddat2[colorrows[[2]],]$area)
sum(tireenddat3$area == tireenddat2[colorrows[[1]],]$area)
tireenddat3$perred <- tireenddat3$redraw/(3*tireenddat3$mean)
tireenddat3$pergreen <- tireenddat3$greenraw/(3*tireenddat3$mean)
tireenddat3$perblue <- tireenddat3$blueraw/(3*tireenddat3$mean)
#
tirestartdat3 <- tirestartdat2[-unlist(colorrowsS),]
tirestartdat3$redraw <- tirestartdat2[colorrowsS[[3]],]$mean
tirestartdat3$greenraw <- tirestartdat2[colorrowsS[[2]],]$mean
tirestartdat3$blueraw <- tirestartdat2[colorrowsS[[1]],]$mean
##CHECK THAT THEY LINE UP, should be 1/4 nrow(tirestartdat2)
sum(tirestartdat3$area == tirestartdat2[colorrowsS[[3]],]$area)
sum(tirestartdat3$area == tirestartdat2[colorrowsS[[2]],]$area)
sum(tirestartdat3$area == tirestartdat2[colorrowsS[[1]],]$area)
tirestartdat3$perred <- tirestartdat3$redraw/(3*tirestartdat3$mean)
tirestartdat3$pergreen <- tirestartdat3$greenraw/(3*tirestartdat3$mean)
tirestartdat3$perblue <- tirestartdat3$blueraw/(3*tirestartdat3$mean)


#run maptowells, generate data on a per-well basis
colnames(tireenddat3)[1:2] <- c("plate","round") #to fit requirements
colnames(tirestartdat3)[1:2] <- c("plate","round") #to fit requirements
tiremapdatend <- MapToWellsT(tireenddat3,tireendmap2,firstcol=9,sumcols=c(5,7),meancols=c(6,12,13,18,19:21,25:27))
#adding variables
tiremapdatend$numrow <- as.numeric(as.factor(tiremapdatend$row))
tiremapdatend$microbe <- trts$microbe
co2s <- as.factor(trts$co2)
co2levs <-  c(400,round(400 + (600/7)*c(1:6)),1000)
levels(co2s) <- co2levs
tiremapdatend$co2 <- as.numeric(as.character(co2s))
tiremapdatstart <- MapToWellsT(tirestartdat3,tirestartmap2,firstcol=8,sumcols=c(5,7),meancols=c(6,12,13,18,19:21,25:27))
tiremapdatstart$numrow <- as.numeric(as.factor(tiremapdatstart$row))
tiremapdatend$temp <- trts$temp
deltapix <- tiremapdatend$area-tiremapdatstart$area # change in pixel area variable
tiremapdatend$deltapix <- deltapix
tiremapdatend$deltaratio <- (deltapix/(tiremapdatend$area+tiremapdatstart$area))
tiremapdatend$tire <- 5*(tiremapdatend$plate - 1)
tirex<- as.factor(trts$tire)
levels(tirex) <- c("0","0.25","0.5")
tiremapdatend$tirex <- as.numeric(as.character(tirex))
tiremapdatend$is.alive <- as.numeric(tiremapdatend$area > 0)
tiremapdatend$is.dead <- 1-tiremapdatend$is.alive
#just double check rows and sorting
sum(paste(oddat$round,oddat$plate,tolower(oddat$row),oddat$col) == paste(tiremapdatend$round, tiremapdatend$plate, tiremapdatend$row, tiremapdatend$col))
dim(tiremapdatend)
tiremapdatend$od600b <- oddat$blankod600
tiremapdatend$flod600b <- log(ifelse(oddat$blankod600 <= 0, min(oddat$blankod600[oddat$blankod600>0]), oddat$blankod600)) # in full dataset this is 59/60 0s for 600/750, but in inoc only dataset it is only 5/7
#store files
write.csv(tiremapdatend,"~/Dropbox/TireEcoExp/tiremapdatend.csv",row.names=F)
write.csv(tiremapdatstart,"~/Dropbox/TireEcoExp/tiremapdatstart.csv",row.names=F)

#ranges and means of temperatures
range(tiremapdatend$temp[tiremapdatend$tire==0])
range(tiremapdatend$temp[tiremapdatend$tire==5])
range(tiremapdatend$temp[tiremapdatend$tire==10])
mean(tiremapdatend$temp[tiremapdatend$tire==0])
mean(tiremapdatend$temp[tiremapdatend$tire==5])
mean(tiremapdatend$temp[tiremapdatend$tire==10])
mean(tiremapdatend$temp[tiremapdatend$tire!=0])

#subset data to different analyzable fractions
nondry <- tiremapdatend[tiremapdatend$dry!="y",]  # this measure of dryness is not fully reliable, so we exclude also wells based on dead duckweed below
nondry$co210 <- nondry$co2 / 100
inoconly <- nondry[nondry$microbe=="y",]
nondrylive <- nondry[nondry$is.alive==1,]
inoconlylive <- inoconly[inoconly$is.alive==1,] #is also excluding dry wells

tempNDL <- nondrylive[nondrylive$temp>18.5 & nondrylive$temp < 27,] #within half a degree
tempNDLinoc <- inoconlylive[inoconlylive$temp>18.5 & inoconlylive$temp < 27,] #within half a degree

mean(tempNDL$temp[tempNDL$tire==0]) #means are only closer together.
mean(tempNDL$temp[tempNDL$tire==5])
mean(tempNDL$temp[tempNDL$tire==10])
mean(tempNDL$temp[tempNDL$tire!=0])


################################
################################
##deltapix plots and models
###############################
###############################

#from most complex to simple; dried wells eliminated. dead wells also eliminated, because identifying a dry well (some may have died at lower but not dry by sticking on side.) not always clear, and dead ones have a threshold on how many pixels they can loose.
##each time, the highest order n.s. term is thrown out, if multiple n.s. terms are of same order, the largest pMCMC term is removed. n.s. terms that are subsets of significant higher order terms are retained.  if the resulting DIC is higher than the more complex model, we retain the more complex model. (doesn't change much)

#reverse stepwise selection
summary(MCMCglmm(deltapix~temp + co210+ microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210+ temp:microbe + co210:microbe + tire:temp:microbe + tire:temp:co210+ temp:co210:microbe + tire:microbe:co210 + temp:co210:microbe:tire ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=500,thin=50))#rm 4way
summary(MCMCglmm(deltapix~temp + co210+ microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210+ temp:microbe + co210:microbe + tire:temp:microbe + tire:temp:co210+ temp:co210:microbe + tire:microbe:co210 ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=500,thin=50))#m 3way
summary(MCMCglmm(deltapix~temp + co210+ microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210+ temp:microbe + co210:microbe + tire:temp:microbe + temp:co210:microbe + tire:microbe:co210 ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=500,thin=50))#rm 3way
summary(MCMCglmm(deltapix~temp + co210+ microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210+ temp:microbe + co210:microbe + tire:temp:microbe + tire:microbe:co210 ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=500,thin=50))#rm 3way
summary(MCMCglmm(deltapix~temp + co210+ microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210+ temp:microbe + co210:microbe + tire:temp:microbe   ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=5000,thin=50)) #rm co210:microbe
summary(MCMCglmm(deltapix~temp + co210+ microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210+ temp:microbe + tire:temp:microbe   ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=5000,thin=50))# temp:co210seemingly not as high pval as temp:tire, but temp:tire required due to sig 3-way
summary(MCMCglmm(deltapix~temp + co210+ microbe + tire + tire:microbe + temp:tire + co210:tire + temp:microbe + tire:temp:microbe   ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=5000,thin=50)) #rm co210:tire. because it is only interaction term left that is both n.s. and not part of higher order sig term.
summary(MCMCglmm(deltapix~temp + co210+ microbe + tire + tire:microbe + temp:tire + temp:microbe + tire:temp:microbe   ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=5000,thin=50)) #dic gain on previous is only slight; but co210 is out
summary(MCMCglmm(deltapix~temp  + microbe + tire + tire:microbe + temp:tire + temp:microbe + tire:temp:microbe   ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=5000,thin=50)) #DONE
#fit best model with more iterations
best.dpx.ndl.LONG <- (MCMCglmm(deltapix~temp + microbe + tire + tire:microbe + temp:tire + temp:microbe + tire:temp:microbe   ,random = ~round, data=nondrylive ,verbose=F,nitt=1000000,burnin=5000,thin=100,pr=T))
summary(best.dpx.ndl.LONG)

#removing extreme temps and refitting
summary(MCMCglmm(deltapix~temp + co210+ microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210+ temp:microbe + co210:microbe + tire:temp:microbe + tire:temp:co210+ temp:co210:microbe + tire:microbe:co210 + temp:co210:microbe:tire ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#rm 4way
summary(MCMCglmm(deltapix~temp + co210+ microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210+ temp:microbe + co210:microbe + tire:temp:microbe + tire:temp:co210+ temp:co210:microbe + tire:microbe:co210 ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))# rm temp microbe co2
summary(MCMCglmm(deltapix~temp + co210+ microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210+ temp:microbe + co210:microbe + tire:temp:microbe + tire:temp:co210 + tire:microbe:co210 ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#rm unclear between tire microbe co2 and tire temp co2; fit both
summary(MCMCglmm(deltapix~temp + co210+ microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210+ temp:microbe + co210:microbe + tire:temp:microbe + + tire:microbe:co210 ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#rm the other, co2 microbe tire
summary(MCMCglmm(deltapix~temp + co210+ microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210+ temp:microbe + co210:microbe + tire:temp:microbe + tire:temp:co210 ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#rm the other, tire temp co2
summary(MCMCglmm(deltapix~temp + co210+ microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210+ temp:microbe + co210:microbe + tire:temp:microbe  ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#rm co2:micr
summary(MCMCglmm(deltapix~temp + co210+ microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210+ temp:microbe + tire:temp:microbe  ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#rm temp:co2
summary(MCMCglmm(deltapix~temp + co210+ microbe + tire + tire:microbe + temp:tire + co210:tire + temp:microbe + tire:temp:microbe  ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#rm co2:tire
summary(MCMCglmm(deltapix~temp + co210+ microbe + tire + tire:microbe + temp:tire + temp:microbe + tire:temp:microbe  ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#rm co2 (gain on previous line's model in DIC was very slight)
summary(MCMCglmm(deltapix~temp + microbe + tire + tire:microbe + temp:tire + temp:microbe + tire:temp:microbe  ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#rm  ()
#same result; actually more significant
best.dpx.ndl.LONGSUB <- (MCMCglmm(deltapix~temp + microbe + tire + tire:microbe + temp:tire + temp:microbe + tire:temp:microbe   ,random = ~round, data=tempNDL ,verbose=F,nitt=1000000,burnin=5000,thin=100,pr=T))
summary(best.dpx.ndl.LONGSUB)

#means, standard errors, and HPDIs for ploting
dpx.m.mns <- sapply(1:3, function(z) tapply(nondrylive$deltapix[nondrylive$plate==z],nondrylive$microbe[nondrylive$plate==z],mean))
dpx.m.ses <- sapply(1:3, function(z) tapply(nondrylive$deltapix[nondrylive$plate==z],nondrylive$microbe[nondrylive$plate==z],std.error))
t.s <- seq(from = min(nondrylive$temp), to = max(nondrylive$temp),length.out=1000)
mod1post <- best.dpx.ndl.LONG$Sol
pred.m.mn <- lapply(c(0,5,10), function(tire) sapply(t.s, function(z)                mean(mod1post[,1] + mod1post[,2]*z + mod1post[,3]*1 + mod1post[,4]*tire + mod1post[,5]*1*tire + mod1post[,6]*z*tire + mod1post[,7]*z*1 + mod1post[,8]*z*tire*1 + rowMeans(mod1post[,9:11])) ) )
pred.m.CI <- lapply(c(0,5,10), function(tire) sapply(t.s, function(z) HPDinterval(as.mcmc(mod1post[,1] + mod1post[,2]*z + mod1post[,3]*1 + mod1post[,4]*tire + mod1post[,5]*1*tire + mod1post[,6]*z*tire + mod1post[,7]*z*1 + mod1post[,8]*z*tire*1 + rowMeans(mod1post[,7:9])),0.9 )) )
pred.n.mn <- lapply(c(0,5,10), function(tire) sapply(t.s, function(z)                mean(mod1post[,1] + mod1post[,2]*z + mod1post[,3]*0 + mod1post[,4]*tire + mod1post[,5]*0*tire + mod1post[,6]*z*tire + mod1post[,7]*z*0 + mod1post[,8]*z*tire*0 + rowMeans(mod1post[,7:9])) ) )
pred.n.CI <- lapply(c(0,5,10), function(tire) sapply(t.s, function(z) HPDinterval(as.mcmc(mod1post[,1] + mod1post[,2]*z + mod1post[,3]*0 + mod1post[,4]*tire + mod1post[,5]*0*tire + mod1post[,6]*z*tire + mod1post[,7]*z*0 + mod1post[,8]*z*tire*0 + rowMeans(mod1post[,7:9])),0.9 )) )


pdf("~/Dropbox/TireEcoExp/results_deltapix_TxTxM_plusTrX.pdf",height=2,width=5)
par(mfrow=c(1,3))
par(oma = c(0,3,2.5,0))
par(mar = c(4,2,0,1))
plot(as.vector(dpx.m.mns)~ rep(c(0,5,10),each=2), pch=16,xlim=c(-2,12), cex=1.5,
		ylim=c(-100,3000),ylab="",xlab="",xaxt="n", col=rgb(0,c(0,0.5),0))
	axis(side=1,at=c(0,5,10),labels=c("0x","0.25x","0.5x") )
 	arrows(x0=rep(c(0,5,10),each=2),y0=as.vector(dpx.m.mns - dpx.m.ses),
 		y1=as.vector(dpx.m.mns +dpx.m.ses),length=0,lwd=2, col=rgb(0,c(0,0.5),0) )
 	mtext("Change in pixel area",side=2,line=3)
	mtext("Tire leachate concentration",side=1,line=2.75)
	mtext("a)",side=3,line=1,adj=-0.25)
	abline(h=0,lty=3,col=rgb(0,0,0,alpha=0.75))
plot(deltapix~temp,data=nondrylive[nondrylive$plate==1 , ],ylim=c(-6000,10000),xlim=c(14,29),ylab="",xlab="",
		pch=1,cex=0.75, col=rgb(0,0.5*(as.numeric(nondrylive$microbe[nondrylive$plate==1])-1),0,alpha=0.5))
	text(14, 9000,labels="0x",adj=0 )
	polygon(c(t.s,rev(t.s)), c(pred.m.CI[[1]][1,], rev(pred.m.CI[[1]][2,])),border=NA, col= rgb(0,0.5,0,alpha=0.25) )
	lines(pred.m.mn[[1]]~t.s,col=rgb(0,0.5,0))
	polygon(c(t.s,rev(t.s)), c(pred.n.CI[[1]][1,], rev(pred.n.CI[[1]][2,])),border=NA, col= rgb(0,0,0,alpha=0.25) )
	lines(pred.n.mn[[1]]~t.s,col=rgb(0,0,0))
	abline(h=0,lty=3,col=rgb(0,0,0,alpha=0.75))
		mtext("b)",side=3,line=1,adj=-0.25)
plot(deltapix~temp,data=nondrylive[nondrylive$plate==3 , ],ylim=c(-6000,10000),xlim=c(14,29),ylab="",xlab="",
		pch=1,cex=0.75,col=rgb(0,0.5*(as.numeric(nondrylive$microbe[nondrylive$plate==3])-1),0,alpha=0.5))
	text(14, 9000,labels="0.5x",adj=0 )
	polygon(c(t.s[200:900],rev(t.s[200:900])), c(pred.m.CI[[3]][1,200:900], rev(pred.m.CI[[3]][2,200:900])),border=NA, col= rgb(0,0.5,0,alpha=0.25) )
	lines(pred.m.mn[[3]][200:900]~t.s[200:900],col=rgb(0,0.5,0))
	polygon(c(t.s[200:900],rev(t.s[200:900])), c(pred.n.CI[[3]][1,200:900], rev(pred.n.CI[[3]][2,200:900])),border=NA, col= rgb(0,0,0,alpha=0.25) )
	lines(pred.n.mn[[3]][200:900]~t.s[200:900],col=rgb(0,0,0))
 	mtext(expression("Temperature"~degree*C),side=1,line=2.75,adj=35)
	abline(h=0,lty=3,col=rgb(0,0,0,alpha=0.75))
	mtext("c)",side=3,line=1,adj=-0.25)
	legend(13,y=-1500,c("with microbiome","without"),fill=c(rgb(0,0.5,0),rgb(0,0,0)),bty="n",horiz=F)
dev.off()

pdf("~/Dropbox/TireEcoExp/results_deltapix_nomicr.pdf",height=2,width=5)
par(mfrow=c(1,3))
par(oma = c(0,3,2.5,0))
par(mar = c(4,2,0,1))
plot(as.vector(dpx.m.mns[1,])~ c(0,5,10), pch=16,xlim=c(-2,12), cex=1.5,
		ylim=c(-100,3000),ylab="",xlab="",xaxt="n")
	axis(side=1,at=c(0,5,10),labels=c("0x","0.25x","0.5x") )
 	arrows(x0=rep(c(0,5,10)),y0=as.vector(dpx.m.mns[1,] - dpx.m.ses[1,]),
 		y1=as.vector(dpx.m.mns[1,] + dpx.m.ses[1,]),length=0,lwd=2)
 	mtext("Change in pixel area",side=2,line=3)
	mtext("Tire leachate concentration",side=1,line=2.75)
	mtext("a)",side=3,line=1,adj=-0.25)
	abline(h=0,lty=3,col=rgb(0,0,0,alpha=0.75))
plot(deltapix~temp,data=nondrylive[nondrylive$plate==1 & nondrylive$microbe=="n"	, ],ylim=c(-6000,10000),xlim=c(14,29),ylab="",xlab="",
		pch=1,cex=0.75, col=rgb(0,0,0,alpha=0.5))
	text(14, 9000,labels="0x",adj=0 )
	polygon(c(t.s,rev(t.s)), c(pred.n.CI[[1]][1,], rev(pred.n.CI[[1]][2,])),border=NA, col= rgb(0,0,0,alpha=0.25) )
	lines(pred.n.mn[[1]]~t.s,col=rgb(0,0,0))
	abline(h=0,lty=3,col=rgb(0,0,0,alpha=0.75))
dev.off()


################################ 
################################
##OD plots and models
###############################
###############################

### OD is  not normal. using logged values, see dataframe building

#inoculation effects
inoc.mod <- (MCMCglmm(flod600b~microbe, random = ~round, data = nondrylive, verbose =F, nitt = 1000000,burnin=5000,thin=100))
summary(inoc.mod)
flod.inoc.mn <- tapply(nondrylive$flod600b, nondrylive$microbe,mean)
exp(flod.inoc.mn)
flod.inoc.se <- tapply(nondrylive$flod600b, nondrylive$microbe,std.error)
flod.inoc.intvs <- rbind(exp(flod.inoc.mn + flod.inoc.se),  exp(flod.inoc.mn - flod.inoc.se))
#
pdf("~/Dropbox/TireEcoExp/inoc_od_ndl.pdf",width=4,height=4)
plot(exp(flod.inoc.mn)~c(1,2),pch=16,ylab="Optical Density 600 nm",xlab="",ylim=bufferX(exp(flod.inoc.mn),1.01), xlim = c(0.5,2.5),xaxt="n",cex=2)
axis(side=1, at=c(1,2), labels = c("uninoculated","inoculated"))
arrows(x0 = c(1,2) , y0 = flod.inoc.intvs[1,], y1=flod.inoc.intvs[2,],length=0,lwd=2)
dev.off()
#             post.mean  l-95% CI  u-95% CI eff.samp  pMCMC   
# (Intercept)  0.019263 -0.028220  0.070962      900  0.264   
# microbey     0.008071  0.002663  0.013167      900 <0.001 **
# ---


#reverse stepwise regression
summary(MCMCglmm((flod600b)~temp + tire + co210+ tire:temp + tire:co210+ temp:co210+ temp:tire:co210,random=~round ,data=inoconlylive,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T)) #last term came out at 0.105 pval, so ran this several times.
summary(MCMCglmm((flod600b)~temp + tire + co210+ tire:temp + tire:co210+ temp:co210,random=~round ,data=inoconlylive,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))# rm temp:co2; note this model often fits slightly worse than previous, but simpler models fit better in DIC by several units
summary(MCMCglmm((flod600b)~temp + tire + co210+ tire:temp + tire:co210 ,random=~round ,data=inoconlylive,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T)) #rm temp:tire
summary(MCMCglmm((flod600b)~temp + tire + co210+ tire:co210 ,random=~round ,data=inoconlylive,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T)) #rm  tire:co2
summary(MCMCglmm((flod600b)~temp + tire + co210,random=~round ,data=inoconlylive,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))#rm co2
summary(MCMCglmm((flod600b)~temp + tire ,random=~round ,data=inoconlylive,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))# stop.
#DIC lowest at last one.
#fitting longer for reporting parameters
mod1.ODsimpleLONG <- (MCMCglmm((flod600b)~temp + tire ,random=~round ,data=inoconlylive,verbose=F,nitt=1000000,burnin=5000,thin=50,pr=T))# 
summary(mod1.ODsimpleLONG)


summary(MCMCglmm((flod600b)~temp + tire + co210+ tire:temp + tire:co210+ temp:co210+ temp:tire:co210,random=~round ,data=tempNDLinoc,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T)) #rm 3-way
summary(MCMCglmm((flod600b)~temp + tire + co210+ tire:temp + tire:co210+ temp:co210,random=~round ,data=tempNDLinoc,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T)) #rm tmp:co2
summary(MCMCglmm((flod600b)~temp + tire + co210+ tire:temp + tire:co210,random=~round ,data=tempNDLinoc,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T)) #rm temp tire or tire co2, close
summary(MCMCglmm((flod600b)~temp + tire + co210+ tire:temp ,random=~round ,data=tempNDLinoc,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T)) #rm temp tire
summary(MCMCglmm((flod600b)~temp + tire + co210 + tire:co210,random=~round ,data=tempNDLinoc,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T)) #rm tire co2
summary(MCMCglmm((flod600b)~temp + tire + co210 ,random=~round ,data=tempNDLinoc,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T)) #rm co2
summary(MCMCglmm((flod600b)~temp + tire  ,random=~round ,data=tempNDLinoc,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T)) #stop
#again no change in best model.
mod1.ODsimpleLONGSUB <- (MCMCglmm((flod600b)~temp + tire ,random=~round ,data=tempNDLinoc,verbose=F,nitt=1000000,burnin=5000,thin=50,pr=T))# 
summary(mod1.ODsimpleLONGSUB)



#generating means, predictions, HPDI, etc for plotting, reporting
odpost <- mod1.ODsimpleLONG$Sol
t.s <- seq(from = min(inoconly$temp), to = max(inoconly$temp),length.out=1000)
pred.od.mn <- lapply(c(0,5,10), function(tire) sapply(t.s, function(z)                mean(odpost[,1] + odpost[,2]*z + odpost[,3]*tire + rowMeans(odpost[,4:6])) ) )
pred.od.CI <- lapply(c(0,5,10), function(tire) sapply(t.s, function(z) HPDinterval(as.mcmc(odpost[,1] + odpost[,2]*z + odpost[,3]*tire + rowMeans(odpost[,4:6])),0.95 )) )
tmean <- mean(t.s)
pred.od.tr.mn <- sapply(c(0,5,10), function(tire)                mean(odpost[,1] + odpost[,2]*tmean + odpost[,3]*tire + rowMeans(odpost[,4:6]) ) )
pred.od.tr.CI <- sapply(c(0,5,10), function(tire)                HPDinterval(odpost[,1] + odpost[,2]*tmean + odpost[,3]*tire + rowMeans(odpost[,4:6])) ) 
tr.s <- c(0,5,10)
exp(pred.od.tr.CI)
#             [,1]        [,2]       [,3]
# [1,] 0.004355133 0.007642494 0.01177040
# [2,] 0.006419395 0.010010944 0.01792584 
exp(pred.od.tr.mn)
#0.005305631 0.008772787 0.014505683
exp(pred.od.CI[[2]][,c(1,length(t.s))])
#             [,1]       [,2]
# [1,] 0.001947476 0.01819768
# [2,] 0.004652139 0.03674911
exp(pred.od.mn[[2]][c(1,length(t.s))])
#[1] 0.003022383 0.025463942
flod.tr.mns <- tapply(inoconlylive$flod600b,inoconlylive$tire,mean)
flod.tr.ses <- tapply(inoconlylive$flod600b,inoconlylive$tire,std.error)
dist.trmean <- range01(abs(inoconlylive$tire-5))
dist.trL <- range01(inoconlylive$tire-5)

pdf("~/Dropbox/TireEcoExp/flod_effs_soloX.pdf",width=6,height=3)
par(mar=c(0,2,0,1))
par(oma=c(4,2,2,0))
par(mfrow=c(1,2))
plot(flod600b~temp,data=inoconly,pch=16,cex=0.8-0.3*dist.trmean,ylab="",xlab="", xlim=c(14,27.75), #,cex=1.25-dist.trmean
		ylim=c(-8,0),col=rgb(0,0,0,alpha=0.8*(1-dist.trmean/1.5)))#, col=rgb(0,0,0, alpha=(2-dist.trmean)/2 ))#,#col=rgb(1,0,0,alpha=0.5))
		polygon(c(t.s,rev(t.s)), c(pred.od.CI[[2]][1,], rev(pred.od.CI[[2]][2,])  ),col=rgb(0,0,0,alpha=0.25))
		lines(pred.od.mn[[2]]~t.s,lwd=2,col=rgb(0,0,0)) 
		mtext("Log optical density",side=2,line=2.5)
		mtext("Temperature",side=1,line=2.5)
		mtext("a)",side=3,line=0.5,adj=-0.25)
plot(flod.tr.mns~c(0,5,10), pch=16,cex=2,ylab="",xlab="",xlim=c(-2,12),ylim=c(-6,-4),xaxt="n")
	axis(side=1,at=c(0,5,10),labels=c("0x","0.25x","0.5x"))
	arrows(c(0,5,10),y0=flod.tr.mns - flod.tr.ses, y1=flod.tr.mns +flod.tr.ses,length=0,lwd=2)
	mtext("Tire leachate concentration",side=1,line=2.5)
	mtext("b)",side=3,line=0.75,adj=-0.25)
dev.off()



################################
################################
##fitness correlation plots and models
###############################
###############################

#we fit a full model of all manipulated variables and round# to each of duckweed and microbe "fitness"
#this accounts for any weak n.s. effects of parameters.
#then we use the residuals to ask if there is residual correlation between the two that varies across our parameter levels

mod1.odf <- (MCMCglmm((flod600b)~temp + tire + co210 + tire:temp + tire:co210 + temp:co210 + temp:tire:co210,random=~round ,data=inoconlylive,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))
modfullinocduck <- MCMCglmm(deltapix~temp + co210 + tire + tire:co210 + tire:temp + temp:co210 + tire:co210:temp,random = ~round, data=inoconlylive,verbose=F,pr=T)
inoconlylive$odres <- predict(mod1.odf) - inoconlylive$flod600b
inoconlylive$duckres <- (predict(modfullinocduck) - inoconlylive$deltapix)

#reverse stepwise
summary(MCMCglmm(duckres~odres +odres:temp + odres:co210+ tire:odres + odres:temp:co210 + odres:tire:temp + odres:tire:co210 + odres:tire:co210:temp ,random=~round ,data=inoconlylive,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))
#rm 4 way 2/2; DIC of simpler below essentially not different.
summary(MCMCglmm(duckres~odres +odres:temp + odres:co210+ tire:odres + odres:temp:co210 + odres:tire:temp + odres:tire:co210 ,random=~round ,data=inoconlylive,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))
#rm od:co2:tire 2/2
summary(MCMCglmm(duckres~odres +odres:temp + odres:co210+ tire:odres + odres:temp:co210  + odres:tire:temp,random=~round ,data=inoconlylive,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))
#rm od:tire:temp, 2/2
summary(MCMCglmm(duckres~odres +odres:temp + odres:co210+   tire:odres + odres:temp:co210 ,random=~round ,data=inoconlylive,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))
#rm od:tem:co2 3/3
summary(MCMCglmm(duckres~odres +odres:temp + odres:co210+   tire:odres  ,random=~round ,data=inoconlylive,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))
#od:temp not sig, but  DIC better than simpler mode below; even with multiple runs of each
summary(MCMCglmm(duckres~odres + odres:co210+   tire:odres  ,random=~round ,data=inoconlylive,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))
#stop with the one before the above model, since fit is better. 
#since it was based on a slight difference, testing all simpler models possible, to see if a simpler model fits better
summary(MCMCglmm(duckres~odres +odres:temp + odres:co210  ,random=~round ,data=inoconlylive,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))
summary(MCMCglmm(duckres~odres +odres:temp + odres:tire  ,random=~round ,data=inoconlylive,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))
summary(MCMCglmm(duckres~odres +odres:co210   ,random=~round ,data=inoconlylive,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))
summary(MCMCglmm(duckres~odres +odres:tire   ,random=~round ,data=inoconlylive,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))
summary(MCMCglmm(duckres~odres +odres:temp   ,random=~round ,data=inoconlylive,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))
summary(MCMCglmm(duckres~odres   ,random=~round ,data=inoconlylive,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))
summary(MCMCglmm(duckres~1  ,random=~round ,data=inoconlylive,verbose=F,nitt=500000,burnin=5000,thin=50,pr=T))
#all worse fitting.
#refit best model at higher iterations
mpFmodb <- (MCMCglmm(duckres~odres +odres:temp + odres:co210+   tire:odres  ,random=~round ,data=inoconlylive,verbose=F,nitt=1000000,burnin=10000,thin=100,pr=T))
summary(mpFmodb)

#repeat with temperature subset data
TNmod1.odf <- (MCMCglmm((flod600b)~temp + tire + co210 + tire:temp + tire:co210 + temp:co210 + temp:tire:co210,random=~round ,data=tempNDLinoc,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))
TNmodfullinocduck <- MCMCglmm(deltapix~temp + co210 + tire + tire:co210 + tire:temp + temp:co210 + tire:co210:temp,random = ~round, data=tempNDLinoc,verbose=F,pr=T)
tempNDLinoc$odres <- predict(TNmod1.odf) - tempNDLinoc$flod600b
tempNDLinoc$duckres <- (predict(TNmodfullinocduck) - tempNDLinoc$deltapix)
summary(MCMCglmm(duckres~odres +odres:temp + odres:co210+ tire:odres + odres:temp:co210 + odres:tire:temp + odres:tire:co210 + odres:tire:co210:temp ,random=~round ,data=tempNDLinoc,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))#rm 4way
summary(MCMCglmm(duckres~odres +odres:temp + odres:co210+ tire:odres + odres:temp:co210 + odres:tire:temp + odres:tire:co210  ,random=~round ,data=tempNDLinoc,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))#rm odres CO2 tire
summary(MCMCglmm(duckres~odres +odres:temp + odres:co210+ tire:odres + odres:temp:co210 + odres:tire:temp   ,random=~round ,data=tempNDLinoc,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))#rm od res temp tire
summary(MCMCglmm(duckres~odres +odres:temp + odres:co210+ tire:odres + odres:temp:co210   ,random=~round ,data=tempNDLinoc,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))#rm last 3-way
summary(MCMCglmm(duckres~odres +odres:temp + odres:co210+ tire:odres   ,random=~round ,data=tempNDLinoc,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))#rm odres:temp NOTE THIS IS SLIGHTLY WORSE IN DIC THAN PREVIOUS; 
summary(MCMCglmm(duckres~odres  + odres:co210+ tire:odres   ,random=~round ,data=tempNDLinoc,verbose=F,nitt=100000,burnin=5000,thin=50,pr=T))#stop ; DIC IS NOW BETTER THAN TWO STEPS BACK
#very similar result
mpFmodbSUB <- (MCMCglmm(duckres~odres + odres:co210+   tire:odres  ,random=~round ,data=tempNDLinoc,verbose=F,nitt=1000000,burnin=10000,thin=100,pr=T))
summary(mpFmodbSUB)



#make objects for plotting, reporting means etc
mpFpost <- mpFmodb$Sol
odres.s <- seq(from=min(inoconlylive$odres),to=max(inoconlylive$odres),length.out=1000)
t.s <- seq(from=min(inoconly$temp),to=max(inoconly$temp),length.out=1000)
co2.s <- seq(from = min(inoconly$co210), to = max(inoconly$co210), length.out=1000)

mT <- mean(inoconly$temp) #use 
mhitemp <- mean( inoconly$temp[ inoconly$temp > mT ] ) # mean of the temps above the mean
mlotemp <- mean( inoconly$temp[ inoconly$temp < mT ] ) #
mco2 <- mean(inoconly$co210) ##
mhico2 <- mean( inoconly$co210[ inoconly$co210 > mco2  ] )#mean of the co2 (rownumber) above the mean
mloco2 <- mean( inoconly$co210[ inoconly$co210 < mco2  ] )
predm.mpF.hitc  <- lapply(c(0,5,10), function(tire) sapply(odres.s, function(z) 			   mean( mpFpost[,1] + mpFpost[,2]*z + mpFpost[,3]*z*mhitemp + mpFpost[,4]*z*mhico2 + mpFpost[,5]*z*tire  + rowMeans(mpFpost[,6:8]) ) ) )
predCI.mpF.hitc <- lapply(c(0,5,10), function(tire) sapply(odres.s, function(z) HPDinterval(as.mcmc( mpFpost[,1] + mpFpost[,2]*z + mpFpost[,3]*z*mhitemp + mpFpost[,4]*z*mhico2 + mpFpost[,5]*z*tire  + rowMeans(mpFpost[,6:8]) ),0.9) ) )
predm.mpF.lotc  <- lapply(c(0,5,10), function(tire) sapply(odres.s, function(z) 			   mean( mpFpost[,1] + mpFpost[,2]*z + mpFpost[,3]*z*mlotemp + mpFpost[,4]*z*mloco2 + mpFpost[,5]*z*tire  + rowMeans(mpFpost[,6:8]) ) ) )
predCI.mpF.lotc <- lapply(c(0,5,10), function(tire) sapply(odres.s, function(z) HPDinterval(as.mcmc( mpFpost[,1] + mpFpost[,2]*z + mpFpost[,3]*z*mlotemp + mpFpost[,4]*z*mloco2 + mpFpost[,5]*z*tire  + rowMeans(mpFpost[,6:8]) ),0.9) ) )
##uses MARGINAL HPDI, for pMCMC < 0.1, as described in text.


#6 scenarios (3 tire levels, temp and CO2 above and below thresholds)
pdf("~/Dropbox/TireEcoExp/plnt_micr_Fcor_scenariosXX.pdf",width=6.5,height=3.75)
layout(matrix(c(1,2,3,7,4,5,6,7),byrow=T,ncol=4),widths=c(2,2,2,2.5))
par(oma=c(2,4,0,0))
par(mar=c(0,0,2,0))
plot(inoconlylive$duckres[inoconlylive$temp<mT & inoconlylive$co2 < mco2*100 & inoconlylive$tire == 0]
	~inoconlylive$odres[inoconlylive$temp<mT & inoconlylive$co2 < mco2*100 & inoconlylive$tire == 0],
		pch=16, ylab="",xlab="",main="",ylim=c(-4000,5000),xlim=c(-4,3.5),col=rgb(0,0,1),xaxt="n")#Scenario 1
    polygon( c(odres.s,rev(odres.s)), c(predCI.mpF.lotc[[1]][1,],rev(predCI.mpF.lotc[[1]][2,]) ), border=NA,col=rgb(0,0,1,alpha=0.25)  ) 
    abline(h=0,lty=2)
	lines(predm.mpF.lotc[[1]] ~ odres.s,col=rgb(0,0,1))
    text(-3.5,4500,labels="a)")
par(mar=c(0,0,2,0))
plot(inoconlylive$duckres[inoconlylive$temp<mT & inoconlylive$co2 < mco2*100 & inoconlylive$tire == 5]
	~inoconlylive$odres[inoconlylive$temp<mT & inoconlylive$co2 < mco2*100 & inoconlylive$tire == 5],
		pch=16, ylab="",xlab="",main="",ylim=c(-4000,5000),xlim=c(-4,3.5),col=rgb(0,0,0.5),xaxt="n",yaxt="n")
    polygon( c(odres.s,rev(odres.s)), c(predCI.mpF.lotc[[2]][1,],rev(predCI.mpF.lotc[[2]][2,]) ), border=NA,col=rgb(0,0,0.5,alpha=0.25)  ) 
    abline(h=0,lty=2)
	lines(predm.mpF.lotc[[2]] ~ odres.s,col=rgb(0,0,0.5))
    text(-3.5,4500,labels="b)")
par(mar=c(0,0,2,0))
plot(inoconlylive$duckres[inoconlylive$temp<mT & inoconlylive$co2 < mco2*100 & inoconlylive$tire == 10]
	~inoconlylive$odres[inoconlylive$temp<mT & inoconlylive$co2 < mco2*100 & inoconlylive$tire == 10],
	pch=16, ylab="",xlab="",main="",ylim=c(-4000,5000),xlim=c(-4,3.5),xaxt="n",yaxt="n")#Scenario 3
    polygon( c(odres.s,rev(odres.s)), c(predCI.mpF.lotc[[3]][1,],rev(predCI.mpF.lotc[[3]][2,]) ), border=NA,col=rgb(0,0,0,alpha=0.25)  ) 
	lines(predm.mpF.lotc[[3]] ~ odres.s)
    abline(h=0,lty=2)
    text(-3.5,4500,labels="c)")
par(mar=c(2,0,0,0))
plot(inoconlylive$duckres[inoconlylive$temp>mT & inoconlylive$co2 > mco2*100 & inoconlylive$tire == 0]
	~inoconlylive$odres[inoconlylive$temp>mT & inoconlylive$co2 > mco2*100 & inoconlylive$tire == 0],
		pch=16,xlab="",ylab="",main="",ylim=c(-4000,5000),xlim=c(-4,3.5),col=rgb(1,0,0))#Scenario 2
    polygon( c(odres.s,rev(odres.s)), c(predCI.mpF.hitc[[1]][1,],rev(predCI.mpF.hitc[[1]][2,]) ), border=NA,col=rgb(1,0,0,alpha=0.25)  ) 
    abline(h=0,lty=2)
	lines(predm.mpF.hitc[[1]] ~ odres.s,col=rgb(1,0,0))
    mtext("residual change in pixel area",side=2, adj = -0.65, line=2.5)
    text(-3.5,4500,labels="d)")
par(mar=c(2,0,0,0))
plot(inoconlylive$duckres[inoconlylive$temp>mT & inoconlylive$co2 > mco2*100 & inoconlylive$tire == 5]
	~inoconlylive$odres[inoconlylive$temp>mT & inoconlylive$co2 > mco2*100 & inoconlylive$tire == 5],
		pch=16,xlab="",ylab="",main="",ylim=c(-4000,5000),xlim=c(-4,3.5),col=rgb(0.75,0,0),yaxt="n")#
    polygon( c(odres.s,rev(odres.s)), c(predCI.mpF.hitc[[2]][1,],rev(predCI.mpF.hitc[[2]][2,]) ), border=NA,col=rgb(0.75,0,0,alpha=0.25)  ) 
    abline(h=0,lty=2)
	lines(predm.mpF.hitc[[2]] ~ odres.s,col=rgb(0.75,0,0))
    mtext("residual log optical density",side=1, line=2.5)
    text(-3.5,4500,labels="e)")
par(mar=c(2,0,0,0))
plot(inoconlylive$duckres[inoconlylive$temp>mT & inoconlylive$co2 > mco2*100 & inoconlylive$tire == 10]
	~inoconlylive$odres[inoconlylive$temp>mT & inoconlylive$co2 > mco2*100 & inoconlylive$tire == 10],
		pch=16, ylab="",xlab="",main="",ylim=c(-4000,5000),xlim=c(-4,3.5),col=rgb(0.5,0,0),yaxt="n")#Scenario 4
    polygon( c(odres.s,rev(odres.s)), c(predCI.mpF.hitc[[3]][1,],rev(predCI.mpF.hitc[[3]][2,]) ), border=NA,col=rgb(0.5,0,0,alpha=0.25)  ) 
    abline(h=0,lty=2)
	lines(predm.mpF.hitc[[3]] ~ odres.s,col=rgb(0.5,0,0))
    text(-3.5,4500,labels="f)")
 par(mar=c(0,0,0,0))
 plot(c(0,4)~c(0,8),pch=NA,bty="n",xaxt="n",yaxt="n",ylab="",xlab="")
 	polygon(c(0,1,1,0),c(2,2,2.25,2.25),col=rgb(0,0,1),border=NA)
  	polygon(c(1,2,2,1),c(2,2,2.25,2.25),col=rgb(0,0,0.5),border=NA)
  	polygon(c(2,3,3,2),c(2,2,2.25,2.25),col=rgb(0,0,0),border=NA)
  	polygon(c(0,1,1,0),c(1.75,1.75,2,2),col=rgb(1,0,0),border=NA)
  	polygon(c(1,2,2,1),c(1.75,1.75,2,2),col=rgb(0.75,0,0),border=NA)
  	polygon(c(2,3,3,2),c(1.75,1.75,2,2),col=rgb(0.5,0,0),border=NA)
  	text(c(0.5,1.5,2.5),y=c(2.35,2.35,2.35),labels=c("0x","0.25x","0.5x"),srt=45,adj=0 )
	text(c(3.1,3.1),y=c(2.125,1.875),labels= c(expression('<700ppm, <21.7'*degree*'C'),expression('>700ppm, >21.7'*degree*'C')),adj=0 )
dev.off()


mpFpostTS <- mpFmodbSUB$Sol
odres.ts <- seq(from=min(tempNDLinoc$odres),to=max(tempNDLinoc$odres),length.out=1000)
predm.mpFts.hitc  <- lapply(c(0,5,10), function(tire) sapply(odres.ts, function(z) 			   mean(   mpFpostTS[,1] + mpFpostTS[,2]*z + mpFpostTS[,3]*z*mhico2 + mpFpostTS[,4]*z*tire  + rowMeans(mpFpostTS[,5:7]) ) ) )
predCI.mpFts.hitc <- lapply(c(0,5,10), function(tire) sapply(odres.ts, function(z) HPDinterval(as.mcmc( mpFpostTS[,1] + mpFpostTS[,2]*z + mpFpostTS[,3]*z*mhico2 + mpFpostTS[,4]*z*tire  + rowMeans(mpFpostTS[,5:7]) ),0.9) ) )
predm.mpFts.lotc  <- lapply(c(0,5,10), function(tire) sapply(odres.ts, function(z) 			   mean(   mpFpostTS[,1] + mpFpostTS[,2]*z + mpFpostTS[,3]*z*mloco2 + mpFpostTS[,4]*z*tire  + rowMeans(mpFpostTS[,5:7]) ) ) )
predCI.mpFts.lotc <- lapply(c(0,5,10), function(tire) sapply(odres.ts, function(z) HPDinterval(as.mcmc( mpFpostTS[,1] + mpFpostTS[,2]*z + mpFpostTS[,3]*z*mloco2 + mpFpostTS[,4]*z*tire  + rowMeans(mpFpostTS[,5:7]) ),0.9) ) )
##uses MARGINAL HPDI, for pMCMC < 0.1, as described in text.


pdf("~/Dropbox/TireEcoExp/plnt_micr_Fcor_scenariosXX_TEMPSUBSET.pdf",width=6.5,height=3.75)
layout(matrix(c(1,2,3,7,4,5,6,7),byrow=T,ncol=4),widths=c(2,2,2,2.5))
par(oma=c(2,4,0,0))
par(mar=c(0,0,2,0))
plot(tempNDLinoc$duckres[tempNDLinoc$temp<mT & tempNDLinoc$co2 < mco2*100 & tempNDLinoc$tire == 0]
	~tempNDLinoc$odres[tempNDLinoc$temp<mT & tempNDLinoc$co2 < mco2*100 & tempNDLinoc$tire == 0],
		pch=16, ylab="",xlab="",main="",ylim=c(-4000,5000),xlim=c(-4,3.5),col=rgb(0,0,1),xaxt="n")#Scenario 1
    polygon( c(odres.ts,rev(odres.ts)), c(predCI.mpFts.lotc[[1]][1,],rev(predCI.mpFts.lotc[[1]][2,]) ), border=NA,col=rgb(0,0,1,alpha=0.25)  ) 
    abline(h=0,lty=2)
	lines(predm.mpFts.lotc[[1]] ~ odres.ts,col=rgb(0,0,1))
    text(-3.5,4500,labels="a)")
par(mar=c(0,0,2,0))
plot(tempNDLinoc$duckres[tempNDLinoc$temp<mT & tempNDLinoc$co2 < mco2*100 & tempNDLinoc$tire == 5]
	~tempNDLinoc$odres[tempNDLinoc$temp<mT & tempNDLinoc$co2 < mco2*100 & tempNDLinoc$tire == 5],
		pch=16, ylab="",xlab="",main="",ylim=c(-4000,5000),xlim=c(-4,3.5),col=rgb(0,0,0.5),xaxt="n",yaxt="n")
    polygon( c(odres.ts,rev(odres.ts)), c(predCI.mpFts.lotc[[2]][1,],rev(predCI.mpFts.lotc[[2]][2,]) ), border=NA,col=rgb(0,0,0.5,alpha=0.25)  ) 
    abline(h=0,lty=2)
	lines(predm.mpFts.lotc[[2]] ~ odres.ts,col=rgb(0,0,0.5))
    text(-3.5,4500,labels="b)")
par(mar=c(0,0,2,0))
plot(tempNDLinoc$duckres[tempNDLinoc$temp<mT & tempNDLinoc$co2 < mco2*100 & tempNDLinoc$tire == 10]
	~tempNDLinoc$odres[tempNDLinoc$temp<mT & tempNDLinoc$co2 < mco2*100 & tempNDLinoc$tire == 10],
	pch=16, ylab="",xlab="",main="",ylim=c(-4000,5000),xlim=c(-4,3.5),xaxt="n",yaxt="n")#Scenario 3
    polygon( c(odres.ts,rev(odres.ts)), c(predCI.mpFts.lotc[[3]][1,],rev(predCI.mpFts.lotc[[3]][2,]) ), border=NA,col=rgb(0,0,0,alpha=0.25)  ) 
	lines(predm.mpFts.lotc[[3]] ~ odres.ts)
    abline(h=0,lty=2)
    text(-3.5,4500,labels="c)")
par(mar=c(2,0,0,0))
plot(tempNDLinoc$duckres[tempNDLinoc$temp>mT & tempNDLinoc$co2 > mco2*100 & tempNDLinoc$tire == 0]
	~tempNDLinoc$odres[tempNDLinoc$temp>mT & tempNDLinoc$co2 > mco2*100 & tempNDLinoc$tire == 0],
		pch=16,xlab="",ylab="",main="",ylim=c(-4000,5000),xlim=c(-4,3.5),col=rgb(1,0,0))#Scenario 2
    polygon( c(odres.ts,rev(odres.ts)), c(predCI.mpFts.hitc[[1]][1,],rev(predCI.mpFts.hitc[[1]][2,]) ), border=NA,col=rgb(1,0,0,alpha=0.25)  ) 
    abline(h=0,lty=2)
	lines(predm.mpFts.hitc[[1]] ~ odres.ts,col=rgb(1,0,0))
    mtext("residual change in pixel area",side=2, adj = -0.65, line=2.5)
    text(-3.5,4500,labels="d)")
par(mar=c(2,0,0,0))
plot(tempNDLinoc$duckres[tempNDLinoc$temp>mT & tempNDLinoc$co2 > mco2*100 & tempNDLinoc$tire == 5]
	~tempNDLinoc$odres[tempNDLinoc$temp>mT & tempNDLinoc$co2 > mco2*100 & tempNDLinoc$tire == 5],
		pch=16,xlab="",ylab="",main="",ylim=c(-4000,5000),xlim=c(-4,3.5),col=rgb(0.75,0,0),yaxt="n")#
    polygon( c(odres.ts,rev(odres.ts)), c(predCI.mpFts.hitc[[2]][1,],rev(predCI.mpFts.hitc[[2]][2,]) ), border=NA,col=rgb(0.75,0,0,alpha=0.25)  ) 
    abline(h=0,lty=2)
	lines(predm.mpFts.hitc[[2]] ~ odres.ts,col=rgb(0.75,0,0))
    mtext("residual log optical density",side=1, line=2.5)
    text(-3.5,4500,labels="e)")
par(mar=c(2,0,0,0))
plot(tempNDLinoc$duckres[tempNDLinoc$temp>mT & tempNDLinoc$co2 > mco2*100 & tempNDLinoc$tire == 10]
	~tempNDLinoc$odres[tempNDLinoc$temp>mT & tempNDLinoc$co2 > mco2*100 & tempNDLinoc$tire == 10],
		pch=16, ylab="",xlab="",main="",ylim=c(-4000,5000),xlim=c(-4,3.5),col=rgb(0.5,0,0),yaxt="n")#Scenario 4
    polygon( c(odres.ts,rev(odres.ts)), c(predCI.mpFts.hitc[[3]][1,],rev(predCI.mpFts.hitc[[3]][2,]) ), border=NA,col=rgb(0.5,0,0,alpha=0.25)  ) 
    abline(h=0,lty=2)
	lines(predm.mpFts.hitc[[3]] ~ odres.ts,col=rgb(0.5,0,0))
    text(-3.5,4500,labels="f)")
 par(mar=c(0,0,0,0))
 plot(c(0,4)~c(0,8),pch=NA,bty="n",xaxt="n",yaxt="n",ylab="",xlab="")
 	polygon(c(0,1,1,0),c(2,2,2.25,2.25),col=rgb(0,0,1),border=NA)
  	polygon(c(1,2,2,1),c(2,2,2.25,2.25),col=rgb(0,0,0.5),border=NA)
  	polygon(c(2,3,3,2),c(2,2,2.25,2.25),col=rgb(0,0,0),border=NA)
  	polygon(c(0,1,1,0),c(1.75,1.75,2,2),col=rgb(1,0,0),border=NA)
  	polygon(c(1,2,2,1),c(1.75,1.75,2,2),col=rgb(0.75,0,0),border=NA)
  	polygon(c(2,3,3,2),c(1.75,1.75,2,2),col=rgb(0.5,0,0),border=NA)
  	text(c(0.5,1.5,2.5),y=c(2.35,2.35,2.35),labels=c("0x","0.25x","0.5x"),srt=45,adj=0 )
	text(c(3.1,3.1),y=c(2.125,1.875),labels= c(expression('<700ppm, <21.7'*degree*'C'),expression('>700ppm, >21.7'*degree*'C')),adj=0 )
dev.off()



##for quick calculations of correlation, reported in text
cor(inoconlylive$flod600b,inoconlylive$deltapix)
cor(inoconlylive$odres,inoconlylive$duckres)
summary(MCMCglmm(deltapix~flod600b,random= ~round, data=inoconlylive, verbose=F,pr=T, nitt=1000000,burnin=500,thin=50))
summary(MCMCglmm(duckres~odres,random= ~round, data=inoconlylive, verbose=F,pr=T, nitt=1000000,burnin=500,thin=50))


################################
################################
##TRAIT: greenness
###############################
###############################

summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210 + temp:microbe + co210:microbe + tire:temp:microbe + tire:temp:co210 + temp:co210:microbe + tire:co210:microbe + temp:co210:microbe:tire ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=500,thin=50))
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210 + temp:microbe + co210:microbe + tire:temp:microbe + tire:temp:co210 + temp:co210:microbe + tire:co210:microbe ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=500,thin=50))
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210 + temp:microbe + co210:microbe + tire:temp:microbe + tire:temp:co210 + temp:co210:microbe ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=500,thin=50))
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210 + temp:microbe + co210:microbe + tire:temp:microbe + tire:temp:co210  ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=500,thin=50))
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210 + temp:microbe + co210:microbe +  tire:temp:co210  ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=500,thin=50))
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210 + temp:microbe + co210:microbe   ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=500,thin=50))
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + tire:microbe + co210:tire + temp:co210 + temp:microbe + co210:microbe   ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=500,thin=50))
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + tire:microbe + co210:tire + temp:co210 + temp:microbe   ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=500,thin=50))
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + co210:tire + temp:co210 + temp:microbe   ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=500,thin=50))
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + co210:tire + temp:microbe   ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=500,thin=50))
 # temp microbey marg sig in some runs but not others.
#since DIC was close, test if continuing to remove terms finds a better model.
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + co210:tire   ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=500,thin=50))
 #DIC worse than previous, but close again, continue one more step
summary(MCMCglmm(pergreen~temp + co210 + tire + co210:tire   ,random = ~round, data=nondrylive ,verbose=F,nitt=100000,burnin=500,thin=50))
#DIC now lower than both preceeding, and no more terms can be rm'd. stop.
best.green.LONG <- MCMCglmm(pergreen~temp + co210 + tire + co210:tire   ,random = ~round, data=nondrylive ,verbose=F,nitt=1000000,burnin=500,thin=50,pr=T)
summary(best.green.LONG)


#for the temperature subsetted data:
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210 + temp:microbe + co210:microbe + tire:temp:microbe + tire:temp:co210 + temp:co210:microbe + tire:co210:microbe + temp:co210:microbe:tire ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#rm 4way
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210 + temp:microbe + co210:microbe + tire:temp:microbe + tire:temp:co210 + temp:co210:microbe + tire:co210:microbe ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#rm co2 micr tire
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210 + temp:microbe + co210:microbe + tire:temp:microbe + tire:temp:co210 + temp:co210:microbe ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#rm temp co2 micr
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210 + temp:microbe + co210:microbe + tire:temp:microbe + tire:temp:co210  ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#rm temp micr tire
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210 + temp:microbe + co210:microbe + tire:temp:co210  ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#rm tire temp co2
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210 + temp:microbe + co210:microbe  ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#rm co2 micr
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + tire:microbe + temp:tire + co210:tire + temp:co210 + temp:microbe   ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#rm micr tire
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + temp:tire + co210:tire + temp:co210 + temp:microbe   ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#rm temp co2
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + temp:tire + co210:tire + temp:microbe   ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#rm temp tire
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + co210:tire + temp:microbe   ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#rm temp micr
summary(MCMCglmm(pergreen~temp + co210 + microbe + tire + co210:tire   ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#rm micr. NOTE DIC SLIGHTLY WORSE THAN PREV
summary(MCMCglmm(pergreen~temp + co210  + tire + co210:tire   ,random = ~round, data=tempNDL ,verbose=F,nitt=100000,burnin=500,thin=50))#stop, and better than two prev models.
#same result.
best.green.LONGSUB <- MCMCglmm(pergreen~temp + co210 + tire + co210:tire   ,random = ~round, data=tempNDL ,verbose=F,nitt=1000000,burnin=500,thin=50,pr=T)
summary(best.green.LONGSUB)


#means, HPDI, predictions, etc across conditions for reporting directly or for figures
grnpost <- best.green.LONG$Sol
t.s <- seq(from = min(nondrylive$temp), to = max(nondrylive$temp),length.out=1000)
co.s <- seq(from = min(nondrylive$co210), to = max(nondrylive$co210),length.out=1000)
mnT <- mean(nondrylive$temp)
mnCO2 <- mean(nondrylive$co210)
pred.grntmp.mnctr <- sapply(t.s, function(z)                mean(grnpost[,1] + grnpost[,2]*z + grnpost[,3]*mnCO2 + grnpost[,4]*5 + grnpost[,5]*5*mnCO2 + rowMeans(grnpost[,6:8])) ) 
pred.grntmp.CIctr <- sapply(t.s, function(z) HPDinterval(as.mcmc(grnpost[,1] + grnpost[,2]*z + grnpost[,3]*mnCO2 + grnpost[,4]*5 + grnpost[,5]*5*mnCO2 + rowMeans(grnpost[,6:8])),0.95 )) 
pred.grntmp.mnctr[c(1,length(t.s))]
#[1] 0.4081532 0.4399450
pred.grntmp.CIctr[,c(1,length(t.s))]
#           [,1]      [,2]
# [1,] 0.4037222 0.4358224
# [2,] 0.4128219 0.4439575
pred.grnco2tr.mnt <- lapply(c(0,5,10), function(tire) sapply(co.s, function(z)                mean(grnpost[,1] + grnpost[,2]*mnT + grnpost[,3]*z + grnpost[,4]*tire + grnpost[,5]*tire*z + rowMeans(grnpost[,6:8])) ) )
pred.grnco2tr.CI <- lapply(c(0,5,10), function(tire) sapply(co.s, function(z) HPDinterval(as.mcmc(grnpost[,1] + grnpost[,2]*mnT + grnpost[,3]*z + grnpost[,4]*tire + grnpost[,5]*tire*z + rowMeans(grnpost[,6:8])),0.95 )) )
pred.grnco2tr.mnt[[1]][c(1,length(co.s))]
#[1] 0.4326884 0.4136973
pred.grnco2tr.CI[[1]][,c(1,length(co.s))]
#           [,1]      [,2]
# [1,] 0.4289504 0.4100834
# [2,] 0.4363188 0.4174099
pred.grnco2tr.mnt[[3]][c(1,length(co.s))]
# [1] 0.4317752 0.4215046
pred.grnco2tr.CI[[3]][,c(1,length(co.s))]
#           [,1]      [,2]
# [1,] 0.4280462 0.4178591
# [2,] 0.4353528 0.4250936
grnmns <- tapply(nondrylive$pergreen,paste(nondrylive$tire,nondrylive$co2,sep="."),mean)
grnses <- tapply(nondrylive$pergreen,paste(nondrylive$tire,nondrylive$co2,sep="."),std.error)
grntire <- tapply(nondrylive$tire,paste(nondrylive$tire,nondrylive$co2,sep="."),mean)
grnco2 <- tapply(nondrylive$co2,paste(nondrylive$tire,nondrylive$co2,sep="."),mean)
distmn.tire <- range01(abs(nondrylive$tire - 5))
distmn.co2 <- range01(abs(nondrylive$co210-mnCO2))
distmn.both <- (distmn.tire+distmn.co2)/2
distL.tire <- range01(nondrylive$tire)
distL.co2 <- range01(nondrylive$co210)
distL.both <- (distL.co2+distL.tire)/2

pdf("~/Dropbox/TireEcoExp/best_greennessX.pdf",width=3,height=6)
par(mfrow=c(2,1))
par(mar=c(4,4,1,1))
par(oma=c(0,0,1,0))
plot(pergreen~temp,pch=16,cex=1-0.75*distmn.both,ylab="",xlab="",data=nondrylive,col=rgb(0,0,0,alpha= 0.75*(1-distmn.both/2) ))# & co2 > 600 & co2 < 800
	polygon( c(t.s,rev(t.s)), c( pred.grntmp.CIctr[1,] , rev(pred.grntmp.CIctr[2,])  ) , col = rgb(0,0,0,alpha=0.25),border=NA)
	lines(pred.grntmp.mnctr~t.s,col=rgb(0,0,0),lwd=2)
	mtext( "greenness",side=2,line=2.5)
	mtext(expression("Temperature"~degree*'C'),side=1,line=2.5)
	mtext("a)", side =3, adj=-0.4, line=0.5)
plot(grnmns[grntire==0]~grnco2[grntire==0],ylim=c(0.405,0.45),xlim=c(380,1040),pch=NA,ylab="",xlab="")
	polygon( c(co.s*100,rev(co.s*100)), c( pred.grnco2tr.CI[[1]][1,] , rev(pred.grnco2tr.CI[[1]][2,])  ) , col = rgb(0,0,1,alpha=0.25),border=NA)
	polygon( c(co.s*100,rev(co.s*100)), c( pred.grnco2tr.CI[[3]][1,] , rev(pred.grnco2tr.CI[[3]][2,])  ) , col = rgb(0,0,0,alpha=0.25),border=NA)
	lines(pred.grnco2tr.mnt[[1]]~I(co.s*100),col=rgb(0,0,1),lwd=2)
	lines(pred.grnco2tr.mnt[[3]]~I(co.s*100),col=rgb(0,0,0),lwd=2)
	points(grnmns[grntire==0]~grnco2[grntire==0],pch=16,col=rgb(0,0,1))
	arrows(grnco2[grntire==0], y0=grnmns[grntire==0] -  grnses[grntire==0], y1=grnmns[grntire==0] +  grnses[grntire==0] ,length=0,col=rgb(0,0,1))
	points(grnmns[grntire==10]~I(grnco2[grntire==10]+10),pch=16,col=rgb(0,0,0))
	arrows(grnco2[grntire==10]+10, y0=grnmns[grntire==10] -  grnses[grntire==10], y1=grnmns[grntire==10] +  grnses[grntire==10] ,length=0,col=rgb(0,0,0))
	mtext("greenness",side=2,line=2.5)
	mtext(expression('CO'[2]~'ppm'),side=1,line=2.5)
	legend(700,y=0.45,legend=c("0x","0.5x"),fill=c(rgb(0,0,1),rgb(0,0,0)),bty="n")	
	mtext("b)", side =3, adj=-0.4, line=0.5)
dev.off()


