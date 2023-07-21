##########################################################################################
##Compaire normalized NND
##R version 3.5.2 (2018-12-20) -- "Eggshell Igloo")
##Hugo Lachuer (last update 21/10/2020)
##########################################################################################

########################################################################################
## 0/ Notes


########################################################################################
## 1/ Library & functions

if("ggplot2" %in% rownames(installed.packages()) == FALSE){
	install.packages("ggplot2")
}

if("lawstat" %in% rownames(installed.packages()) == FALSE){
	install.packages("lawstat")
}

if("dunn.test" %in% rownames(installed.packages()) == FALSE){
	install.packages("dunn.test")
}

library(tcltk)	#User interface
library(ggplot2)	#Plot
library(lawstat)	#Levene's test
library(dunn.test)	#Dunn testing

#Computing NND
ComputeNND <- function(Pattern){
	NND <- vector(, length=nrow(Pattern))
	for(i in 1:nrow(Pattern)){
		Index <- 1:nrow(Pattern)
		Index <- Index[-i]
		NND[i] <- min(sqrt((Pattern[i,1] - Pattern[Index,1])^2 + (Pattern[i,2] - Pattern[Index,2])^2 + (Pattern[i,3] - Pattern[Index,3])^2))
	}
	return(mean(NND))
}

#Computing IOD
ComputeIOD <- function(Pattern){
	IOD <- matrix(, ncol=nrow(Pattern), nrow=nrow(Pattern)-1)
	for(i in 1:nrow(Pattern)){
		Index <- 1:nrow(Pattern)
		Index <- Index[-i]
		k <- 1
		for(j in Index){
			IOD[k,i] <- sqrt((Pattern[i,1] - Pattern[j,1])^2 + (Pattern[i,2] - Pattern[j,2])^2 + (Pattern[i,3] - Pattern[j,3])^2)
			k <- k + 1
		}
	}
	return(mean(IOD))
}

#Computing SD of NND
SD_NND <- function(Pattern){
	NND <- vector(, length=nrow(Pattern))
	for(i in 1:nrow(Pattern)){
		Index <- 1:nrow(Pattern)
		Index <- Index[-i]
		NND[i] <- min(sqrt((Pattern[i,1] - Pattern[Index,1])^2 + (Pattern[i,2] - Pattern[Index,2])^2 + (Pattern[i,3] - Pattern[Index,3])^2))
	}
	return(sd(NND))
}

#Computing SD of IOD
SD_IOD <- function(Pattern){
	IOD <- matrix(, ncol=nrow(Pattern), nrow=nrow(Pattern)-1)
	for(i in 1:nrow(Pattern)){
		Index <- 1:nrow(Pattern)
		Index <- Index[-i]
		k <- 1
		for(j in Index){
			IOD[k,i] <- sqrt((Pattern[i,1] - Pattern[j,1])^2 + (Pattern[i,2] - Pattern[j,2])^2 + (Pattern[i,3] - Pattern[j,3])^2)
			k <- k + 1
		}
	}
	return(sd(IOD))
}

#Test and plot
SpatialTesting <- function(DF, YName){

	#Check homoscedasticity and normality
	pLevene <- levene.test(DF[,2],DF[,1],location='median')$p.value	#Homoscedasticity test
	LM <- aov(DF[,2] ~ DF[,1], data=DF)
	Residuals <- residuals(object = LM)
	pShapiro <- shapiro.test(x = Residuals )$p.value	#Normality test
	
	
	#Test
	if(pLevene >= 10^-4){
		print(paste0("Levene's test for homoscedasticity : p-value = ", round(pLevene, 4)))
	} else {
		print("Levene's test for homoscedasticity : p-value < 0.0001")
	}	

	if(pShapiro >= 10^-4){
		print(paste0("Shapiro-Wilk test for normality : p-value = ", round(pShapiro, 4)))
	} else {
		print("Shapiro-Wilk test for normality : p-value < 0.0001")
	}

	if( (pShapiro>=0.05) & (pLevene>=0.05) ){
		print("Normality and homoscedasticity assumptions meet. Statistical testing by ANOVA with Tukey post-hoc test")
		p.value <- summary(LM)[[1]][["Pr(>F)"]][1]
		TukeyHSD(LM)
	} else {
		print("Normality and homoscedasticity assumptions not meet. Statistical testing by Kruskal-Wallis with Dunn-Sidak correction")
		DF[,1] <- paste0(substr(DF[,1], 1, 4), " ", substr(DF[,1], 11, nchar(DF[,1])))
		DF[,1] <- factor(DF[,1], ordered = TRUE)
		p.value <- kruskal.test(DF[,2] ~ DF[,1], data=DF)$p.value
		dunn.test(x=DF[,2], g=DF[,1], method = "sidak", altp=TRUE)
	}
	
	p <- ggplot(DF, aes(x=DF[,1], y=DF[,2], fill="Gray")) + geom_boxplot(outlier.shape = NA, fill='#A4A4A4', color="black") + geom_jitter(color="black", size=2, alpha=1, width = 0.25)
	p <- p + theme(axis.text.x = element_text(size=18))+theme(axis.title.x = element_text(size=18))+theme(axis.text.y = element_text(size=18))+theme(axis.title.y = element_text(size=18))
	p <- p + theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
	p <- p + theme(plot.title = element_text(size=22)) + theme(legend.position = "none")
	if(p.value < 10^-4){
		p <- p + ggtitle("p.value < 0.0001") + xlab("Condition") + ylab(YName) + theme(plot.title = element_text(size=22))
	} else {
		p <- p + ggtitle(paste0("p.value = ", round(p.value,4))) + xlab("Condition") + ylab(YName) + theme(plot.title = element_text(size=22))
	}
	return(p)
}


########################################################################################
## 2/ User interface and data importation


NCond <- tclVar("3")
NNDcomput <- tclVar("YES")
NormalizeDistances <- tclVar("YES")
IODcomput <- tclVar("YES")
DistanceToPoint <- tclVar("YES")
BarycenterComp <- tclVar("YES")
PixelSize <- tclVar("0.1075")
dz <- tclVar("0.20")
Nstruc <- tclVar("100")

tt <- tktoplevel()
tkwm.title(tt,"parameters")

NCond.entry <- tkentry(tt, textvariable=NCond)
NNDcomput.entry <- tkentry(tt, textvariable=NNDcomput)
NormalizeDistances.entry <- tkentry(tt, textvariable=NormalizeDistances)
IODcomput.entry <- tkentry(tt, textvariable=IODcomput)
DistanceToPoint.entry <- tkentry(tt, textvariable=DistanceToPoint)
BarycenterComp.entry <- tkentry(tt, textvariable=BarycenterComp)
PixelSize.entry <- tkentry(tt, textvariable=PixelSize)
dz.entry <- tkentry(tt, textvariable=dz)
Nstruc.entry <- tkentry(tt, textvariable=Nstruc)


reset <- function() {
	tclvalue(NCond)<-""
	tclvalue(NNDcomput)<-""
	tclvalue(NormalizeDistances)<-""
	tclvalue(IODcomput) <- ""
	tclvalue(DistanceToPoint)<-""
	tclvalue(BarycenterComp)<-""
	tclvalue(Magnification)<-""
	tclvalue(PixelSize)<-""
	tclvalue(dz)<-""
	tclvalue(Nstruc)<-""
}
reset.but <- tkbutton(tt, text="Reset", command=reset)
done <- tclVar(0)
OK <- function(){
	tclvalue(done)<-1
	tkdestroy(tt)
}
OK.but <- tkbutton(tt,text="OK",command=OK)
tkgrid(tklabel(tt,text="Enter parameters (at least 2)"),columnspan=3, pady = 10)
tkgrid(tklabel(tt,text="Enter the number of condition"), NCond.entry, pady= 10, padx= 10)
tkgrid(tklabel(tt,text="Do you want to normalize distances (YES or NO)"), NormalizeDistances.entry, pady= 10, padx= 10)
tkgrid(tklabel(tt,text="Do you want to compute NND (Nearest Neighbor Distance) (YES or NO)"), NNDcomput.entry, pady= 10, padx= 10)
tkgrid(tklabel(tt,text="Do you want to compute IOD (Inter Organelle Distance) (YES or NO)"), IODcomput.entry, pady= 10, padx= 10)
tkgrid(tklabel(tt,text="Do you want to compute distance to arbitrary point (cell center for example) (YES or NO)"), DistanceToPoint.entry, pady= 10, padx= 10)
tkgrid(tklabel(tt,text="Use barycenter as arbitrary point for previous distance (YES or NO)"), BarycenterComp.entry, pady= 10, padx= 10)
tkgrid(tklabel(tt,text="Enter the pixel size (edge of one pixel) (in µm)"), PixelSize.entry, pady= 10, padx= 10)
tkgrid(tklabel(tt,text="Enter the z shift (in µm)"), dz.entry, pady = 10, padx = 10)
tkgrid(tklabel(tt,text="Minimal number of structures required to analyze one cell"), Nstruc.entry, pady = 10, padx = 10)
tkgrid(reset.but, OK.but, pady= 10, padx= 10)
tkwait.variable(done)

NCond <- as.numeric(tclvalue(NCond))
NNDcomput <- tclvalue(NNDcomput)
ifelse(NNDcomput == "YES", NNDcomput <- TRUE, NNDcomput <- FALSE)
IODcomput <- tclvalue(IODcomput)
ifelse(IODcomput == "YES", IODcomput <- TRUE, IODcomput <- FALSE)
DistanceToPoint <- tclvalue(DistanceToPoint)
ifelse(DistanceToPoint == "YES", DistanceToPoint <- TRUE, DistanceToPoint <- FALSE)
NormalizeDistances <- tclvalue(NormalizeDistances)
ifelse(NormalizeDistances == "YES", NormalizeDistances <- TRUE, NormalizeDistances <- FALSE)
BarycenterComp <- tclvalue(BarycenterComp)
ifelse(BarycenterComp == "YES", BarycenterComp <- TRUE, BarycenterComp <- FALSE)
PixelSize <- as.numeric(tclvalue(PixelSize))
dz <- as.numeric(tclvalue(dz))
Nstruc <- as.numeric(tclvalue(Nstruc))

path <- tk_choose.dir("Choose directory to save results")
setwd(path)


MetaList <- list()
Label <- vector()
Ncell <- 0
M <- 0
for(i in 1:NCond){
	NRepetition <- tclVar("3")
	tt <- tktoplevel()
	tkwm.title(tt,"parameters")
	NRepetition.entry <- tkentry(tt, textvariable=NRepetition)
	reset.but <- tkbutton(tt, text="Reset", command=reset)
	done <- tclVar(0)
	OK.but <- tkbutton(tt,text="OK",command=OK)
	tkgrid(tklabel(tt,text=paste0("Number of repetition (=files) for condition ",i)), NRepetition.entry, pady= 10, padx= 10)
	tkgrid(reset.but, OK.but, pady= 10, padx= 10)
	tkwait.variable(done)
	NRepetition <- as.numeric(tclvalue(NRepetition))
	
	#Pool in a list
	List <- list()
	k <- 1
	for(j in 1:NRepetition){
	
		#Importation
		Segmentation <- tk_choose.files(default = "", caption = paste0("Select segmentation .txt file number ", j, " for condition ", i))
		Data <- read.csv(Segmentation, sep="\t", header=FALSE)
		Data <- rbind(Data, rep(NA, ncol(Data)), c("Index", rep(NA, ncol(Data)-1)))
		Pos <- which(Data[,1] == "Index")

		for(l in 1:(length(Pos)-1)){
			A1 <- Data[(Pos[l] + 1):(Pos[l+1]-2), c(3, 13, 12, 11)]	#CentroidX=z, centroidY=y, CendroidZ=x
			A1 <- as.matrix(A1)
			A1 <- as.numeric(A1)
			A1 <- matrix(A1, ncol=4)
			colnames(A1) <- c("Volume", "X", "Y", "Z")
			A1[,1] <- A1[,1] * dz*PixelSize^2	#Convert in µm^3
			A1[,2:3] <- A1[,2:3] * PixelSize	#Convert in µm
			A1[,4] <- A1[,4] * dz
			List[[k]] <- A1
			k <- k + 1
			
			Label[Ncell+1] <- paste0("Condition ", i, " Repetition ", j, " Cell ", l)
			Ncell <- Ncell + 1
		}
	}
	MetaList[[i]] <- List
}

if(BarycenterComp==FALSE){
	CenterList <- list()
	M <- 0
	for(i in 1:NCond){
		NRepetition <- tclVar("3")
		tt <- tktoplevel()
		tkwm.title(tt,"parameters")
		NRepetition.entry <- tkentry(tt, textvariable=NRepetition)
		reset.but <- tkbutton(tt, text="Reset", command=reset)
		done <- tclVar(0)
		OK.but <- tkbutton(tt,text="OK",command=OK)
		tkgrid(tklabel(tt,text=paste0("Number of repetition (=files) for condition ",i)), NRepetition.entry, pady= 10, padx= 10)
		tkgrid(reset.but, OK.but, pady= 10, padx= 10)
		tkwait.variable(done)
		NRepetition <- as.numeric(tclvalue(NRepetition))
		
		#Pool in a list
		for(j in 1:NRepetition){
			M <- M + 1
			A1 <- tk_choose.files(default = "", caption = paste0("Select center position .csv file number ", j, " for condition ", i))
			A2 <- read.csv(A1)
			if(ncol(A2)<3){
				A2 <- read.delim(A1)
			}
			CenterList[[M]] <- A2
		}
	}


	CenterMatrix <- CenterList[[1]]
	CenterMatrix <- CenterMatrix[,which(colnames(CenterMatrix) == "X" | colnames(CenterMatrix) == "Y" | colnames(CenterMatrix) == "Z")]
	for(i in 2:length(CenterList)){
		A1 <- CenterList[[i]]
		A2 <- A1[,which(colnames(A1) == "X" | colnames(A1) == "Y" | colnames(A1) == "Z")]
		CenterMatrix <- rbind(CenterMatrix, A2)
	}
	CenterList <- CenterMatrix
	CenterList[,1:2] <- CenterList[,1:2] * PixelSize
	CenterList[,3] <- (CenterList[,3]-1)*dz
}


########################################################################################
## 3/ Minimal number of structure and remove cells with < 100 structures

#Number quantification
Nmin <- 10^6
for(i in 1:length(MetaList)){
	List <- MetaList[[i]]
	for(j in 1:length(List)){
		Nmin <- min(Nmin, nrow(List[[j]]))
	}
}

Nmin <- max(Nstruc, Nmin)	#No less than 100 structures
print(paste0("minimal number of structure : ", Nmin))

#Remove cells <Nmin
NRemoved <- 0
NCell <- 0
m <- 0
M <- length(Label)
for(i in length(MetaList):1){
	List <- MetaList[[i]]
	for(j in length(List):1){
		if(nrow(List[[j]])<Nmin){
			List[[j]] <- NULL
			NRemoved <- NRemoved + 1
			Label <- Label[-(M-m)]
		} else {
			NCell <- NCell + 1
		}
		m <- m + 1
	}
	MetaList[[i]] <- List
}

print(paste0(NRemoved, " removed cells (because number of structure < 100)"))


########################################################################################
## 4/ Compute Number of structures and average volume

Result <- as.data.frame(matrix(, ncol=10, nrow=NCell))
colnames(Result) <- c("Label", "Number of structure", "Average volume", "NND", "IOD", "Distance to point", "Inertia", "Barycenter", "Normalized", "Nmin")

Result[,1] <- Label
if(NormalizeDistances){
	Result[,9] <- rep("YES", nrow(Result))
} else {
	Result[,9] <- rep("NO", nrow(Result))
}
Result[,10] <- rep(Nmin, nrow(Result))

M <- 0
for(i in 1:length(MetaList)){
	List <- MetaList[[i]]
	for(j in 1:length(List)){
			A1 <- List[[j]]
			M <- M + 1
			Result[M,2] <- nrow(A1)
			Result[M,3] <- mean(A1[,1])
	}
}

########################################################################################
## 5/ Compute NND and IOD

m <- 1
for(i in 1:length(MetaList)){
	List <- MetaList[[i]]
	for(j in 1:length(List)){
		A1 <- List[[j]]
		A1 <- A1[,-1]
			
		if(NormalizeDistances){
			NND_Est <- IOD_Est <- 0
			A2 <- A1[sample(1:nrow(A1), Nmin, replace=FALSE),]
			N1 <- ((SD_NND(A2)/sqrt(Nmin))/(0.01*ComputeNND(A2)))^2	#SEM²/µ² = 1%
			N2 <- ((SD_IOD(A2)/sqrt(Nmin))/(0.01*ComputeIOD(A2)))^2	#SEM²/µ² = 1%
			N <- max(N1, N2)
			N <- round(N)
			for(k in 1:N){
				A2 <- A1[sample(1:nrow(A1), Nmin, replace=FALSE),]
				NND_Est <- NND_Est + ComputeNND(A2)
				IOD_Est <- IOD_Est + ComputeIOD(A2)
			}
			Result[m,4] <- NND_Est/N
			Result[m,5] <- IOD_Est/N
		} else{
			Result[m,4] <- ComputeNND(A1)
			Result[m,5] <- ComputeIOD(A1)
		}
		m <- m + 1
	}
}

########################################################################################
## 5/ Distance to point, Inertia and barycenter

m <- 1
for(i in 1:length(MetaList)){
	List <- MetaList[[i]]
	for(j in 1:length(List)){
		A1 <- List[[j]]
		A1 <- A1[,-1]
		Barycenter <- c(mean(A1[,1]), mean(A1[,2]), mean(A1[,3]))	#Barycenter
		if(BarycenterComp){
			Point <- Barycenter
			Result[m,8] <- "YES"
		} else {
			Point <- CenterList[m,]
			Point <- as.vector(Point)
			Point <- as.numeric(Point)
			Result[m,8] <- "NO"
		}
		
		V1 <- sqrt((A1[,1] - Point[1])^2 + (A1[,2] - Point[2])^2 + (A1[,3] - Point[3])^2)
		V2 <- (A1[,1] - Barycenter[1])^2 + (A1[,2] - Barycenter[2])^2 + (A1[,3] - Barycenter[3])^2
		Result[m,6] <- mean(V1)
		Result[m,7] <- mean(V2)
		m <- m + 1
	}
}

setwd(path)
write.table(Result, file="Result.txt", sep="\t", row.names=FALSE)

########################################################################################
## 6/ Statistics

#Condition vector
CondVector <- Label
for(i in 1:length(CondVector)){
	A1 <- CondVector[i]
	k <- 13
	while(substr(A1,k,k+9)!="Repetition"){
		k <- k + 1
	}
	CondVector[i] <- substr(A1, 1, k-2)
}

colnames(Result) <- c("Label", "Number of structure", "Average volume", "NND", "IOD", "Distance to point", "Inertia", "Barycenter", "Normalized", "Nmin")


print("Test on the number of structures")
Number_DF <- cbind(CondVector, Result[,2])
Number_DF <- as.data.frame(Number_DF)
Number_DF[,2] <- Result[,2]
Number_DF[,1] <- CondVector
p <- SpatialTesting(DF=Number_DF, YName= "Number of structures")
p
ggsave("Number analysis.pdf", plot = p)


print("Test on the volume of structures")
Volume_DF <- cbind(CondVector, Result[,3])
Volume_DF <- as.data.frame(Volume_DF)
Volume_DF[,2] <- Result[,3]
Volume_DF[,1] <- CondVector
p <- SpatialTesting(DF=Volume_DF, YName= "Volume of structures (µm^3)")
p
ggsave("Volume analysis.pdf", plot = p)


print("Test on the 	NNDs")
NND_DF <- cbind(CondVector, Result[,4])
NND_DF <- as.data.frame(NND_DF)
NND_DF[,2] <- Result[,4]
NND_DF[,1] <- CondVector
p <- SpatialTesting(DF=NND_DF, YName= "NND (µm)")
p
ggsave("NND analysis.pdf", plot = p)


print("Test on the 	IODs")
IOD_DF <- cbind(CondVector, Result[,5])
IOD_DF <- as.data.frame(IOD_DF)
IOD_DF[,2] <- Result[,5]
IOD_DF[,1] <- CondVector
p <- SpatialTesting(DF=IOD_DF, YName= "IOD (µm)")
p
ggsave("IOD analysis.pdf", plot = p)


print("Test on the 	Distance to point")
DistanceToPoint_DF <- cbind(CondVector, Result[,6])
DistanceToPoint_DF <- as.data.frame(DistanceToPoint_DF)
DistanceToPoint_DF[,2] <- Result[,6]
DistanceToPoint_DF[,1] <- CondVector
p <- SpatialTesting(DF=DistanceToPoint_DF, YName= "Distance to point (µm)")
p
ggsave("Distance to point analysis.pdf", plot = p)


print("Test on Inertia")
Inertia_DF <- cbind(CondVector, Result[,7])
Inertia_DF <- as.data.frame(Inertia_DF)
Inertia_DF[,2] <- Result[,7]
Inertia_DF[,1] <- CondVector
p <- SpatialTesting(DF=Inertia_DF, YName= "Inertia")
p
ggsave("Inertia analysis.pdf", plot = p)