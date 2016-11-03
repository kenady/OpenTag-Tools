##-----------------------Calculate Dive Statistics for OpenTag Data---------------------##

divestats <- function(ptmp,ulim,blim,sn,fs){
# ptmp is a data frame of your time depth data from OpenTag. Create this using otload.R
# ulim is the upper limit of possible surface values (positive value). This is used in the zero-correction depth process.
# blim in the lower limit of possible surface values (negative value). This is used in the zero-correction depth process.
# sn is the surface interval within which surface noise should be zero (i.e., values less than 1 should be 0)
# dr is the proportion of max depth below which bottom phase begins
# fs is the sampling frequency of your data

# This function will calculate the following dive statistics for OpenTag pressure sensor data:
  # Dive number: divenum
  # Dive start time: start
  # Dive end time: end
  # Maximum dive depth (m): maxdepth
  # Dive duration (sec): Tdown
  # Post-dive surface duration (sec): Tup
  # Duration of time spent at bottom (sec): bottom.duration # bottom is calculated as the time between the first and last local minimum depth
  # Number of wiggles at bottom: wig.n
  # Descent duration (sec): descent
  # Ascent duration (sec): ascent
  # Descent rate (m/sec): descent.rate
  # Ascent rate (m/sec): ascent.rate

# Example: dives <- divestats(ptmp,2,-2,0.2,1)
# Optional plotting option an the end of the code. If left in, each dive will be plotted, with 
# red points representing the descent and ascent portion of a dive, blue points representing the 
# bottom phase of a dive, and green vertical lines indicating the start and end points of the bottom 
# duration. Comment this section out if you do not want to plot the dives (this will up the function).

# Modified from code written by Mark Hindell, Mark.Hindell@utas.edu.au, for dealing with dive 
# data collected from Wildlife Computer satellite tags 

# Reny Tyson, rtyson@mote.org                                                            
# 26 October 2016                                                                         

############################################################################################
#Step 1: Zero correct depth and identify individual dives

# Depth needs to be in POSIX structure
ptmp$datetime <- strptime(ptmp$datetime,"%Y-%m-%d %H:%M:%S")  
    
# Depth needs to be positive values
ptmp$depth <- ptmp$depth * -1
  
# Generate a zocing function based on modal surface depth 
# Corrects all depth by the zoc value. 
# Includes the depth range from which to develop the depth histogram
# Optional plotting option
zoc <- function(x,ulim,blim){
  h <- hist(x[x < ulim & x > blim], breaks = seq(blim, ulim, .1), plot=FALSE) # To plot say plot=TRUE
  #abline(v=h$breaks[which.max(h$counts)], col="red") # To plot turn off comment
  #print(h$breaks[which.max(h$counts)]) # To plot turn off comment
  ##scan("", 1) # To plot turn off comment
  x - h$breaks[which.max(h$counts)]
}
ptmp$depth_zoc <- zoc(ptmp$depth,ulim,blim)

# Set surface interval within the surface noise value to 0
ptmp$depth_nzoc <- ptmp$depth_zoc
for (i in 1:nrow(ptmp)){
  if(ptmp$depth_zoc[i] < sn){ptmp$depth_nzoc[i] <- 0}
}

# Assigns numbers to each dive based on a user defined "surface noise" (sn) value 
# Gives each surface or dive section a consecutive number
num <- ceiling(cumsum(abs(c(0, diff(ptmp$depth_nzoc == 0))))/2)
ptmp$num <- num

# Optional plot to look at how dives are broken up (this can help you set your sn)
#ggplotly(qplot(datetime,-depth_nzoc,data=ptmp,geom='line',colour=factor(num)))

# Add max depth to dives
dives <- ptmp[,c('datetime','depth_nzoc','num')]
max.d <- aggregate(depth_nzoc~num, data=dives, max) #find maximum depth of each dive
dives <- merge(dives, max.d, by.x="num", by.y="num", all.x=T)
colnames(dives) <- c("num", "datetime", "depth", "max.dep") #rename the variables
dives <- dives[dives$num > 0,] 

###############################################################################
# Step 2. extract dive paramaters

# Define the unique dives in the file
dive.list <- unique(dives$num)
n <- length(dive.list)## the number of dives in the record

# Create an empty file into which to write the dive data
dout <- data.frame("divenum"=rep(0,n), "start"=0, "end"=0,"maxdepth"=0,"Tdown"=0, "Tup"=0, "bottom.duration"=0,"wig.n"=0,'descent'=0,'ascent'=0,'descent.rate'=0,'ascent.rate'=0)

# Set time structure of empty file
dout$start <- strptime(dout$start,"%Y-%m-%d %H:%M:%S")
dout$end <- strptime(dout$end,"%Y-%m-%d %H:%M:%S")

for(i in 1:length(dive.list)){
  
  # Identify start and end of dive...
  # Need to add a zero to first record
  first <- dives[dives$num==dive.list[i],][1,]
  first["datetime"] <- first["datetime"]-fs
  first["depth"] <- 0  
  
  # ...and to the last record
  last <- dives[dives$num==dive.list[i],][nrow(dives[dives$num==dive.list[i],]),]
  last["datetime"] <- last["datetime"]+fs
  last["depth"] <- 0
  x. <-rbind(first, dives[dives$num==dive.list[i],], last) # so you bind a zero to begining and end of dive with the sampling time offest (+/- four included)
  
  # Add a seconds column
  for (k in 2:nrow(x.)){
    x.$sec[1] <- 0
    x.$sec[k] <- x.$sec[k-1] + fs
  }
  
  # Find start and end of dive
  start <- which(diff(x.$depth > 0)==1) # which depths have lagged differences greater than mind (4)? just sets the begining of dive since the first value should be greater than 4 m (everything shallower we set to zero)
  end <- which(diff(x.$depth > 0)==-1)+1 # sets the end point as the point before 
  x <- x.[start:end,] # cuts off the surface intervals
  
  # Identify the points within your "bottom range" calculated as the first and last local minimum
  # Find the first local minimum
  inf1 <- c()
  for (j in 2:nrow(x)){
        g <- x$depth[j] > x$depth[j-1]
       inf1 <- c(inf1,g)
  }
  inf1 <- c(0,diff(inf1))
  dp <- which(inf1 == -1)[1] # This is the row of the first minimum
    
  # Find the last local minimum
  x2 <- x[order(x$datetime, decreasing=T),] # dive data in reverse
  inf2 <- c()
  for (k in 2:nrow(x)){
      g <- x2$depth[k] > x2$depth[k-1]
      inf2 <- c(inf2,g)
  }
  inf2 <- c(0,diff(inf2))
  ap <- which(inf2== -1)[1]   # This is the row of the last minimum
  ap <- which(x$datetime == x2$datetime[ap])
  if(ap == 1) {ap <- which(inf2== -1)[1]}
  
  bottom <- x[c(dp:ap),] # This is your bottom range
  
  # Have to deal with really shallow dives with no bottom time
  if(bottom$datetime[1] == x$datetime[1]){
     bottom <- data.frame(0)
      
     inf1 <- c(0,0,inf1)
     dp <- which(inf1 == -1)[1] 
      
     inf2 <- c(0,0,inf2)
     ap <- which(inf2 == -1)[1] 
      
  } 
  
  #####Overall Statistics #############
  dout$divenum[i] <- min(x$num)
  dout$start[i] <- x$datetime[1]
  dout$end[i] <-   x$datetime[nrow(x)]
  dout$maxdepth[i] <- max(x$depth)
  dout$Tdown[i]  <- x$sec[nrow(x)]
  dout$Tup[i]  <- x.$sec[nrow(x.)] - x$sec[nrow(x)]
  
    if (nrow(bottom) > 1){
      dout$bottom.duration[i] <- as.numeric(bottom$datetime[nrow(bottom)]) - as.numeric(bottom$datetime[1])
      dout$descent[i] <- as.numeric(bottom$datetime[1]) - as.numeric(x$datetime[1])
      dout$ascent[i]  <- as.numeric(x$datetime[nrow(x)]) - as.numeric(bottom$datetime[nrow(bottom)]) 
      dout$descent.rate[i] <- mean( x[1:which(x$datetime == bottom$datetime[1]),]$depth.diff / fs )
      dout$ascent.rate[i]  <- mean( x[which(x$datetime == bottom$datetime[nrow(bottom)]):nrow(x),]$depth.diff / fs )
    }
    
    if (nrow(bottom) == 1){
      dout$descent[i] <- as.numeric(x$datetime[dp]) - as.numeric(x$datetime[1])
      dout$ascent[i]  <- as.numeric(x$datetime[nrow(x)]) - as.numeric(x$datetime[ap]) 
      dout$descent.rate[i] <- mean( x[1:which(x$datetime == x$datetime[dp]),]$depth.diff / fs )
      dout$ascent.rate[i]  <- mean( x[which(x$datetime == x$datetime[ap]):nrow(x),]$depth.diff / fs )  
    }
  
  ##################wiggles ############
  ##set up empty variables
  diff.d <- NULL
  diff.2 <- NULL
  diff.3 <- NULL
  
  ##find records with a change in direction
  # Loop through each record  1 to the max record
  for(j in 1: nrow(x)) {
    # Take for example, second record from the first record
    diff.d[j+1] <- x$depth[j+1]-x$depth[j] ## gets diff from one record to the next
    
    # Set Diff.2 Variable
    diff.2[j+1] <- ifelse(diff.d[j+1]<0,  -1, diff.d[j+1]) ## sets all negatives to -1
    diff.2[j+1] <- ifelse(diff.2[j+1]>0,  1, diff.2[j+1])  ##sets all positives to +1
    diff.2[j+1] <- ifelse(diff.2[j+1]==0, diff.2[j], diff.2[j+1])
    diff.3[j+1] <- diff.2[j+1]-diff.2[j]  ##relaces 0s with the precding value
  }
  if(length(which(diff.3!=0))< 3) diff.3<- NULL    ##excludes points that are not real wiggles
  wig <- x[which(diff.3!=0)-1,]  ##number of wiggles
  dout$wig.n[i] <- nrow(wig)
  
  ##-------------------Optional Plot of each dive------------------------##
  # Comment this section out to speed up code. 
  # Option to save plots
    
  if(nrow(bottom) > 1){
     #plot.new()  # Comment if you only want to print, not save plots
     #par(mfrow=c(1,1),mar=c(4,4,2,2))  # Comment if you only want to print, not save plots
     #tiff(paste("Dive", bottom$num[1],".tiff"),res=100)  # Comment if you only want to print, not save plots
     tit1 <- paste("Dive", bottom$num[1])
      
     plot(x$datetime, x$depth, ylim=c(max(bottom$depth),0), t="l", ylab="Depth", lwd=2, axes=T, main=tit1,xlab='Time')
     points(bottom$datetime, bottom$depth, t="b", lwd=4, col="red")
     #points(wig$datetime, wig$depth, pch=19, cex=1.5, col="blue")
     abline(v=bottom$datetime[1],col='green',lwd=3)
     abline(v=bottom$datetime[nrow(bottom)],col='green',lwd=3)
     #dev.off()  # Comment if you only want to print, not save plots
    }
    
  if(nrow(bottom) == 1){
     #plot.new()  # Comment if you only want to print, not save plots
     #par(mfrow=c(1,1),mar=c(4,4,2,2))  # Comment if you only want to print, not save plots
     #tiff(paste("Dive", x$num[1],".tiff"),res=100)  # Comment if you only want to print, not save plots
     tit1 <- paste("Dive", x$num[1]) 
      
     plot(x$datetime, x$depth, ylim=c(max(x$depth),0), t="l", ylab="Depth", lwd=2, axes=T, main=tit1,xlab='Time')
     #points(wig$datetime, wig$depth, pch=19, cex=1.5, col="blue")
     abline(v=x$datetime[dp],col='green',lwd=3)
     abline(v=x$datetime[ap],col='green',lwd=3)
     #dev.off()  # Comment if you only want to print, not save plots
    }  
  }


  # We want to return the new ptmp datafrme with the zero corrected depth variable and dive numbers 
  # as well as the dive statistics (dout)
  return(list(dout,ptmp))


}
