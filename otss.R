##------------------ Subsample OpenTag sensor data----------##

# Use this function to subsample OpenTag sensor data (Loggerhead
# Instruments, loggerhead.com). For example, you can use this 
# function to subsample your raw data collected at 100 Hz to 5 Hz. 

# Reny Tyson, rtyson@mote.org
# 5 October 2016

##############################################################
otss <- function(d,w,order,Hz,ss){

require(signal)
  
# d is your raw data  
# W should be at least half of your Nygest frequency
# order is the desired order of the butterworth filter
# Hz is the sampling frequency of your raw data
# ss is your desired sampling frequency  

#Filter your data
bf <- butter(order,w,type='low')
filt <- filter(bf,d)
s.s <- seq(1,length(d),Hz/ss)

#subset your data
ssd <- d[s.s]

return(ssd)
}
