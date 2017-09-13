##--------------------------------------------------------------------------------
##
## Final Project: Comparing rRNA depletion methods
##
## Class: PCE Data Science Methods Class
##
## Name: Ben Larson
##
##--------------------------------------------------------------------------------

# Change working directory to project directory
# setwd("/data/basespace/submit/")

## Load Libraries----
library(logging)
library(dplyr)
library(grofit)
library(stringr)
library(boot)
library(pbapply)

## Declare functions----

# Set up a logger that goes to a file but suppresses console output
createLog <- function(logFile, loggingLevel) {
    # Remove the test file (invisibly)
    suppressWarnings(invisible(file.remove(logFile)))
    
    # Reset the logging
    logReset()
    
    # Create the logger
    x <- getLogger()
    
    # Set the level of the logger
    setLevel(loggingLevel, x)
    
    # Create a handler, set it's level, and direct it to the logger
    # NOTE: Change first argument to writeToConsole to show on screen
    addHandler(writeToFile, file = logFile, level = loggingLevel, logger = x)
    
    # Return the logger
    return(x)
}

# Get K-S statistic from two different pile-ups
ks_stat <- function(pup1, pup2, makePlot = FALSE) {
    # Create standard x vector that goes from 1 to the highest base position of 
    # the two vectors
    x_seq = seq(1, max(max(pup1$pos), max(pup2$pos)) ,len=1000)
    
    # Cumuluative distribution function for first pile-up
    pup1.tot <- sum(pup1$tot.count)
    y_cdf1 <- sapply(x_seq, function(base.pos) {
        # NOTE: Can have more than one count at a position, and interesed in
        #       counts here so need to subset the count column for all positions 
        #       less than current position. Since positions are sorted, might be
        #       pup1 more efficient way to do this, but this will work and be more
        #       error proof
        sum(pup1$tot.count[pup1$pos<base.pos])/pup1.tot
    })
    
    # Cumuluative distribution function for second pile-up
    pup2.tot <- sum(pup2$tot.count)
    y_cdf2 <- sapply(x_seq, function(base.pos) {
        sum(pup2$tot.count[pup2$pos<base.pos])/pup2.tot
    })
    
    # Now have distributions standardized over the x distance.
    ks.stat <- max(abs(y_cdf1-y_cdf2))
    
    if (makePlot) {
        # Plot the first cdf with labels
        plot(x_seq,y_cdf1, col='blue', pch=16,
             xlab = "Base Position", ylab = "CDF")
        # Add points for the second cdf
        points(x_seq,y_cdf2,col='red', pch=16)
        # Add a legend
        legend("bottom",
               c(deparse(substitute(pup1)),deparse(substitute(pup2))),
               lty = c(1,1), lwd=c(2.5,2.5), col=c("blue","red"))
        
        # where does max value occur?
        k_index = which.max(abs(y_cdf1-y_cdf2))
        # translate the index to a value
        k_s_x = x_seq[k_index]
        # Add a vertical line to the plot
        lines(c(k_s_x,k_s_x), c(y_cdf1[k_index],y_cdf2[k_index]),
              col='black', lwd=2)
    }
    
    # Return the statistic
    return(ks.stat)
}

simulate_pileup <- function(seed.num, xlocation) {
    # Simulate a pile-up for testing the ks_stat function 
    set.seed(seed.num);
    # Generate random histogram of counts from 1 to 285 (max of a representative 
    # pileup) 4000 times. Remember that lapply sets the anonymous function
    # argument to the first argument in the lapply command, whether this
    # information is needed or not.
    x_seq <- seq(from = 1, to = 4000, by = 25)
    s <- sapply(x_seq, function(x) {sample(x = c(1:285), size = 1)})
    
    # Condition hisotgram to look like a rRNA pile-up
    # NOTE: Selected parameters to give a function that would
    #       cover the range of x-values (base positions) in the data
    #       and span a vertical range from 0 to 1 
    xlogit <- dlogis(x_seq, location = xlocation,
                     scale = 700, log = FALSE)*10^3.44
    s <- ceiling(xlogit*s)
    
    # Show some results
    # plot(x_seq, s)
    # Unset the seed
    set.seed(Sys.time())
    # Construct Result
    data.frame(pos = as.integer(x_seq), tot.count = as.integer(s))
}

# Function for plotting simulation results and calculating a p-value
empirical.p.value <- function(ksh, kst, numSim) {
    # Calculate a p-value for the distribution of observations
    emp <- sum(ksh > kst)/numSim
    return(emp)
}

compare.pups <- function(pup1, pup2, ks.stat.sample) {
    # Calculate empirical p value
    emp <- empirical.p.value(ks_hypothesis$t, ks.stat.sample, num.sim)
    # Make the histogram for the null
    hist(ks_hypothesis$t, main = "K-S Stat Distribution for Exptl. Replicates",
         xlab = "K-S Statistic")
    abline(v=ks.stat.sample, col='red', lwd=2)
    # Log results
    loginfo(paste("Scenario:",
                 deparse(substitute(pup1)), "vs.",
                 deparse(substitute(pup2))))
    loginfo(paste("K-S Stat:", ks.stat.sample))
    if(emp==0) {
        loginfo(paste("Empirical p value < ", 1/num.sim))
    } else {
        loginfo(paste("Empirical p value = ", emp))
    }
}

## Declare unit tests----
test_logger <- function() {
    # Call the logging function to be tested
    logMe <- createLog("unitTest.log", "DEBUG")
    
    # Create a log entry
    logdebug("Result correctly logged")
    
    # Make sure the entry was correctly logged
    # NOTE: both dplyr and stringr libraries must be loaded for this to work
    stopifnot(readLines("unitTest.log") %>% str_detect("Result correctly logged"))
    
    # Reset the logger
    logReset()
    
    # Remove the test file (invisibly)
    invisible(file.remove("unitTest.log"))
}

test_ks_stat <- function() {
    # Simulate first pileup, cenetered around position 1000
    pup1 <- simulate_pileup(315, 1000)
    # Simulate second pileup, cenetered around position 1700
    pup2 <- simulate_pileup(513, 1700)
    # Return the K-S statistic
    ks.stat <- round(ks_stat(pup1, pup2, makePlot = FALSE), digits = 3)
    # Check the result
    stopifnot(ks.stat==0.203)
}

## Main Script----
if(interactive()) {
    # Conduct unit tests
    test_logger()
    test_ks_stat()
    
    
    #Create Logger (overwrites exsiting files)
    log.me <- createLog("basespace.local.log", "DEBUG")
    
    logdebug("***ATTENTION***")
    logdebug(paste("You ***MUST*** be in the directory where larson.Rdata was",
                   "downloaded"))
    
    tryCatch(
        {
            # Load the data generated on AWS and downloaded locally, should
            # contain a list of data frame pile-ups
            load("larson.Rdata")
        },
        warning = function(cond) {
            cat(paste0("Could not find file:","\"","larson.Rdata","\"\n"))
            cat("Make sure it was downloaded from remote host")
        }
    )
    
    loginfo("Data sucessfully loaded")
    
    # Clean trailing underscore in names of list elements
    names(pup) <- lapply(names(pup),
                         function(x) {str_sub(x, start = 1, end = -2)})
    
    # Show a pile-up
    plot(pup$L11_L1$pos, pup$L11_L1$tot.count, 
         xlab = "Count", ylab = "Position",
         main =  "Mapped Read Pile-up")
        
    # Generate K-S statistics for NULL distribution
    # General Idea, compare all lanes within any single run, where method was
    # constant and pile-ups should be similar
    # Get 1 combination each for L11, L12, and L14 (which only had 2 lanes)
    # Get 6 combinations each for Red07 and Red08 (which each had 4 lanes)
    # Total ks.stats for NULL is 15,
    
    # First create a list with two elements containing the indices of the pile-up
    # data frames that can be compared with one another
    index.null <- list(i1 = c(1, 3, 11,
                              combn(c(4, 5, 6, 13), 2)[1,],
                              combn(c(7, 8, 9, 14), 2)[1,]),
                       i2 = c(2, 10, 12,
                              combn(c(4, 5, 6, 13), 2)[2,],
                              combn(c(7, 8, 9, 14), 2)[2,]))
    
    # Apply ks_stat to combinations qualified for null to get a base set of 15
    ks.stat.null <- mapply(function(a, b) {ks_stat(pup[[a]], pup[[b]])},
                           index.null$i1, index.null$i2)
    
    # Bootstrap the set of 15
    # Create boot function
    mean_trim_fun = function(data, ix) mean(data[ix])
    num.sim <- 5000
    ks_hypothesis = boot(ks.stat.null, mean_trim_fun, R=num.sim, sim="ordinary")
    # Show some results
    hist(ks_hypothesis$t)
    
    # Get K-S stats for depleted samples relative to a control (e.g., L11_L1)
    # Order of depletion is L12 -> Red7 -> Red8 (See Table)
    
    loginfo(paste("Percent RNA shown in following log entries comes from",
                  "*_Log.final.out stored on the remote host in",
                  "/data/basespace/star_output"))

    
    # Compare lanes in an experiment
    # NOTE: ks_stat and compare.pups not combined to make plot labeling easier
    # L11_L1 v L11_L2
    ks.stat.sample <- ks_stat(pup$L11_L1, pup$L11_L2, makePlot = TRUE)
    compare.pups(pup$L11_L1, pup$L11_L2, ks.stat.sample)
    
    # Compare 3 different depleted vs. the undepleted control
    # L11_L1 v L12_L1
    ks.stat.sample <- ks_stat(pup$L11_L1, pup$L12_L1, makePlot = TRUE)
    compare.pups(pup$L11_L1, pup$L12_L1, ks.stat.sample)
    
    # L11_L1 v Red7_L1
    ks.stat.sample < ks_stat(pup$L11_L1, pup$Red7_L1, makePlot = TRUE)
    compare.pups(pup$L11_L1, pup$Red7_L1, ks.stat.sample)
    
    # L11_L1 v Red8_L1
    ks.stat.sample <- ks_stat(pup$L11_L1, pup$Red8_L1, makePlot = TRUE)
    compare.pups(pup$L11_L1, pup$Red8_L1, ks.stat.sample)
}