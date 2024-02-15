# RRR
# try to use the strucchange package to find a breakpoint

library(data.table)
library(strucchange)

broke <- data.table("x" = 1:70, "y" = c(rnorm(n = 20, mean = 5, sd = 1), rnorm(n = 30, mean = 10, sd = 1), rnorm(n = 20, mean = 5, sd = 1)))
plot(y ~ x, data = broke)

?breakpoints()
breakpoints(broke$y ~ 1)

fs <- Fstats(broke$y ~ 1)
plot(fs)
fs$breakpoint
abline(v = fs$breakpoint / fs$nobs, col = "green") # calls peak whether signif or not
p.val <- sctest(fs)
p.val$p.value
breakpoints(fs)
lines(breakpoints(fs))

# # notes on what's what:
# plot(my.fstats)
# lines(boundary(my.fstats), col = "green")
# lines(my.fstats$Fstats, col = "green")
# abline(v = my.fstats$breakpoint / my.fstats$nobs, col = "green") # my.fstats$breakpoint is the index
# plot(linear.interp.dist.first$dist ~ linear.interp.dist.first$time.approx.year, type = "l")
# abline(v = linear.interp.dist.first$time.approx.year[my.fstats$breakpoint], col = "green")

# so it ID'ed the index of the LAST point BEFORE the break
# F test is only testing for a single breakpoint
# but breakpoints


# ok ok this is gonna be messy, but jsut paste in the efp stuff and the breakpoints stuff here, went with fstat instead in the real script

breakpoints.obj <- get.breakpoints.strucchange(linear.interp.10day =  linear.interp.10day)
breakpoints.table <- get.breakpoints.table(package.output = breakpoints.obj, linear.interp.10day = linear.interp.10day)
step.change.stats <- classify.if.step.change(break.tab = breakpoints.table)
step.change.stats <- adjust.breakpoint.dates(linear.interp.10day = linear.interp.10day, change.stats = step.change.stats)
step.change.stats <- get.step.change.date(change.stats = step.change.stats, first.date = dist.first$date.1[1])

# for paper: I identified significant breakpoints using a recursive cumulative sum fluctuation process 
# (efp function in the strucchange R package) and identified breakpoints with a confidence level of alpha = 0.05. 
# I calculated the length of the change based on the length of time until the confidence level was crossed again.
# Since any cumulative sum-based changepoint detector has a lag, I identified the true breakpoint location by 
# choosing the nearest breakpoint identified by a least squares fit with a Bayesian information criterion penalty 
# (breakpoint function in strucchange R package). In this way I used efp() to identify statistically significant breakpoints,
# and breakpoint() to identify the exact break location. If breakpoints were not called by both functions, I considered them absent. 
# If a breakpoint resulted in change that lasted for at least XX years (i.e., no additional breakpoints in that timeframe) 
# and occurred after at least 1 year (i.e. not counting any breakpoints right at the start) I considered it a step change pattern. 

get.breakpoints.strucchange <- function(linear.interp.10day){
  # efp gives a confidence interval, and when the line crosses the confidence interval it's a significant breakpoint
  # z <- efp(linear.interp.10day$dist ~ linear.interp.10day$time.approx.year) # this allows slanted lines
  z <- efp(linear.interp.10day$dist ~ 1) # this requires flat lines with zero slope
  # plot(z) # this is the package's plot, I used the default test which is
  return(z)
}

get.breakpoints.table <- function(package.output, linear.interp.10day){
  # y <- data.table("abs.efp" = abs(package.output$process), "conf.bound" = boundary(package.output), "x.vals" = linear.interp.10day$time.approx.year[-1]) # need the [-1] if fitting sloped lines
  y <- data.table("abs.efp" = abs(package.output$process), "conf.bound" = boundary(package.output), "x.vals" = linear.interp.10day$time.approx.year) # not if fitting flat lines
  y[ ,above.cutoff := abs.efp >= conf.bound]
  y[ ,above.cutoff := as.numeric(above.cutoff)]
  y$change.point <- FALSE
  for (r in 2:nrow(y)){
    if (y[r, above.cutoff] != y[r - 1, above.cutoff]){
      y[r, change.point := TRUE]
    }
  }
  
  return(y)
}

classify.if.step.change <- function(break.tab){
  
  # call the last data point a change point as well if it's above the cutoff to calculate the length of change
  if (break.tab[nrow(break.tab), above.cutoff] == 1){
    break.tab[nrow(break.tab), change.point := TRUE]
  }
  
  # get maximum time above the red confidence line, if it's ever above it
  num.change.points <- nrow(break.tab[change.point == TRUE])
  if (num.change.points > 0){
    break.tab.changes <- break.tab[change.point == TRUE, x.vals]
    
    # remove breeakpoints called right at the beginning, say it's got to be at least a year in to call a breakpoint
    index.remove <- which(break.tab.changes < 1)
    if (length(index.remove) %% 2 != 0){ # if it's odd
      index.remove <- c(index.remove, max(index.remove) + 1)
      break.tab.changes <- break.tab.changes[-index.remove]
    }
    num.change.points <- length(break.tab.changes)
    
    if (num.change.points > 0){
      max.length <- NA
      change.loc <- break.tab.changes[1]
      for (n in seq.int(from = 1, to = (num.change.points - 1), by = 2)){ # step by 2 to only tally chunks ABOVE the red line
        length.change <- break.tab.changes[n + 1] - break.tab.changes[n]
        if (length.change > max.length & !is.na(max.length)){
          change.loc <- break.tab.changes[n]
        }
        max.length <- max(max.length, length.change, na.rm = T)
      }
    }else{
      max.length <- NA
      change.loc <- NA
    }
  }else{
    max.length <- NA
    change.loc <- NA
  }
  
  return(data.table("step.change.length" = max.length, "step.change.loc" = change.loc))
}

adjust.breakpoint.dates <- function(linear.interp.10day, change.stats){
  if(!is.na(change.stats$step.change.length)){
    # my.obj <- breakpoints(linear.interp.10day$dist ~ linear.interp.10day$time.approx.year)
    my.obj <- breakpoints(linear.interp.10day$dist ~ 1)
    my.tab <- linear.interp.10day[my.obj$breakpoints]
    
    break.start <- my.tab[time.approx.year <= change.stats$step.change.loc, time.approx.year]
    if (length(break.start) == 0){
      break.start <- NA
    }else{
      break.start <- break.start[which((change.stats$step.change.loc - break.start) == min((change.stats$step.change.loc - break.start)))]
      cat("adjusting breakpoint start from",change.stats$step.change.loc,"to",break.start,"\n")
    }
  }else{
    break.start <- NA
  }
  change.stats <- change.stats[ ,breakpoint.loc := break.start]
  return(change.stats)
}

get.step.change.date <- function(change.stats, first.date){
  start.year <- parse_date_time(first.date, "ymd") |>
    year()
  
  if (!is.na(change.stats[ ,step.change.loc])){
    change.date <- date_decimal(change.stats[ ,step.change.loc] + start.year) |>
      round_date(unit = "day") |>
      as.character()
    change.stats[ ,step.change.date := change.date]
  }else{
    change.stats[ ,step.change.date := NA]
  }
  if (!is.na(change.stats[ ,breakpoint.loc])){
    change.date <- date_decimal(change.stats[ ,breakpoint.loc] + start.year) |>
      round_date(unit = "day") |>
      as.character()
    change.stats[ ,breakpoint.date := change.date]
  }else{
    change.stats[ ,breakpoint.date := NA]
  }
  
  return(change.stats)
}