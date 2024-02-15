# RRR

# use fft- fast fourrier transform
# to ID if there's an annual oscillation

# toy data

my.x <- 0:999
my.y <- sin(2*pi*50/1000 * (my.x)) # 1000 / 50 = period of 20
plot(my.x, my.y, type = "l")

# check is period as expected? yes peaks are 20 apart
plot(my.x, my.y, type = "l", xlim = c(0,100))
abline(v = 40, col = "Red")
abline(v = 20, col = "blue")

# do fft
my.fft <- fft(z = my.y)
my.fft <- abs(my.fft) # always do this, gets rid of the i's, it's just how you do it.

# check is the peak where we expect?
plot(my.fft, type = "l")
# note there's always 2 peaks because it's mirrored, so just look at the first half
index.half <- 1:(length(my.fft) / 2)
plot(my.fft[index.half], type = "l")
index.peak <- which(my.fft[index.half] %in% max(my.fft[index.half]))
index.peak # 51
my.fft[index.peak] # 500
my.period <- length(my.fft) / (index.peak - 1)
my.period # 20 as expected


# add gradual change ----
my.line <- .01*my.x + 30 # y = mx + b
my.y <- my.y + my.line
plot(my.y ~ my.x, type = "l")

# fit the line
my.lm <- lm(my.y ~ my.x)
abline(my.lm)
my.lm <- summary(my.lm)
my.slope <- my.lm$coefficients[2,1]
my.intercept <- my.lm$coefficients[1,1]

# subtract the line from the data (also centers it at zero)
my.fit <- my.slope * my.x + my.intercept
my.y <- my.y - my.fit
plot(my.y, type = "l")

# do fft
my.fft <- fft(z = my.y)
my.fft <- abs(my.fft) # always do this, gets rid of the i's, it's just how you do it.

# check is the peak where we expect?
plot(my.fft, type = "l")
# note there's always 2 peaks because it's mirrored, so just look at the first half
index.half <- 1:(length(my.fft) / 2)
plot(my.fft[index.half], type = "l")
index.peak <- which(my.fft[index.half] %in% max(my.fft[index.half]))
index.peak # 51
my.fft[index.peak] # 499.8804 close
my.period <- length(my.fft) / (index.peak - 1)
my.period # 20 as expected


# how to calculate the index of your expected period
expected.period <- 20
# period = length of vector / index of spike
expected.index <- length(my.fft) / expected.period + 1 # plus 1 because indexing is not zero-based
expected.index # 51 as expected


# use periodogram test instead
library(TSA)

my.data <- periodogram(my.y, log = "yes", plot = T)

# ok this package is confusing
