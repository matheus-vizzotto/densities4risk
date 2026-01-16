set.seed(123) # for reproducibility
data <- rt(n = 100, df = 3) 
plot(data)

library(kdensity)

gausst <- list(
  kernel = function(y, x, h) {
    stats::dt((y - x) / h, df = 3)
  },
  sd = sqrt(5 / (5 - 2)),          # finite variance, df > 2
  support = c(-Inf, Inf)
)
# bandwidth choice (example Silverman or custom)
h <- 1.06 * sd(data) * length(data)^(-1/5)

kde_t <- kdensity(
  data,
  kernel = gausst,
  bw = h
)

kde_n <- kdensity(
  data,
  bw = h
)

plot(kde_n)
lines(kde_t, col='red')
