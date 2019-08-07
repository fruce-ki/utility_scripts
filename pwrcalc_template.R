library(pwr)

## https://www.r-bloggers.com/validating-type-i-and-ii-errors-in-a-b-tests-in-r/

# How many measurements do I need in each group (of equal size), to measure a given difference between treatments A and B.

cr_a <- 0.25     # the expected success rate for group A
mde <- 0.1       # minimum detectable effect (%difference between success rates A and B)
alpha <- 0.05    # the false positive rate
power <- 0.80    # 1-false negative rate (beta)

ptpt <- pwr.2p.test(h = ES.h(p1 = cr_a, p2 = (1+mde)*cr_a), 
                    sig.level = alpha, 
                    power = power
)
n_obs <- ceiling(ptpt$n)

## !!! The p-value of a test does not decrease monotonicaly with sample size increase !!!
## At smaller sample sizes it can fluctuate randomly, with high probability of mislassification (much higher than the specified one).
## p-value stabilises as the sample size meets and exceeds the size determined by the power calculation.
