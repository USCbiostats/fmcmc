## code to prepare `lifeexpect` dataset goes here

set.seed(123)
n <- 1000
lifeexpect <- data.frame(
  smoke  = sample.int(2, n, replace = TRUE) - 1L,
  female = sample.int(2, n, replace = TRUE) - 1L
)

# According to 
# Arias, E., Tejada-Vera, B., & Ahmad, F. (2021). Provisional life expectancy estimates for January through June, 2020.
# url https://www.cdc.gov/nchs/data/vsrr/VSRR10-508.pdf
# Difference in 2020 in the US was of 5.4 years.
#
# Smokers have < 10 years of life expectancy
# according to:
# Jha P, Ramasundarahettige C, Landsman V, Rostrom B, Thun M, Anderson RN, McAfee T, Peto R. 21st Century
# Hazards of Smoking and Benefits of Cessation in the United States    [PDF–738 KB]external icon.
# New England Journal of Medicine, 2013;368(4):341–50 [accessed 2015 Aug 17].
females <- 80.5
males   <- 75.1

lifeexpect$age <- with(lifeexpect,  -10 * smoke + (females - males) * female + rnorm(n, males + 5, sd = 2)) 
MASS::truehist(lifeexpect$age[lifeexpect$female == 1])

summary(lm(age ~ smoke + female, data = lifeexpect))

usethis::use_data(lifeexpect, overwrite = TRUE, compress = "gzip")
