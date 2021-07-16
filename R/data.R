#' Life expectancy in the US (2020)
#' 
#' A simulated data set based on statistics of life expectancy in the US at birth.
#' 
#' @details The data was generated using official statistics from the CDC and a
#' study of life expectancy of smokers in the US published in the New England
#' Journal of Medicine (see references).
#' 
#' According to the CDC, data from 2020 indicates that the average life expectancy
#' of females in the US is 80.5 years vs 75.1 years for males (which declined
#' with respect to 2019 after COVID hit the US). In Jha et al. (2013), evidence
#' is presented indicating that individuals who smoke have at least ten years
#' left of life expectancy compared to non-smokers.
#' 
#' The parameter estimates for the data generating process where:
#' 
#' - An average of life expectancy of 80.1.
#' - Smokers live 10 years less than non-smokers.
#' - Females live 5.4 years longer than males.
#' 
#' @format 
#' A 1,000 rows data frame with three columns:
#' 
#' - `age`: Years lived.
#' - `smoke`: 0/1 variable equal to 1 if the individual smokes.
#' - `female`: 0/1 variable equal to 1 if the individual is a female at birth.
#' 
#' @references 
#' Jha, P., Ramasundarahettige, C., Landsman, V., Rostron, B., Thun, M., Anderson, R. N.,
#' ... Peto, R. (2013). 21st-Century Hazards of Smoking and Benefits of Cessation
#' in the United States. New England Journal of Medicine, 368(4), 341–350.
#' \doi{10.1056/NEJMsa1211128}
#' 
#' Arias, E., Tejada-Vera, B., & Ahmad, F. (2021). Provisional life expectancy
#' estimates for January through June, 2020.
#' \url{https://www.cdc.gov/nchs/data/vsrr/VSRR10-508.pdf}
"lifeexpect"