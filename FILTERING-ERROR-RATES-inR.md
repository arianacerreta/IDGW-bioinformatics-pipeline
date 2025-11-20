# Final filtering and error rates in R

## 1. Download and install [SNPfiltR](https://devonderaad.github.io/SNPfiltR/) (DeRaad 2022)
1. You can use this [install_SNPfiltR.R](utility_files/R_code/install_SNPfiltR.R) code to help install this package. You can also refer to the package's github linked above.

## 2. Filter your data using read values from PCRneg for each run
To account for potentially contamination in the lab, we filter usuable loci my comparing the read depth of a given locus for an individual to the read depth found in the negative control for each plate. We only keep loci for an individual that have >2x read depth compared to the PCR negative. This is often not too important.
1. Download [filter_with_PCRneg.R](utility_files/R_code/filter_with_PCRneg.R). Requires ```vcfR```, ```tidyverse```, ```ggplot2```, and ```SNPfiltR```.
2. Edit path in Line 8 to your file
3. Lines 9-10 will be part of your naming convention in the final line of code. Adjust as you see fit.
4. Lines 11-19, 31-35, and 59 deal with our in-house naming conventions. You may need to adjust for your case use.
5. Edit path in line 69 for your outputs.

## 3. Calculate error rates by filtering scheme
1. Download [error_rates_GWAdapt.R](utility_files/R_code/error_rates_GWAdapt.R). Requires ```vcfR```, ```tidyverse```, ```ggplot2```, and ```SNPfiltR```.
2. Edit line 8 with your path. Edit line 11 with the name of your panel.
