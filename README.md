# safegraph_DiD_piecewise_likelihood_MLE
This is a basic Julia script for difference-in-differences estimation that incorporates the unique censoring and truncation found in publically available Safegraph patterns data


This is not a plug-and-play function, the script itself must be edited for the particular needs of your project. This is the first script I've ever written in Julia, so I apologize for clunkiness and the less than ideal user experience. 

The code requires an n by 4 .rda file (dataframe) as an input. column 1 is the outcome variable, column 2 is the "post" indicator. column 3 is the "treatment group" indicator, and column 4 is the "post by treated" indicator. 

It estimates the coefficient values based on a custom piecewise likelihood function that accounts for the censoring and truncation of the data generating process for Safegraph patterns data. This is discussed [here](https://docs.safegraph.com/docs/monthly-patterns).

To quote that page: " We do not report data unless at least 2 visitors are observed from that group. If there are between 2 and 4 visitors this is reported as 4. "

The code does not necessarily need to be for DiD estimation, the likelihood function is essentially the same (my purposes were DiD, and I do not know enough Julia to make this more adaptable)

It will also bootstrap vectors of coefficient estimates, and take their standard deviatons for estimates of standard errors. 

I hope this is helpful in some capacity.
