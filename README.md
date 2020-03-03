# Quantile extraction from compound distributions

The program can be used to compute any quantiles from a compound distribution where there is no closed form formulae to describe the distribution. 

A randomly stopped sums or a compound distribution is the distribution of the sum of independently distributed random variables, where the number of such variables is also random.
As an example, in a auto insurance portfolio, the occurence of car accidents is random and can be described as a Poisson distribution. If a random number of accidents occurs, the interest is in the sum of the losses from each individual accidents. 
The sum then is a convolution and described as a compound distribution. 

Fourier transform is a better alternative to the simulation approach in terms of computation time. It is also a good alternative to Panjer recursion. 

The method is applicable to property and casualty insurance, pyschology studies, and operational risk management where aggregate/compound distributions are used. 

In these areas, researchers are interested in either mixed effect regressions or distributions without covariates.  

## References
+ *Stuart A. Klugman, Gordon E. Willmot, Harry Panjer. Loss Models: From Data to Decisions. 3rd ed, 1998. §9.3*
+ *Michael Smithson and Yiyun Shou. Randomly stopped sums: models and psychological applications. Front. Psychol., 10 November 2014*
+ *https://arxiv.org/pdf/1211.5802.pdf*
