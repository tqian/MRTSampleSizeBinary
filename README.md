# mrtbincalc

This package provides a sample size calculator for Micro-randomized Trials with binary outcomes. The sample size formula is developed in Sample Size Considerations for Micro-Randomized Trials with Binary Outcome 


To use this package you need some required inputs, listed below under subtiles numbered from one to four.
1. Study Setup, which includes the duration of the study in days, the number of decision time points within each day, the randomization probability, i.e. the probability of assigning the treatment at a decision time point. 

2. Availability
Treatment can only be provided when an individual is available; The expected availability is the probability a person is available to receive the intervention at the decision times.
You need to select a time-varying pattern for the expected availability. There are tow patterns you can choose from: constant, or linear over decision points.

3. Success Probability Null Curve
The Success Probability Null Curve at each decision time point is defined as the probability of the proximal outcome equal to 1 for available individuals who are not assigned treatment.
You need to provide the trend of success probability null curve. This could be either constant or linear over decision points.

4. Proximal Treatment Effect
The Proximal Treatment Effect at each decision time point is defined as the mean difference in the proximal outcome between available people who are assigned a treatment versus available people who are not assigned treatment. In this work, we only consider the binary treatment. You need to provide the trend of proximal treatment effects. This could be either constant or time-varying, e.g. linear or quadratic over decision points.

The outputs one can get using this package are mainly two values, power and sample size.
In both cases, you will need to input the desired significance level.
