# Bayesian adaptive psychometric estimation #

IMPORTANT NOTE:
This work is based partially on previous code, acquired informally through the academic grapevine. Attributions are forthcoming.

The basic estimation procedure in this code follows Kontsevich and Tyler, 1999, the citation for which is (APA):

Kontsevich, L. L., & Tyler, C. W. (1999). Bayesian adaptive estimation of psychometric slope and threshold. Vision Research, 39(16), 2729–2737. https://doi.org/10.1016/S0042-6989(98)00285-5

This procedure has been modified to incorporate a stickiness parameter, which re-weights the "Bayesian" update to be more conservative than a than a truly Bayesian scheme, by increasing the contribution of the prior and decreasing the contribution of the posterior (relative to Bayes rule) on each step.

Testing and parameter recovery code, as well as code to manually and randomly generate "mouse" multi-session behavioral (incl. learning) profiles have also been introduced.

## The code ##

- The adaptive estimation code is in `PsyRoutine.m`
- The script `testPsiRoutine.m` illustrates how to use the PsiRoutine class and its components
- The functions `get[Basic | Random]MouseMeta()` generate test "mice" for parameter recovery analysis.

## Other notes ##
The Weibull CDF is given by

![equation](https://latex.codecogs.com/svg.latex?p%20%3D%20%281-l-g%29%281-%5Ctext%7Bexp%7D%5E%7B-%5Cfrac%7Bx%7D%7Bt%7D%5Es%7D%29%29%20&plus;%20g)

so that its slope at threshold is

![equation](https://latex.codecogs.com/svg.latex?%5Cfrac%7Bdp%7D%7Bdx%7D%7C_%7Bx%7D%20%3D%20%281-l-g%29%5Cfrac%7Bs%7D%7Bt%7De%5E%7B-1%7D)

The parameter s is estimated by the adaptive algorithm, and is called the "slope" parameter (although it is technically the "shape" parameter).

The parameter recovery code seems to indicate that with moderate hardcoded guess and lapse rates, the model is fairly robust to mis-specification.

Additional sampling modifications, e.g. to avoid directional bias in sampling, introduce catch trials, etc will be included shortly.
