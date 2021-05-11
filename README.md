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
The slope estimation either has a bug, or is just very poorly constrained. The code as I recieved it had a severe bug in that component, which I fixed, so it is likely there is some additional error there. This doesn't appear to affect parameter recovery of threshold, so I haven't attempted to fix it yet.

The parameter recovery code seems to indicate that with moderate hardcoded guess and lapse rates, the model is fairly robust to mis-specification.

Additional sampling modifications, e.g. to avoid directional bias in sampling, introduce catch trials, etc will be included shortly.
