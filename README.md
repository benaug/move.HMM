move.HMM
========

Flexible, user-friendly Hidden (Semi) Markov Models for animal movement data (Langrock et al. 2012)

This package is intended to take the maximum likelihood methodology
for HMM and HSMM animal movement models of Langrock et al. (2012) 
and make it accessible to ecologists.  Further, it works for an 
arbitrary number of hidden states, arbitrary number of observation
distributions, and arbitrary observation distributions.  This aids
model comparison via AIC and there are several functions for graphical
goodness of fit checking.

The code is also now more efficient.

Before installing the package, install it's dependencies:
numDeriv
pscl
CircStats
psych
VGAM

If bootstrapping, install
foreach
snow
doSNOW


Things to add:

1.  Profile likelihood CIs?
2.  Random effects, maybe via ADMB
3.  covariate-dependent state transitions
4.  mid-residuals for discrete distributions

