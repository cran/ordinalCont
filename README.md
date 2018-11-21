Ordinal regression analysis for continuous scales
==================

Ordinal regression analysis is a convenient tool for analyzing ordinal response variables in the presence of covariates. We extend this methodology to the case of continuous self-rating scales such as the Visual Analog Scale (VAS) used in pain assessment, or the Linear Analog Self-Assessment (LASA) scales in quality of life studies. Subjects are typically given a linear scale of 100 mm and asked to put a mark where they perceive themselves. These scales measure subjects' perception of an intangible quantity, and cannot be handled as ratio variables because of their inherent nonlinearity. Instead we treat them as ordinal variables, measured on a continuous scale. We express the likelihood in terms of a function (the “g function”) connecting the scale with an underlying continuous latent variable. In the current version the g function is expressed with monotone increasing I-splines (Ramsey 1988). The link function is the inverse of the CDF of the assumed underlying distribution of the latent variable. Currently the logit link, which corresponds to a standard logistic distribution, is implemented. (This implies a proportional odds model.) The likelihood is maximized using the MI algorithm (Ma, 2010). Fixed- and mixed-effects models are implemented in the function *ocm*.

### References

- Manuguerra M, Heller GZ (2010). Ordinal Regression Models for Continuous Scales, *The International Journal of Biostatistics*: 6(1), Article 14.

- Heller, GZ, Manuguerra M, Chow R (2016). How to analyze the Visual Analogue Scale: Myths, truths and clinical relevance, *Scandinavian Journal of Pain*, Volume 13, 67 - 75
- Ma, J. (2010). Positively Constrained Multiplicative Iterative Algorithm for Maximum Penalized Likelihood Tomographic Reconstruction, *Nuclear Science* 57 (1): 181-92.

- Ramsay, J. O. (1988). Monotone regression splines in action. *Statistical science*, 425-441.

