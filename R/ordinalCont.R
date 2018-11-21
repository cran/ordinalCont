#' ordinalCont-package
#'
#' @name ordinalCont-package
#' @docType package
#' @details Ordinal regression analysis is a convenient tool for analyzing ordinal response variables 
#' in the presence of covariates. We extend this methodology to the case of continuous self-rating 
#' scales such as the Visual Analog Scale (VAS) used in pain assessment, or the Linear Analog 
#' Self-Assessment (LASA) scales in quality of life studies. Subjects are
#' typically given a linear scale of 100 mm and asked to put a mark where they perceive
#' themselves. These scales  measure subjects' 
#' perception of an intangible quantity, and cannot be handled as ratio variables because of their 
#' inherent nonlinearity.  Instead we treat them as ordinal variables, measured on a continuous scale. We 
#' express  the likelihood in terms of a function (the ``g function'') connecting the  
#' scale with an underlying continuous latent  variable. In the current version the g function 
#' is expressed with monotone increasing I-splines (Ramsey 1988).
#' The link function is the inverse of the CDF of the assumed underlying distribution of the 
#' latent variable. Currently 
#' the logit link, which corresponds to a standard logistic distribution, is implemented. 
#' (This implies a proportional odds model.)  The likelihood is 
#' maximized using the \code{MI} algorithm (Ma, 2010). Fixed- and mixed-effects models are implemented 
#' in the function \code{\link{ocm}}. 
#' @references Manuguerra M, Heller GZ, Ma J (2017). Semi-parametric Ordinal Regression Models for Continuous 
#'  Scales, \emph{Proceedings of the 32nd International Workshop on Statistical Modelling}. July 3-7, 2017, Groningen, Netherlands.
#' @references   Manuguerra M, Heller GZ (2010). Ordinal Regression Models for Continuous 
#'  Scales, \emph{The International Journal of Biostatistics}: 6(1), Article 14.
#' @references Heller, GZ, Manuguerra M, Chow R (2016). How to analyze the Visual Analogue Scale: 
#' Myths, truths and clinical relevance, \emph{Scandinavian Journal of Pain}, Volume 13, 67 - 75
#' @references Ma, J. (2010). Positively Constrained Multiplicative Iterative Algorithm for Maximum 
#' Penalized Likelihood Tomographic Reconstruction, \emph{Nuclear Science} 57 (1): 181-92.
#' @references Ramsay, J. O. (1988). Monotone regression splines in action. \emph{Statistical science}, 425-441.
#' @author Maurizio Manuguerra, Gillian Heller
NULL


#' @title ANZ0001 trial
#' @description The complete ANZ0001 trial data set
#' @details  The ANZ0001 trial, conducted by the ANZ Breast Cancer Trials Group, is an unblinded, multi-centre, 
#' randomized trial with three chemotherapy treatment arms, concluded in 2005 (Stockler et al 2007). 
#' Health-related quality of life measures (Overall quality of life, Physical Well-Being, Mood, Pain, Nausea and 
#' Vomiting, Appetite) are assessed at each chemotherapy treatment cycle, from randomization until 
#' disease progression, when treatment is interrupted. 
#' The treatments Intermittent Capecitabine (IC) and Continuous Capecitabine (CC) are compared with the 
#' standard combination treatment CMF, each with its own protocol. 
#' There is no maximum duration of treatment, but it is interrupted on disease progression, or when patient 
#' intolerance or unacceptable toxicity are recorded.
#' The data set is extracted from the ANZ0001 trial and contains information from 292 patients with complete 
#' quality of life measurements.
#'
#' The variables are as follows:
#'
#' \tabular{ll}{
#' \code{randno}\tab patient ID number \cr
#' \code{cycleno}\tab chemotherapy cycle number \cr
#' \code{age}\tab age of patient at entry to study \cr
#' \code{bsa}\tab Body Surface Area (m\eqn{^2}) \cr
#' \code{treatment}\tab  treatment received by  patient (1,2,3)\cr
#' \code{overall}\tab Overall quality of life as recorded by the patient on a LASA scale,  normalized to  (0, 1)\cr
#' \code{phys}\tab Physical Well-Being as recorded by the patient on a LASA scale,  normalized to  (0, 1)\cr
#' \code{mood}\tab Mood as recorded by the patient on a LASA scale,  normalized to  (0, 1)\cr
#' \code{pain}\tab Pain as recorded by the patient on a LASA scale,  normalized to  (0, 1)\cr
#' \code{nausvom}\tab Nausea and Vomiting as recorded by the patient on a LASA scale,  normalized to  (0, 1)\cr
#' \code{appetite}\tab Appetite as recorded by the patient on a LASA scale,  normalized to  (0, 1)
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ANZ0001
#' @usage data(ANZ0001)
#' @references Stockler, M., T. Sourjina, P. Grimison, V. Gebski, M. Byrne, V. Harvey, P. Francis et al. 
#' ``A randomized trial of capecitabine (C) given intermittently (IC) rather than continuously (CC) compared 
#' to classical CMF as first-line chemotherapy for advanced breast cancer (ABC).'' In \emph{ASCO Annual Meeting 
#' Proceedings}, vol. 25, no. 18_suppl, p. 1031. 2007.
#' @format A data frame with 2473 rows and 11 variables
NULL

#' @title ANZ0001 trial subset
#' @description A subset from the ANZ0001 trial data set
#' @details  The ANZ0001 trial, conducted by the ANZ Breast Cancer Trials Group, is an unblinded, multi-centre, 
#' randomized trial with three chemotherapy treatment arms, concluded in 2005 (Stockler et al 2007). 
#' Health-related quality of life measures (Overall quality of life, Physical Well-Being, Mood, Pain, Nausea and 
#' Vomiting, Appetite) are assessed at each chemotherapy treatment cycle, from randomization until disease 
#' progression, when treatment is interrupted. 
#' The treatments Intermittent Capecitabine (IC) and Continuous Capecitabine (CC) are compared with the 
#' standard combination treatment CMF, each with its own protocol. 
#' There is no maximum duration of treatment, but it is interrupted on disease progression, or when 
#' patient intolerance or unacceptable toxicity are recorded.
#' The data set is extracted from the ANZ0001 trial and contains information from a subset of 292 patients 
#' with complete quality of life measurements, limited to cycle numbers 0 and 5.
#'
#' The variables are as follows:
#'
#' \tabular{ll}{
#' \code{randno}\tab patient ID number\cr
#' \code{cycleno}\tab chemotherapy cycle number, either 0 (initial assessment) or 1 (fifth cycle). \cr
#' \code{age}\tab age of patient at entry to study \cr
#' \code{bsa}\tab Body Surface Area (m\eqn{^2}) \cr
#' \code{treatment}\tab  treatment received by  patient (1,2,3)\cr
#' \code{overall}\tab Overall quality of life as recorded by the patient on a LASA scale,  normalized to  (0, 1)\cr
#' \code{phys}\tab Physical Well-Being as recorded by the patient on a LASA scale,  normalized to  (0, 1)\cr
#' \code{mood}\tab Mood as recorded by the patient on a LASA scale,  normalized to  (0, 1)\cr
#' \code{pain}\tab Pain as recorded by the patient on a LASA scale,  normalized to  (0, 1)\cr
#' \code{nausvom}\tab Nausea and Vomiting as recorded by the patient on a LASA scale,  normalized to  (0, 1)\cr
#' \code{appetite}\tab Appetite as recorded by the patient on a LASA scale,  normalized to  (0, 1)
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ANZ0001.sub
#' @usage data(ANZ0001.sub)
#' @references Stockler, M., T. Sourjina, P. Grimison, V. Gebski, M. Byrne, V. Harvey, P. Francis et al. 
#' ``A randomized trial of capecitabine (C) given intermittently (IC) rather than continuously (CC) compared to 
#' classical CMF as first-line chemotherapy for advanced breast cancer (ABC).'' In \emph{ASCO Annual Meeting 
#' Proceedings}, vol. 25, no. 18_suppl, p. 1031. 2007.
#' @format A data frame with 428 rows and 11 variables
NULL

#' @title Neck pain data set
#' @description A subset from an Australian chronic neck pain study 
#' @details  A randomized, double-blind, placebo-controlled study of low-level laser therapy (LLLT) in 88 
#' subjects with chronic neck pain was conducted with the aim of determining the efficacy of 300 mW, 
#' 830 nm laser in the management of chronic neck pain. Subjects were randomized to receive a course of 
#' 14 treatments over 7 weeks with either active or sham laser to tender areas in the neck. The primary 
#' outcome measure was change in a 10 cm Visual Analogue Scale (VAS) for pain.
#' Measurements were taken at baseline, at the end of 7 weeks\' treatment and 12 weeks from baseline.
#'
#' The variables are as follows:
#'
#' \tabular{ll}{
#' \code{id}\tab patient ID number\cr
#' \code{vas}\tab Neck pain as recorded by the patient on a VAS scale,  normalized to (0, 1)\cr
#' \code{laser}\tab laser treatment received by  patient, either 1 (active) or 2 (placebo)\cr
#' \code{time}\tab the measurement time, either 1 (initial assessment), 2 (after 7 weeks) or 3 (after 12 weeks). \cr
#' }
#'
#' @docType data
#' @keywords datasets
#' @name neck_pain
#' @usage data(neck_pain)
#' @references Chow RT, Heller GZ, Barnsley L (2006). ``The effect of 300 mW, 830 nm laser on chronic 
#' neck pain: a double-blind, randomized, placebo-controlled study.'' Pain, 124(1-2), 201-10. doi:16806710.
#' @format A data frame with 264 rows and 4 variables
NULL
