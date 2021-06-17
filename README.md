# MALD
Use multivariate asymmetric laplace distribution to improve jump modeling in cryptocurrency.

* Empirical_data_study.R: 2014-2020 BTC, S&P SVMALD, SVMVN, SVLD, SVIND.
 - Output: keeps_SVMALD, keeps_SVMVN, keeps_SVLD, keeps_SVIND rds files saved to keeps_long folder
* Empirical_data_study_meme.R: 2020 (BTC, GME, AMC, DOGE), S&P SVMALD, SVMVN, SVLD, SVIND
 - Output: 16 rds files saved to keeps_short folder
* Empirical_results.R: plots and tables showing results of MCMC resulting from above 2 programs
 - Outpue: 