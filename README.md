# neutral_imagefit_routine
Neutral experiment base cloud analysis routine

AnalysisVariables is the main configuration file for the routine and imagefit_Analysis_Runner is a wrapper for the three main functions described below.
1. imagefit_Backgrounds_PCA - background subtractino algorithm which utilizes principal component analysis on sets of images to better approximate the atomic background.
2. imagefit_NumDistFit - non-linear regression routine to fit the atom cloud image to the desired functional form
3. imagefit_ParamEval - plotting routine which reads in the results of the background subtraction and cluod fitting then makes all the relevant plots.
