## Script to generate density measure based on history of BI-RADS density
## Authors: Emma Atakpa, Adam Brentnall
## Last update 17th November 2021

##1. Load parameters

load("model_output.RData")

##2. Function to calc long. density

EmmaDensity <- function(indata, sigma, B, betahat){

  ## indata
  ## columns age, bmi, density
  ## sigma - residual sd
  ## B - varcov matrix for random effects
  ## betahat - fixed effect estimates
  ##"EmmaDensity" function uses "sigma", "B", and "betahat" from "model_output.RData" and "age", "bmi", and "density" from the data that you input ("indata")
  
  ##order by age (ties will be left in their original order)
  rDat <- indata[order(indata$age),]

  ##winsorise bmi (so that, in risk modeling, morbidly obese women are taken to have the same risk due to adiposity as obese women)
  rDat$bmi <- pmin(rDat$bmi, 35)
  rDat$bmi <- pmax(rDat$bmi, 15)
  
  ##centre & scale age, centre bmi, create age polynomials, create intercept column
  rDat$myage <- rDat$age-40
  rDat$myage10 <- rDat$myage/10
  rDat$mybmi <- rDat$bmi-25
  rDat$sqmyage10 <- rDat$myage10^2
  rDat$cbmyage10 <- rDat$myage10^3
  rDat$fourthmyage10 <- rDat$myage10^4
  rDat$intercept <- 1
  
  ## covariates in the model - fixed effects
  x<-rDat[, c("intercept", "myage10", "mybmi", "sqmyage10", "cbmyage10", "fourthmyage10")]

  x$myage10_mybmi<-rDat$myage10*rDat$mybmi

  Xi<-as.matrix(x)
  
  ## random effects
  Zi<-as.matrix(x[,1:2])

  ## density
  yi<-rDat[, c("density")]
  
  ## number of points
  myni<-length(yi)
  
  ## Empirical bayes calculation
  Ei<-sigma^2 * diag(myni)

  SIGMAix<-Zi %*% B %*% (t(Zi))+ Ei

  SIGMAix.inv<-solve(SIGMAix)

  myre<-B %*% t(Zi) %*% SIGMAix.inv %*% (yi - Xi %*% betahat) ## random effects prediction (empirical bayes)

  mypred<-Xi %*% betahat + Zi %*% myre ##prediction
  
  ## longitudinal density also converted to an eight-category variable to help interpretation with BI-RADS density
  mypredcat<-as.integer(cut(mypred[myni,], c(-9999,1.54407238,1.95699767,2.19156056,2.60191337,2.89512705,3.24643613,3.59640892,9999)))

  ## output: longitudinal density prediction at last age (continuous),
  ##         longitudinal density prediction at last age (categorical),
  ##         random effects prediction at last age (3rd element of vector is random intercept, 4th element of vector is random slope).
  c(as.double(mypred[myni,]), mypredcat, myre)
}

## 3. Demo for one woman

rDat<-data.frame(age=c(50,52,55), bmi=c(25,26,24), density=c(3,3,2))

EmmaDensity(rDat, sigma, B, betahat)

## 4. Example data

load("exampledata.RData")

  ##mammo_label: number assigned to each mammogram
  ##studyid: study ID for each woman
  ##N: total number of mammograms per woman
  ##n: number assigned to each mammogram per woman (ordered chronologically)
  ##age: age at mammogram (years)
  ##density: BI-RADS density at mammogram
  ##bmi: BMI at mammogram (kg/m^2)

EmmaDensity(exampledata[exampledata$studyid=="ID1", c("age","bmi","density")], sigma, B, betahat)

  ## [1]  2.20830737  4.00000000 -0.30304190  0.01617053

  ## 1st element: Longitudinal density prediction at last mammogram for ID1 (continuous)
  ## 2nd element: Longitudinal density prediction at last mammogram for ID1 (categorical) - this is an integer from 1 to 8, but displays as a numerical
  ## 3rd element: Random intercept prediction at last mammogram for ID1
  ## 4th element: Random slope prediction at last mammogram for ID1

  ##substitute "ID1" with each studyid to get the output for each woman

