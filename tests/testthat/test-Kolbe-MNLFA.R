test_that("Testing MNLFA from ", {
  # The following is copied directly from Laura Kolbe, Terrence D. Jorgensen, Suzanne Jak, and Dylan Molenaar
  # at https://osf.io/527zr

  set.seed(123)

  ### Laura Kolbe
  ### Last updated: 22 February 2022
  ### MNLFA in OpenMx with mokken::DS14 data.
  ### Supplementary material


  ## -------------------------------
  ## Step 1: Install and load OpenMx
  ## -------------------------------

  library(OpenMx)


  ## -----------------------------
  ## Step 2: Load and prepare data
  ## -----------------------------

  library(mokken)

  ## Load data
  data("DS14", package="mokken")

  ## Save as data frame
  DS14 <- data.frame(DS14)

  ## Recode negatively worded items
  DS14$Si1 <- 4 - DS14$Si1.
  DS14$Si3 <- 4 - DS14$Si3.

  ## Standardize age
  DS14$Age <- (DS14$Age - mean(DS14$Age))/sd(DS14$Age)
  # mean-centering age is another option, but caused convergence difficulties with
  # this dataset

  ## Change order of variables
  DS14 <- DS14[,c("Male","Age","Si1","Si3","Si6","Si8","Si10","Si11","Si14",
                  "Na2","Na4","Na5","Na7","Na9","Na12", "Na13")]

  ## Check data
  head(DS14)

  ## Create data object
  mxdata1 <- mxData(observed=DS14, type="raw")

  ## Indicate names and number of indicators
  manVars <- colnames(DS14[,-c(1,2)])
  nv <- length(manVars)


  ## ------------------------------------
  ## Step 3: Test full scalar invariance
  ## ------------------------------------

  ## Specify matrices for configural model
  matT0 <- mxMatrix(type="Full", nrow=1, ncol=nv,
                    free=TRUE,
                    values=1,
                    name="matT0")
  matB1 <- mxMatrix(type="Full", nrow=1, ncol=nv,
                    free=TRUE,
                    values=0,
                    name="matB1")
  matB2 <- mxMatrix(type="Full", nrow=1, ncol=nv,
                    free=TRUE,
                    values=0,
                    name="matB2")
  matL0 <- mxMatrix(type="Full", nrow=nv, ncol=2,
                    free=c(rep(c(TRUE,FALSE),7),rep(c(FALSE,TRUE),7)),
                    values=c(rep(c(1,0),7),rep(c(0,1),7)),
                    byrow=TRUE,
                    name="matL0")
  matC1 <- mxMatrix(type="Full", nrow=nv, ncol=2,
                    free=c(rep(c(TRUE,FALSE),7),rep(c(FALSE,TRUE),7)),
                    byrow=TRUE,
                    values=0,
                    name="matC1")
  matC2 <- mxMatrix(type="Full", nrow=nv, ncol=2,
                    free=c(rep(c(TRUE,F),7),rep(c(FALSE,TRUE),7)),
                    byrow=TRUE,
                    values=0,
                    name="matC2")
  matE0 <- mxMatrix(type="Diag", nrow=nv, ncol=nv,
                    free=TRUE,
                    values=1,
                    name="matE0")
  matD1 <- mxMatrix(type="Diag", nrow=nv, ncol=nv,
                    free=TRUE,
                    values=0,
                    name="matD1")
  matD2 <- mxMatrix(type="Diag", nrow=nv, ncol=nv,
                    free=TRUE,
                    values=0,
                    name="matD2")
  matP0 <- mxMatrix(type="Symm", nrow=2, ncol=2,
                    free=c(FALSE,TRUE,TRUE,FALSE),
                    values=c(1,0,0,1),
                    name="matP0")
  matH1 <- mxMatrix(type="Symm", nrow=2, ncol=2,
                    free=c(FALSE,TRUE,TRUE,FALSE), # to identify the model
                    values=0,
                    name="matH1")
  matH2 <- mxMatrix(type="Symm", nrow=2, ncol=2,
                    free=c(FALSE,TRUE,TRUE,FALSE),  # to identify the model
                    values=0,
                    name="matH2")
  matA0 <- mxMatrix(type="Full", nrow=2, ncol=1,
                    free=FALSE,
                    values=0,
                    name="matA0")
  matG1 <- mxMatrix(type="Full", nrow=2, ncol=1,
                    free=FALSE, # to identify the model
                    values=0,
                    name="matG1")
  matG2 <- mxMatrix(type="Full", nrow=2, ncol=1,
                    free=FALSE, # to identify the model
                    values=0,
                    name="matG2")

  matV1 <- mxMatrix(type="Full", nrow=1, ncol=1,
                    free=FALSE,
                    labels="data.Male",
                    name = "Male")
  matV2 <- mxMatrix(type="Full", nrow=1, ncol=1,
                    free=FALSE,
                    labels="data.Age",
                    name = "Age")
  matIa <- mxMatrix(type="Diag", nrow=2, ncol=2,
                    free=FALSE,
                    values=1,
                    name="matIa")
  matIb <- mxMatrix(type="Full", nrow=2, ncol=2,
                    free=FALSE,
                    values=c(0,1,1,0),
                    name="matIb")

  ## Specify algebra for the dependent parameters
  matT <- mxAlgebra(expression=matT0+matB1*Male+matB2*Age,
                    name="matT")
  matL <- mxAlgebra(expression=matL0+matC1*Male+matC2*Age,
                    name="matL")
  matE <- mxAlgebra(expression=matE0*exp(matD1*Male+matD2*Age),
                    name="matE")
  matA <- mxAlgebra(expression=matA0+matG1*Male+matG2*Age,
                    name="matA")

  ## Specify algebra for covariance matrix of factors (transformed to ensure positive definite matrices)
  matVar <- mxAlgebra(expression=(matP0*exp(matH1*Male+matH2*Age)),
                      name="matVar")
  matR <- mxAlgebra(expression=(exp(2*(matP0+matH1*Male+matH2*Age))-1)/
                      (exp(2*(matP0+matH1*Male+matH2*Age))+1),
                    name="matR")
  matCov <- mxAlgebra(expression=(matIa*sqrt(matVar))%*%matR%*%(matIa*sqrt(matVar)),
                      name="matCov")
  matP <- mxAlgebra(expression=matIa*matVar+matIb*matCov,
                    name="matP")

  ## Specify model-implied matrices
  matC <- mxAlgebra(expression=matL%*%matP%*%t(matL)+matE,
                    name="matC")
  matM <- mxAlgebra(expression=matT+t(matL%*%matA),
                    name="matM")

  ## Specify expectation and fit function
  expF <- mxExpectationNormal(covariance="matC",
                              means="matM",
                              dimnames=manVars)
  fitF <- mxFitFunctionML()

  ## Make mxModel object and run the model
  modConfig <- mxModel(model="Configural",
                       matT, matT0, matB1, matB2,
                       matL, matL0, matC1, matC2,
                       matE, matE0, matD1, matD2,
                       matP, matP0, matH1, matH2,
                       matA, matA0, matG1, matG2,
                       matIa, matIb, matV1, matV2,
                       matVar, matR, matCov, matM, matC,
                       expF, fitF, mxdata1)
  fitConfig <- mxRun(modConfig)
  summary(fitConfig)

  #### Same model with mxsem ####
  library(mxsem)
  mnlfa_syntax <- "
==== MNLFA ====

SI =~ {lSi_1  := lSi0_1  + lSi1_1*data.Age  + lSi2_1*data.Male }*Si1  +
      {lSi_3  := lSi0_3  + lSi1_3*data.Age  + lSi2_3*data.Male }*Si3  +
      {lSi_6  := lSi0_6  + lSi1_6*data.Age  + lSi2_6*data.Male }*Si6  +
      {lSi_8  := lSi0_8  + lSi1_8*data.Age  + lSi2_8*data.Male }*Si8  +
      {lSi_10 := lSi0_10 + lSi1_10*data.Age + lSi2_10*data.Male}*Si10 +
      {lSi_11 := lSi0_11 + lSi1_11*data.Age + lSi2_11*data.Male}*Si11 +
      {lSi_14 := lSi0_14 + lSi1_14*data.Age + lSi2_14*data.Male}*Si14

NA =~ {lNa_2  := lNa0_2  + lNa1_2*data.Age  + lNa2_2*data.Male }*Na2  +
      {lNa_4  := lNa0_4  + lNa1_4*data.Age  + lNa2_4*data.Male }*Na4  +
      {lNa_5  := lNa0_5  + lNa1_5*data.Age  + lNa2_5*data.Male }*Na5  +
      {lNa_7  := lNa0_7  + lNa1_7*data.Age  + lNa2_7*data.Male }*Na7  +
      {lNa_9  := lNa0_9  + lNa1_9*data.Age  + lNa2_9*data.Male }*Na9  +
      {lNa_12 := lNa0_12 + lNa1_12*data.Age + lNa2_12*data.Male}*Na12 +
      {lNa_13 := lNa0_13 + lNa1_13*data.Age + lNa2_13*data.Male}*Na13

SI ~~ 1*SI
NA ~~ 1*NA + {cov := cov0  + cov1*data.Age  + cov2*data.Male }*SI

Si1  ~~ {vSi_1  := exp(vSi0_1  + vSi1_1*data.Age  + vSi2_1*data.Male )}*Si1
Si3  ~~ {vSi_3  := exp(vSi0_3  + vSi1_3*data.Age  + vSi2_3*data.Male )}*Si3
Si6  ~~ {vSi_6  := exp(vSi0_6  + vSi1_6*data.Age  + vSi2_6*data.Male )}*Si6
Si8  ~~ {vSi_8  := exp(vSi0_8  + vSi1_8*data.Age  + vSi2_8*data.Male )}*Si8
Si10 ~~ {vSi_10 := exp(vSi0_10 + vSi1_10*data.Age + vSi2_10*data.Male)}*Si10
Si11 ~~ {vSi_11 := exp(vSi0_11 + vSi1_11*data.Age + vSi2_11*data.Male)}*Si11
Si14 ~~ {vSi_14 := exp(vSi0_14 + vSi1_14*data.Age + vSi2_14*data.Male)}*Si14

Na2  ~~ {vNa_2  := exp(vNa0_2  + vNa1_2*data.Age  + vNa2_2*data.Male )}*Na2
Na4  ~~ {vNa_4  := exp(vNa0_4  + vNa1_4*data.Age  + vNa2_4*data.Male )}*Na4
Na5  ~~ {vNa_5  := exp(vNa0_5  + vNa1_5*data.Age  + vNa2_5*data.Male )}*Na5
Na7  ~~ {vNa_7  := exp(vNa0_7  + vNa1_7*data.Age  + vNa2_7*data.Male )}*Na7
Na9  ~~ {vNa_9  := exp(vNa0_9  + vNa1_9*data.Age  + vNa2_9*data.Male )}*Na9
Na12 ~~ {vNa_12 := exp(vNa0_12 + vNa1_12*data.Age + vNa2_12*data.Male)}*Na12
Na13 ~~ {vNa_13 := exp(vNa0_13 + vNa1_13*data.Age + vNa2_13*data.Male)}*Na13

Si1  ~ {iSi_1  := iSi0_1  + iSi1_1*data.Age  + iSi2_1*data.Male }*1
Si3  ~ {iSi_3  := iSi0_3  + iSi1_3*data.Age  + iSi2_3*data.Male }*1
Si6  ~ {iSi_6  := iSi0_6  + iSi1_6*data.Age  + iSi2_6*data.Male }*1
Si8  ~ {iSi_8  := iSi0_8  + iSi1_8*data.Age  + iSi2_8*data.Male }*1
Si10 ~ {iSi_10 := iSi0_10 + iSi1_10*data.Age + iSi2_10*data.Male}*1
Si11 ~ {iSi_11 := iSi0_11 + iSi1_11*data.Age + iSi2_11*data.Male}*1
Si14 ~ {iSi_14 := iSi0_14 + iSi1_14*data.Age + iSi2_14*data.Male}*1

Na2  ~ {iNa_2  := iNa0_2  + iNa1_2*data.Age  + iNa2_2*data.Male }*1
Na4  ~ {iNa_4  := iNa0_4  + iNa1_4*data.Age  + iNa2_4*data.Male }*1
Na5  ~ {iNa_5  := iNa0_5  + iNa1_5*data.Age  + iNa2_5*data.Male }*1
Na7  ~ {iNa_7  := iNa0_7  + iNa1_7*data.Age  + iNa2_7*data.Male }*1
Na9  ~ {iNa_9  := iNa0_9  + iNa1_9*data.Age  + iNa2_9*data.Male }*1
Na12 ~ {iNa_12 := iNa0_12 + iNa1_12*data.Age + iNa2_12*data.Male}*1
Na13 ~ {iNa_13 := iNa0_13 + iNa1_13*data.Age + iNa2_13*data.Male}*1
"

  library(mxsem)
  mnlfa_model <- mxsem(model = mnlfa_syntax,
                       data = DS14,
                       # we scaled the latent variables manually,
                       # so we will set all automatic scalings to FALSE:
                       scale_loadings = FALSE,
                       scale_latent_variances = FALSE)

  mnlfa_model <- mxRun(mnlfa_model)

  #### Comparison ####

  testthat::expect_equal(logLik(fitConfig), logLik(mnlfa_model), tolerance = 0.001)

  # Intercepts of regression coefficients:
  testthat::expect_lt(max(abs(fitConfig$matT0$values - c(coef(mnlfa_model)[grepl("^iSi0_", names(coef(mnlfa_model)))],
                                                         coef(mnlfa_model)[grepl("^iNa0_", names(coef(mnlfa_model)))]))), 0.001)

  # Regression coefficients for Age
  testthat::expect_lt(max(abs(fitConfig$matB2$values -
                                c(coef(mnlfa_model)[grepl("^iSi1_", names(coef(mnlfa_model)))],
                                  coef(mnlfa_model)[grepl("^iNa1_", names(coef(mnlfa_model)))]))), 0.001)

  # Regression coefficients for Gender
  testthat::expect_lt(max(abs(fitConfig$matB1$values -
                                c(coef(mnlfa_model)[grepl("^iSi2_", names(coef(mnlfa_model)))],
                                  coef(mnlfa_model)[grepl("^iNa2_", names(coef(mnlfa_model)))]))), 0.001)

  # Expectations
  testthat::expect_lt(max(abs(mxGetExpected(fitConfig, "covariance") -
                                mxGetExpected(mnlfa_model, "covariance"))), 0.001)

  testthat::expect_lt(max(abs(mxGetExpected(fitConfig, "means") -
                                mxGetExpected(mnlfa_model, "means"))), 0.001)

  # Model matrices
  testthat::expect_lt(max(abs(mxEval(A, mnlfa_model, compute=TRUE)[1:14, c("SI", "NA")] -
                                mxEval(matL, fitConfig, compute=TRUE))), 0.001)

  testthat::expect_lt(max(abs(mxEval(S, mnlfa_model, compute=TRUE)[1:14, 1:14] -
                                mxEval(matE, fitConfig, compute=TRUE))), 0.001)

  testthat::expect_lt(max(abs(mxEval(M, mnlfa_model, compute=TRUE)[,1:14] -
                                mxEval(matM, fitConfig, compute=TRUE))), 0.001)

})
