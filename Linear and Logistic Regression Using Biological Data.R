# Jesus Zeno

# Part 1 Linear Regression Using Parkinson's data

#Let's start by writing a function that will install a group of packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# Going to turn it into a different construction and call in libraries
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# Send the package names to a vector to be installed at once later
packages <- c("Hmisc","corrplot","PerformanceAnalytics","dominanceanalysis","pscl")

# call function to install the list of packages. Install if they haven't been yet.
ipak(packages)

# call in all the libraries
library(Hmisc)
library(corrplot)
library(PerformanceAnalytics)
library(pscl)
library(dominanceanalysis)

# read in the Parkinson's dataset with the headers as TRUE to dictate column 
# names.
parkinsons = read.csv(file.choose(),header=TRUE)
names(parkinsons)

# This will get the pairwise correlation of all the columns in the dataset.
correlationMatrix=cor(parkinsons,method="pearson",use="complete.obs")

#This will tell us if they are significantly correlated or not.
correlationMatrixSignificance=rcorr(as.matrix(parkinsons))

# This will give the pvalues for that.
pvals=correlationMatrixSignificance$P

# Another look at the correlation matrix. 
flattenedMatrix=flattenCorrMatrix(correlationMatrixSignificance$r,correlationMatrixSignificance$P)

# useful visual
symnum(correlationMatrix,abbr.colnames=FALSE)

# Call a correlation plot to analyze the above correlation numbers
corrplot(correlationMatrix,type="upper",order="original",tl.col="black",tl.srt=45)
dev.off()

# Conclusion: There is s high positive correlation between motor_UPDRS
# and total_UPDRS. There is also a high positive correlation within
# the Jitter group and the Shimmer group. There is a slightly lesser
# positive correlation between the Jitter groups and Shimmer groups.
# There is also a high positive correlation between NHR and the Jitter
# groups or Shimmer groups. In addition, there is a strong negative
# correlation HNR and the Jitter group, Shimmer group, NHR, RPDE, and 
# PPE. RPDE and PPE show much weaker positive correlations to other
# groups such as the Jitter group, Shimmer group, and NHR. RPDE and 
# PPE also showed positive correlation between each other, which makes
# sense given their similar correlation to other groups. 

# Let's create a performance plot to visualize the correlation matrix. 
# Shows the scatter plot for all variables against one another.
# Down diagonal there is scatter plot of individual variables.
# Due to it showing relation among all the variables, we have a lot of
# plots packed into a very tight space here. It's not the most visually
# appealing and might be better suited if it were done in smaller
# groups
chart.Correlation(parkinsons,histogram=TRUE,pch=19)
dev.off()

# Let's predict the total UPDRS from the motor. 
total_linear_predict = lm(total_UPDRS~motor_UPDRS, data = parkinsons)
# Show us what we need to know about the object
summary(total_linear_predict)
# Interpretation: the intercept is a non-zero value and has a p-value
# <.001 which means it is a significant intercept and the standard error
# is about 1/8. The intercept matters to this model. 
# Motor_UPDRS has a slope different than zero which implies there
# could be some pattern. The p-value is <.001 as well. 
# Adjusted R-squared value is 0.8972 which is close to 1 and a good 
# value.
abline(coef=coef(total_linear_predict))
plot(total_linear_predict)
# Looking at the residuals vs. fitted plot, there are a lot of values 
# that don't fit in a type of belt. They don't really have a cohesive
# pattern. Because there are less data points in a belt, that means
# I can expect a smaller spread of deviation. 
# Q-Q plot shows that the points are not a perfectly fitted
# to normal distribution dataset. We see that it starts to go above the 
# line at about a little larger than one standard deviation. 

# Let's predict the motor UPDRS from the total 
motor_linear_predict = lm(motor_UPDRS~total_UPDRS, data = parkinsons)
# Show us what we need to know about the object
summary(motor_linear_predict)
# Interpretation: looking at the intercept, it's a non-zero value but
# it's also less than one so I'm not convinced the intercept is 
# important to the prediction despite it having a p-value of <.001
# and a standard error of .098284. 
# Total_UPDRS has a slope different than zero which implies there
# could be some pattern. The p-value is <.001 as well. 
# Adjusted R-squared value is 0.8972 which is close to 1 and a good 
# value.

abline(coef=coef(motor_linear_predict))
plot(motor_linear_predict)
# Looking at the residuals vs. fitted plot, there are a lot of values 
# that don't fit in a type of belt. They don't really have a cohesive
# pattern. Because there are less data points in a belt, that means
# I can expect a smaller spread of deviation. 
# Q-Q plot shows that the points that there is not a perfectly fitted
# to normal distribution dataset. We see that it starts to go above the 
# line at about a little larger than two standard deviations and below
# the line as it goes past the -1 standard deviation. 

# Let's create a linear regression model to predict the motor and 
# total UPDRS scores from all the other variables. Let's do motor first.
motor_prediction = glm(as.factor(motor_UPDRS) ~ subject.+age+sex+test_time
                    +total_UPDRS+Jitter...+Jitter.Abs.+Jitter.RAP+
                      Jitter.PPQ5+Jitter.DDP+Shimmer+Shimmer.dB.+
                      Shimmer.APQ3+Shimmer.APQ5+Shimmer.APQ11+
                      Shimmer.DDA+NHR+HNR+RPDE+DFA+PPE, 
                    data=parkinsons, family=binomial(link='logit'))
summary(motor_prediction)
# Analysis: the intercept is a positive number (5.183e+02) which means
# that there is a positive correlation between motor_UPDRS and the 
# other variables. Unfortunately, the none of the p-values are 
# significant which leads me to believe this isn't a good model. 

# Now let's do total_UPDRS
total_prediction = glm(as.factor(total_UPDRS) ~ subject.+age+sex+
                         test_time+motor_UPDRS+Jitter...+Jitter.Abs.+
                         Jitter.RAP+Jitter.PPQ5+Jitter.DDP+Shimmer+
                         Shimmer.dB.+Shimmer.APQ3+Shimmer.APQ5+
                         Shimmer.APQ11+Shimmer.DDA+NHR+HNR+RPDE+DFA+PPE, 
                       data=parkinsons, family=binomial(link='logit'))
summary(total_prediction)
# Analysis: This is similar to the previous motor prediction model.
# The positive intercept is showing a positive correlation. However,
# the model doesn't seem good since none of the p-values are significant.

# Conclusion: Neither model is great to me, but I would have to choose
# the linear regression models using all of the other categories to
# predict the total or motor. We will always have more explanatory power
# when using more variables for prediction. I would probably use 
# dominance analysis to see which categories would be best to predict
# motor and total. This way we can hopefully see significant values
# when using a linear regression model with multiple categories for
# prediction.

# Part 2 Logistic Expression with breast cancer data
# Read in the breast cancer csv file
breast_cancer = read.csv(file.choose(),header=TRUE)

# Let's look at the column names in the dataset
names(breast_cancer)

# Let's setup a logistical regression model to try and predict if
# a cell is malignant (4) or benign (2). We will be using the cell
# columns from the data. We will be excluding the ID column from 
# the analysis since it is just the identifier of the sample. 
cancer_prediction = glm(as.factor(BenignMalignant)~ClumpThickness+
                          CellSizeUniformity+CellShapeUniformity+
                          MarginalAdhesion+SingleEpithelialCell+
                          BareNuclei+BlandChromatin+NormalNucleoli+
                          Mitoses, data=breast_cancer, 
                        family=binomial(link='logit'))
summary(cancer_prediction)
# Conclusion: The intercept, ClumpThickness, MarginalAdhesion, BareNuclei,
# and BlandChromatin are all significant terms at predicting if a 
# mass is benign or malignant since they have p-values <.05. 
# However, since this would be used in a medical setting
# for possible patient diagnosis, I would personally only consider
# the following as significant terms since their p-value is <.001:
# intercept, ClumpThickness, and BareNuclei. I would want a high degree
# of confidence when assessing a tumor as malignant or benign for 
# the reasons stated in class about how false negatives and false
# positives would be highly undesirable in this case.


# Let's setup a dominance analysis regarding the columns contributing
# to the logistic regression. 
cancer_dom_analysis<-dominanceAnalysis(cancer_prediction)


# General rubric for dominance matrices are as follows: 
# variables with values of 1 added to the model
# mean they add significant value to predicting the outcome. Variables
# with value of 0 mean they are submissive and don't contribute. 
# Values of .5 can go either way. It's a toss up. 

# Setup a complete dominance matrix
dominanceMatrix(cancer_dom_analysis, type="complete",fit.functions = "r2.m" )
# Clump thickness with 4 1s, BareNuclei with 5 1s, BlandChromatin with 2 1s, 
# and CellShapeUniformity with one 1. The remaining variables are .5 or 0. 
# This means the 4 aforementioned variables add more to the model than other
# variables. 

# Let's setup a conditional dominance matrix. Generally, we expect to see more 
# 1s in this model than the complete dominance one. 
dominanceMatrix(cancer_dom_analysis, type="conditional",fit.functions = "r2.m")
# Clump thickness with 4 1s, BareNuclei with 6 1s, BlandChromatin with 4 1s,
# and CellShapeUniformity with one 1. BareNuclei was the only one that increased
# the amount of 1s. This means that given other terms had been added to the 
# model this is a beneficial addition. 

# Let's setup a general dominance matrix. The contribution overall of something
# to the correlation is on average better. 
dominanceMatrix(cancer_dom_analysis, type="general",fit.functions = "r2.m")
# BareNuclei seems to have a general dominance over everything since the scores
# are always .5 or greater. Conversely, we see that mitosis always has values
# of .5 or less, which means that mitosis adds no value to the outcome. 
