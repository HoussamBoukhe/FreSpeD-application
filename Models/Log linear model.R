setwd("C:/Users/labo-maths/Desktop/FreSpeD-application")

data <- read.csv("Output/windowLen250-Welch2.csv")

df <- data[c("subject","region", "epoch","eyeState", "nCp_alpha")]
df$region <- relevel(as.factor(df$region), ref = "1") # region 1 as reference 
df$eyeState <- relevel(as.factor(df$eyeState), ref = "0") # eyes closed as reference

# Fit the Poisson regression model
model <- glm(nCp_alpha ~ region + epoch + eyeState, 
             data = df, family = poisson(link = "log"))

# Model summary
summary(model)

########
# The estimate of eyeState is negative meaning the number of change points 
# decreases when eyes are closed 
# The number of change points of region 7 increases with exp(0.91) comared 
# to region 1 

# estimated lambdas of the Poisson model 
log_lambda <- predict(model, type = "link") 
lambda <- exp(log_lambda)
lambda_o <- lambda[df$eyeState == 1]  
lambda_c <- lambda[df$eyeState == 0] 

plot(density(lambda_o), col = "blue", xlim = c(0,10),lwd = 2, main = "Density of lambda_o (Eyes Open) vs lambda_c (Eyes Closed)",
     xlab = "Expected lambda", ylab = "Density")
lines(density(lambda_c), col = "red", lwd = 2)
legend("topright", legend = c("lambda_o (Eyes Open)", "lambda_c (Eyes Closed)"),
       col = c("blue", "red"), lwd = 2)



############# Mixed effect model
library(lme4)

df <- data[c("subject","region", "epoch","eyeState", "nCp_alpha")]
df$region <- as.factor(df$region) 
df$eyeState <- as.factor(df$eyeState)


# Fit a Poisson mixed-effects model with random intercept for subjects
mixed_model <- glmer(nCp_alpha ~ region + epoch + eyeState + (1 | subject), 
               data = df, 
               family = poisson(link = "log"))

# View the model summary
summary(mixed_model)

# The standard deviation of the mixed effect model is 0.67 which suggests 
# variability among subjects 
