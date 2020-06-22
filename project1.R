library(lme4)
lamb =read.csv("lamb_lab1.csv")
#preprocessing data
lamb$sire = factor(lamb$sire)
lamb$line = factor(lamb$line)
lamb$damage = factor(lamb$damage)

#Maximum Likelihood method
lamb_mle = lmer(weight ~ line + damage - 1 + (1|sire), data = lamb, REML = F)
summary(lamb_mle)
fixef(lamb_mle)
ranef(lamb_mle)

#Asymptotic covariance matrix
sigma_e_mle = sigma(lamb_mle)
sigma_s_mle = sqrt(unlist(VarCorr(lamb_mle)))
R_mle = sigma_e_mle^2*diag(62)
Z_mle = getME(lamb_mle, "Z")
Z_mle = as.matrix(Z_mle)
V_mle = sigma_s_mle^2*Z_mle %*% t(Z_mle) + R_mle
 

#fisher information matrix of MLE method
V_mle_1 = Z_mle %*% t(Z_mle)
V_mle_2 = diag(62)
fisher_mle_11 = 1/2*sum(diag(solve(V_mle) %*% V_mle_1 %*% solve(V_mle) %*% V_mle_1))
fisher_mle_12 = 1/2*sum(diag(solve(V_mle) %*% V_mle_1 %*% solve(V_mle) %*% V_mle_2))
fisher_mle_21 = 1/2*sum(diag(solve(V_mle) %*% V_mle_2 %*% solve(V_mle) %*% V_mle_1))
fisher_mle_22 = 1/2*sum(diag(solve(V_mle) %*% V_mle_2 %*% solve(V_mle) %*% V_mle_2))
fisher_mle = matrix(c(fisher_mle_11,fisher_mle_12,fisher_mle_21,fisher_mle_22),2)
sigma_s_mle_asymptotic = sqrt(solve(fisher_mle)[1,1])
sigma_e_mle_asymptotic = sqrt(solve(fisher_mle)[2,2])

######method2: bootstrap for mle######
mySumm = function(mod) {
  c(sigma_e_mle_boot = sigma(mod)^2, sigma_s_mle_boot = unlist(VarCorr(mod)))
}
booted_mle = bootMer(lamb_mle, mySumm, nsim = 100, seed = 2047)
booted_mle
library(boot) # for nice print-out
booted_mle



######RMEL method######
lamb_reml = lmer(weight ~ line + damage - 1 + (1|sire), data = lamb, REML = T)
summary(lamb_reml)
fixef(lamb_reml)
ranef(lamb_reml)

#Asymptotic covariance matrix
sigma_e_reml = sigma(lamb_reml)
sigma_s_reml = sqrt(unlist(VarCorr(lamb_reml)))
R_reml = sigma_e_reml^2*diag(62)
Z_reml = getME(lamb_reml, "Z")
Z_reml = as.matrix(Z_reml)
X_reml = getME(lamb_reml, "X")
X_reml = as.matrix(X_reml)
V_reml = sigma_s_reml^2*Z_reml %*% t(Z_reml) + R_reml
P_reml = solve(V_reml) - solve(V_reml) %*% X_reml %*% solve(t(X_reml)%*%solve(V_reml)%*%X_reml) %*% t(X_reml) %*% solve(V_reml)
V_reml_1 = Z_reml %*% t(Z_reml)
V_reml_2 = diag(62)

fisher_reml_11 = 1/2*sum(diag(P_reml %*% V_reml_1 %*% P_reml %*% V_reml_1))
fisher_reml_12 = 1/2*sum(diag(P_reml %*% V_reml_1 %*% P_reml %*% V_reml_2))
fisher_reml_21 = 1/2*sum(diag(P_reml %*% V_reml_2 %*% P_reml %*% V_reml_1))
fisher_reml_22 = 1/2*sum(diag(P_reml %*% V_reml_2 %*% P_reml %*% V_reml_2))

fisher_reml = matrix(c(fisher_reml_11,fisher_reml_12,fisher_reml_21,fisher_reml_22),2)
sigma_s_reml_asymptotic = sqrt(solve(fisher_reml)[1,1])
sigma_e_reml_asymptotic = sqrt(solve(fisher_reml)[2,2])

#bootstrap(REML)
mySumm = function(mod) {
  c(sigma_e_reml_boot = sigma(mod), sigma_s_reml_boot = sqrt(unlist(VarCorr(mod))))
}
booted_reml = bootMer(lamb_reml, mySumm, nsim = 100, seed = 2047)
booted_reml
library(boot) # for nice print-out
booted_reml

#delete damage fiexed effect#
lamb_reml_noage = lamb_reml = lmer(weight ~ line - 1 + (1|sire), data = lamb, REML = F)
anova(lamb_reml_noage, lamb_mle)


