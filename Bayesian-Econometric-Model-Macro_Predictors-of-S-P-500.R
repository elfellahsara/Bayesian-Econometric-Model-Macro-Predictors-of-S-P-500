#======================================
# Bayesian Econometric Model: Macro Predictors of S&P 500 
# Author: Sara El Fellah 
# Description: OLS vs Bayesian regression (rstanarm + brms) 
#======================================

#1. Install Packages
packages <- c("tidyverse", "fredr","rstanarm","brms",
              "bayesplot","loo","lubridate","forecast",
              "Metrics", "patchwork", "scales")
installed <- rownames(installed.packages())
to_install <- packages[!packages %in% installed]
if (length(to_install)) install.packages(to_install)
library(tidyverse)
library(fredr)
library(rstanarm)
library(brms)
library(bayesplot)
library(loo)
library(lubridate)
library(forecast)
library(Metrics)
library(patchwork)
library(scales)
library(broom)

#2. FRED API KEY 
# Free key accesss at https://fred.stlouisfed.org/docs/api/api_key.html
fredr_set_key("954c0e284ebcadc7e8a3bb40c10b9db7")

#3. Pull Data from FRED 
start_date <- as.Date("2000-01-01")
end_date <- as.Date("2024-12-01")

pull_fred <- function(series_id,col_name){
  fredr(series_id = series_id,
        observation_start = start_date,
        observation_end = end_date,
        frequency = "m") %>%
    select(date,value) %>%
    rename(!!col_name:=value)
}

sp500 <- pull_fred("SP500", "sp500") #S&P 500 Index
cpi <- pull_fred("CPIAUCSL", "cpi")  #CPI (inflation proxy)
unrate <- pull_fred("UNRATE", "unemployment") #Unemployment rate 
fedfunds <- pull_fred("FEDFUNDS", "interest_rate") #Fed Funds Rate
indpro <- pull_fred("INDPRO", "indpro") #Industrial Production

#4. Compute Monthly Returns & YoY Changes 
sp500_ret <- sp500%>%
arrange(date) %>%  
mutate(sp500_return = (sp500/lag(sp500)-1)*100) %>%
  drop_na()
cpi_chg <- cpi %>%
  arrange(date) %>%
  mutate(inflation = (cpi/lag(cpi,12)-1)*100) %>%
  drop_na()
indpro_chg <- indpro %>% 
  arrange(date) %>% 
  mutate(indpro_growth = (indpro/lag(indpro,12)-1)*100) %>%
  drop_na()

#5. Merge Into Master Dataset
df <- sp500_ret %>%
  select(date,sp500_return) %>% 
  inner_join(cpi_chg    %>% select(date,inflation), by="date") %>%
  inner_join(unrate    %>% select(date, unemployment),by="date") %>%
  inner_join(fedfunds  %>% select(date,interest_rate),by="date") %>%
  inner_join(indpro_chg %>% select(date,indpro_growth),by="date") %>%
  drop_na() %>%
  arrange(date)

cat("Dataset: ", nrow(df), "monthly observations from", 
    format(min(df$date), "%b %Y"), "to", format(max(df$date),"%b %Y"),"\n")

# Exploratory Data Analysis
summary(df)
#Correlation heatmap 
cor_matrix <- df %>%
  select(-date) %>%
  cor () %>%
  round(2)

cor_long <- as.data.frame(as.table(cor_matrix)) %>%
  rename(Var1 = Var1, Var2 = Var2, Correlation=Freq)

p_cor <- ggplot(cor_long, aes(Var1, Var2, fill=Correlation))+
  geom_tile(color="white")+
  geom_text(aes(label=Correlation),size=3)+
  scale_fill_gradient2(low="steelblue",mid="white",
                       high="firebrick",midpoint=0)+
  labs(title="Correlation Matrix: S&P 500 Returns & Macro Variables",
       x=NULL, y=NULL)+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=30,hjust=1))
print(p_cor)

#7. Train/Test Split (80/20)
n <- nrow(df)
train_n <- floor(0.8*n)
df_train <- df[1:train_n,]
df_test <- df[(train_n +1):n, ]
cat("Training obs:",nrow(df_train),
    "| Test obs:",nrow(df_test),"\n")
#8. Model 1 - OLS Regression 
ols_model <- lm(
  sp500_return ~ inflation + unemployment + interest_rate + indpro_growth, 
  data=df_train
)
summary(ols_model)
#OLS Predictions on test set
ols_preds <- predict(ols_model, newdata=df_test, interval="prediction")
ols_rmse <- rmse(df_test$sp500_return,ols_preds[, "fit"])
ols_mae <- mae(df_test$sp500_return, ols_preds[,"fit"])
cat("OLS RMSE:", round(ols_rmse,4),
    "|MAE:", round(ols_mae,4),"\n")
#9. Model 2A - Bayesian Regression (rstanarm)
#Weakly informative priors (normal(0,2.5)by default in rstanarm)
bayes_stan <- stan_glm(
  sp500_return ~ inflation + unemployment + interest_rate + indpro_growth,
  data = df_train,
  family=gaussian(),
  prior=normal(location=0,scale=2.5,autoscale=TRUE),
  prior_intercept = normal(0,10),
  chains =4, 
  iter=4000,
  warmup=1000,
  seed=42
  
)
  
print(summary(bayes_stan),digits=4)

#Posterior summaries
posterior_interval(bayes_stan,prob=0.95)

#10. Model 2B - Bayesian Regression (brms)
bayes_brms <- brm(
  sp500_return ~ inflation + unemployment + interest_rate + indpro_growth,
  data = df_train,
  family=gaussian(),
  prior = c(
    prior(normal(0,2.5), class=b),
    prior(normal(0,10), class=Intercept),
    prior(exponential(1),class=sigma)
  ),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  seed = 42,
  cores = parallel::detectCores()-1
)
summary(bayes_brms)
#11. MCMC Diagnostics 
#Trace plots (rstanarm)
p_trace <- mcmc_trace(
  as.array(bayes_stan),
  pars=c("inflation","unemployment","interest_rate","indpro_growth")
)
print(p_trace)
#R-hat convergence check 
rhats <- rhat(bayes_stan)
cat("R-hat values(should be <1.01):\n")
print(round(rhats,4))
#Effective sample size
cat("Effective Sample Sizes:\n")
print(round(neff_ratio(bayes_stan),4))
#12. Posterior Distributions 
posterior_samples <- as.data.frame(bayes_stan)
coef_names <- c("inflation", "unemployment", "interest_rate","indpro_growth")
posterior_long <- posterior_samples %>%
  select(all_of(coef_names)) %>% 
  pivot_longer(everything(), names_to = "variable", values_to ="value")
p_posterior <- ggplot(posterior_long,aes(x=value, fill=variable))+
  geom_density(alpha=0.7, color=NA)+
  facet_wrap(~variable,scales="free")+
  geom_vline(xintercept=0,linetype="dashed", color="black")+
  scale_fill_brewer(palette="Set2")+
  labs(title = "Posterior Distributions of Macro Coefficients (rstanarm)",
       x="Coefficient Vlue", y="Density") +
  theme_minimal()+ 
  theme(legend.position="none")
print(p_posterior)
#13. Posterior Predictive Check 
p_ppc <- pp_check(bayes_stan,nreps=100)+
  labs(title="Posterior Predictive Check rstanarm)",
       x="S&P 500 Monthly Return (%)")
print(p_ppc)
#14. Predictions & Credible Intervals 
# rstanarm posterior predictions on test set 
bayes_pred_dist <- posterior_predict(bayes_stan,newdata=df_test)
bayes_fit <- colMeans(bayes_pred_dist)
bayes_lower <- apply(bayes_pred_dist,2,quantile,probs=0.025)
bayes_upper <- apply(bayes_pred_dist,2,quantile,probs=0.975)
bayes_rmse <- rmse(df_test$sp500_return,bayes_fit)
bayes_mae <- mae(df_test$sp500_return,bayes_fit)
cat("Bayes RMSE:",round(bayes_rmse,4),
    "|MAE:",round(bayes_mae,4),"\n")
#15. Forecast Comparison Plot 
forecast_df <- df_test %>%
select(date, sp500_return) %>% 
  mutate(
    ols_pred = ols_preds[, "fit"],
    ols_lower = ols_preds[, "lwr"],
    ols_upper = ols_preds[, "upr"], 
    bayes_pred = bayes_fit, 
    bayes_lower = bayes_lower, 
    bayes_upper = bayes_upper
  )
p_forecast <- ggplot(forecast_df, aes(x=date))+
  #Bayesian credible interval
  geom_ribbon(aes(ymin=bayes_lower,ymax=bayes_upper),
              fill="steelblue",alpha=0.2)+
  #OLS prediction interval 
  geom_ribbon(aes(ymin=ols_lower,ymax=ols_upper), 
              fill="firebrick",alpha=0.15)+
  #Actual Returns 
  geom_line(aes(y=sp500_return,color="Actual"),linewidth = 0.8)+
  #OLS forecast
  geom_line(aes(y=ols_preds,color="Bayesian"),linewidth=0.8,linetype="dotdash")+
  scale_color_manual(values=c("Actual"="black",
                              "OLS" = "firebrick",
                              "Bayesian"="steelblue"))+
  labs(title = "S&P 500 Monthly Return Forecasts: OLS vs. Bayesian", 
       subtitle="Shaded bands=95% prediction/credible intervals",
       x=NULL, y="Monthly Return(%)", color=NULL)+
  theme_minimal()+ 
  theme(legend.position = "bottom")
print(p_forecast)
#16. Coefficient Comparison: OLS vs Bayes 
ols_coefs <- broom::tidy(ols_model,conf.int=TRUE) %>%
filter(term != "(Intercept)") %>% 
  mutate(model ="OLS")
bayes_coefs <- as.data.frame(bayes_stan) %>%
  select(all_of(coef_names)) %>%
  pivot_longer(everything(), names_to="term", values_to="value") %>%
  group_by(term) %>% 
  summarise(
    estimate = mean(value), 
    conf.low = quantile(value, 0.025), 
    conf.high = quantile(value, 0.975),
    .groups = "drop"
  ) %>%
  mutate(model="Bayesian")
coef_df <- bind_rows(
  ols_coefs %>% select(term, estimate, conf.low, conf.high, model),
  bayes_coefs %>% select(term, estimate, conf.low, conf.high, model)
)
p_coefs <- ggplot(coef_df, aes(x = term, y = estimate,
                               color = model, shape = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(position = position_dodge(0.4), size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                position = position_dodge(0.4), width = 0.25, linewidth = 0.8) +
  scale_color_manual(values = c("OLS" = "firebrick", "Bayesian" = "steelblue")) +
  labs(title = "Coefficient Estimates: OLS vs Bayesian (95% Intervals)",
       x = "Predictor", y = "Coefficient", color = NULL, shape = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

cat("p_coefs exists:", exists("p_coefs"), "\n")
print(p_coefs)

#17. Model Comparison Table 
results_table <- tibble(
  Model = c("OLS", "Bayesian(rstanarm"), 
  RMSE = round(c(ols_rmse, bayes_rmse), 4), 
  MAE = round(c(ols_mae, bayes_mae), 4)
)
print(results_table)

#LOO Cross-Validation for Bayesian Model 
loo_result <- loo(bayes_stan)
print(loo_result)
#18. Save Outputs
ggsave("correlation_heatmap.png", p_cor, width=8, height=6, dpi=150)
ggsave("posterior_distributions.png", p_posterior, width = 10, height = 6, dpi = 150)
ggsave("posterior_predictive_check.png", p_ppc,  width = 8, height = 5, dpi = 150)
ggsave("forecast_comparison.png",  p_forecast, width = 10, height = 5, dpi = 150)
ggsave("coefficient_comparison.png", p_coefs,  width = 8, height = 5, dpi = 150)

write_csv(results_table, "model_performance.csv")
write_csv(forecast_df,   "forecast_results.csv")

cat("\nAll outputs saved.\n")
