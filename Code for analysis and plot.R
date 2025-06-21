#Data preparation
library(readxl)
data <- read_excel('data.xlsx', sheet = "Sheet1")

# Kolmogorov-Smirnov normality test of the outcome
column_data <- data$Congenital_birth_defects
ks_test_result <- ks.test(column_data, "pnorm", mean = mean(column_data), sd = sd(column_data))
print(ks_test_result)

#Test whether there is overdispersion in the Poisson distribution
library(qcc)
outcome <- c('Congenital_birth_defects','Orofacial_clefts','Congenital_musculoskeletal_and_limb_anomalies',
             'Neural_tube_defects','Digestive_congenital_anomalies','Klinefelter_syndrome',
            'Congenital_heart_anomalies','Down_syndrome','Urogenital_congenital_anomalies','Turner_syndrome')
for (var in outcome) {
  test_result <- qcc.overdispersion.test(data[[var]], type='poisson')
  print(test_result)
} #All p-values are 0, hence the data fit the negative binomial distribution.


#Spearman correlation
library(corrplot)
subset_data <- data[c('SDI','DryFrequency','DryDuration','DrySeverity',
                      'DryImpact_year_','Temp','Congenital_birth_defects')]
correlation_matrix <- cor(subset_data, method = "spearman")
p_values <- cor.mtest(subset_data, method = "spearman")$p
color_scheme <- colorRampPalette(c("#4477AA","#77AADD","#FFFFFF","#EE9988","#BB4444"))
corrplot(correlation_matrix, method = "color", col = color_scheme(200), 
         addCoef.col = "black", 
         number.cex = 1, 
         tl.col = "black",
         tl.srt = 30, 
         p.mat = p_values, 
         sig.level = 0.05, 
         insig = "blank" 
)

#Boxplot of the annual variation of drought data
library(ggplot2)
library(dplyr)
p1 <- ggplot(data, aes(x = factor(Year), y = DryImpact_year_)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  labs(title = "Boxplot of drought regions fraction by Year",
       x = "Year",
       y = "Drought regions fraction") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "gray80"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_blank())
print(p1)


#Negative binomial regression (main model)
library(MASS)
library(broom)
library(writexl)
response_var <- "Congenital_birth_defects"
covariate <- c('Temp',"SDI_rank","Year")
predictors <- c("DryFrequencyLag1", "DryDurationLag1", "DrySeverityLag1","DryImpact_year_Lag1")
run_NBRA_model <- function(predictor) {
  formula <- as.formula(paste(response_var, "~", predictor, "+", paste(covariate, collapse = "+")))
  model <- glm.nb(formula, data = data)
  summary_model <- summary(model)
  model_summary_df <- tidy(model)
  return(model_summary_df)
}
results <- lapply(predictors, run_NBRA_model)
results_list <- lapply(1:length(results), function(i) {
  data.frame(results[[i]], Model = paste("Model", i))
})
results_df <- do.call(rbind, results_list)
write_xlsx(list("Model_Results" = results_df), "Model_Results.xlsx")


#RCS
library(rms)
dd <- datadist(data)
options(datadist='dd')
fit <- ols(Congenital_birth_defectsP ~ rcs(DryFrequency,5)
           + SDI_rank +Year
           +rcs(Temp,3)
           ,data=data)
summary(fit)
out <- anova(fit)
OLS<-Predict(fit,DryFrequency)
plot(Predict(fit,DryFrequency),anova=out,pval=T,
     xlab = "Drought frequency", ylab = "Congenital birth defects"
     )

AIC(fit)



#Subgroup analysis based on SDI
library(dplyr)
data_SDI <- data %>% filter(SDI_rank == "High")
#Establish a loop
response_var <- "Congenital_birth_defects"
covariate <- c('Temp',"Year")
predictors <- c("DryFrequency", "DryDuration", "DrySeverity", "DryImpact_year_")
run_NBRA_model <- function(predictor) {
  formula <- as.formula(paste(response_var, "~", predictor, "+", paste(covariate, collapse = "+")))
  model <- glm.nb(formula, data = data_SDI)
  summary_model <- summary(model)
  model_summary_df <- tidy(model)
  return(model_summary_df)
}
results <- lapply(predictors, run_NBRA_model)
results_list <- lapply(1:length(results), function(i) {
  data.frame(results[[i]], Model = paste("Model", i))
})
results_df <- do.call(rbind, results_list)
write_xlsx(list("newModel_Results_High" = results_df), "newModel_Results_High.xlsx")

#Subgroup analysis of CBDs subtypes
response_var <- "Urogenital_congenital_anomalies"
covariate <- c('Temp',"SDI_rank","Year")
predictors <- c("DryFrequency", "DryDuration", "DrySeverity", "DryImpact_year_")
run_NBRA_model <- function(predictor) {
  formula <- as.formula(paste(response_var, "~", predictor, "+", paste(covariate, collapse = "+")))
  model <- glm.nb(formula, data = data)
  summary_model <- summary(model)
  model_summary_df <- tidy(model)
  return(model_summary_df)
}
results <- lapply(predictors, run_NBRA_model)
results_list <- lapply(1:length(results), function(i) {
  data.frame(results[[i]], Model = paste("Model", i))
})
results_df <- do.call(rbind, results_list)
write_xlsx(list("08.Urogenital_congenital_anomalies" = results_df), 
           "08.Urogenital_congenital_anomalies.xlsx")



#Global Burden of Disease Analysis
#Distribution of ASIR across all countries (world map)
library(readr)
library(dplyr)
library(sf)  
library(patchwork) 
library(ggplot2)
library(tidyverse)
data_2021 <- data %>% 
  filter(Year == "2021")
GBD <- data_2021

#World map
worldData <- map_data('world')
unique_values <- unique(worldData$region)  
unique_values_df <- data.frame(unique_values)  
main_map_data <- full_join(worldData,GBD,by = c('region'='location')) %>%   
  filter(Congenital_birth_defects != "NA")
dim(GBD)
head(GBD)

#Main figure plotting
fig <- main_map_data %>% 
  ggplot()+
  geom_polygon(aes(x = long, y = lat,group = group,fill=Congenital_birth_defects),
               colour="black",size=0.5) +
  theme_bw()+
  scale_fill_distiller(palette="Spectral",
                       name="ASIR(/10^5)") + 
  labs(x="",y="",title="")+ 
  theme(legend.position = c(0.1,0.2),
        legend.title = element_text(color="black", 
                                    size = 10, 
                                    #family = "A",
                                    #face = "bold"
        ),
        plot.title = element_text(color="black", 
                                  size = 14, 
                                  #family = "A",
                                  #face = "bold"
        ),
        legend.text = element_text(color="black", 
                                   size = 12, 
                                   #family = "A",
                                   #face = "bold"
        ),
        panel.grid=element_blank(),
        panel.border = element_blank(), 
        #legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
  )
fig

#Zoom in on specific areas
fig_central_african_republic <- main_map_data %>%
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group, fill = Congenital_birth_defects),
               colour = "black", size = 0.5) +
  theme_bw() +
  scale_fill_distiller(palette = "Spectral", name = "ASIR(/10^5)") +
  labs(x = "", y = "", title = "Central African Republic") +
  theme(
    legend.position = "none",  
    plot.title = element_text(color = "black", size = 14),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA), 
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  coord_sf(xlim = c(15, 30), ylim = c(2, 12))  

#Composite map
library(patchwork)
fig_combined <- fig + 
  (fig_central_african_republic | fig_brunei | fig_tajikistan | fig_cook_islands | fig_bermuda | fig_malaysia) +  # 使用 | 来并排放置小图
  plot_layout(heights = c(5, 2))  
print(fig_combined)


# EAPC calculation for all countries
library(tidyverse)
library(dplyr)
EAPC <- data  %>% 
  .[,c(3,17,60)]
head(EAPC)
EAPC_cal <- data.frame(location=unique(EAPC$location),
                       EAPC=rep(0,times=length(unique(EAPC$location))),
                       LCI=rep(0,times=length(unique(EAPC$location))),
                       UCI=rep(0,times=length(unique(EAPC$location))))
for (i in 1:length(unique(EAPC$location))){
  country_cal <- as.character(EAPC_cal[i,1])
  a <- subset(EAPC, EAPC$location==country_cal)
  a$y <- log(a$Congenital_birth_defects)
  mod_simp_reg<-lm(y~Year,data=a)
  estimate <- (exp(summary(mod_simp_reg)[["coefficients"]][2,1])-1)*100#per100,000
  low <- (exp(summary(mod_simp_reg)[["coefficients"]][2,1]-1.96*summary(mod_simp_reg)[["coefficients"]][2,2])-1)*100
  high <- (exp(summary(mod_simp_reg)[["coefficients"]][2,1]+1.96*summary(mod_simp_reg)[["coefficients"]][2,2])-1)*100
  EAPC_cal[i,2] <- estimate
  EAPC_cal[i,3] <- low
  EAPC_cal[i,4] <- high
}
EAPC_cal <- EAPC_cal %>% mutate(EAPC=round(EAPC,3),
                                LCI=round(LCI,3),
                                UCI=round(UCI,3))
EAPC_cal <- EAPC_cal %>% mutate(EAPC_CI = paste(EAPC, LCI,sep = '\n(')) %>% 
  mutate(EAPC_CI = paste(EAPC_CI, UCI,sep = ' to ')) %>% 
  mutate(EAPC_CI = paste0(EAPC_CI, ')'))
write_xlsx(list("EAPC_cal" = EAPC_cal), "EAPC_cal.xlsx")

#plot map of EAPC
main_map_data <- full_join(worldData,EAPC_cal,by = c('region'='location')) %>%   
  filter(EAPC != "NA")
dim(EAPC_cal)
head(EAPC_cal)

fig_EAPC <- main_map_data %>% 
  ggplot()+
  geom_polygon(aes(x = long, y = lat,group = group,fill=EAPC),
               colour="black",size=0.5) +
  theme_bw()+
  scale_fill_distiller(palette="Spectral",
                       name="EAPC(%)") + 
  labs(x="",y="",title="")+ 
  theme(legend.position = c(0.1,0.2),
        legend.title = element_text(color="black", 
                                    size = 10, 
                                    #family = "A",
                                    #face = "bold"
        ),
        plot.title = element_text(color="black", 
                                  size = 14, 
                                  #family = "A",
                                  #face = "bold"
        ),
        legend.text = element_text(color="black", 
                                   size = 12, 
                                   #family = "A",
                                   #face = "bold"
        ),
        panel.grid=element_blank(),
        panel.border = element_blank(),
        #legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
  )
fig_EAPC


#Zoom in on specific areas
fig_spain <- main_map_data %>%
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group, fill = EAPC),
               colour = "black", size = 0.5) +
  theme_bw() +
  scale_fill_distiller(palette = "Spectral", name = "EAPC(%)") +
  labs(x = "", y = "", title = "Spain") +
  theme(
    legend.position = "none",  
    plot.title = element_text(color = "black", size = 14),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),  
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  coord_sf(xlim = c(-10, 5), ylim = c(35, 44))

fig_combined <- fig_EAPC + 
  (fig_spain | fig_georgia | fig_american_samoa | fig_serbia | fig_saudi_arabia | fig_qatar) +  # 使用 | 来并排放置小图
  plot_layout(heights = c(5, 2))  
print(fig_combined)


#The distribution of the four types of droughts across all countries
GBD <- data_2021
worldData <- map_data('world')
main_map_data <- full_join(worldData,GBD,by = c('region'='location')) %>%   
  filter(DryFrequency != "NA")
dim(GBD)
head(GBD)
fig_Dryfrequency<- main_map_data %>% 
  ggplot()+
  geom_polygon(aes(x = long, y = lat,group = group,fill=DryFrequency),
               colour="black",size=0.5) +
  theme_bw()+
  scale_fill_distiller(palette="Spectral",
                       name="Drought frequency\n(times per year)") + 
  labs(x="",y="",title="")+ 
  theme(legend.position = c(0.1,0.2),
        legend.title = element_text(color="black", 
                                    size = 10, 
                                    #family = "A",
                                    #face = "bold"
        ),
        plot.title = element_text(color="black", 
                                  size = 14, 
                                  #family = "A",
                                  #face = "bold"
        ),
        legend.text = element_text(color="black", 
                                   size = 12, 
                                   #family = "A",
                                   #face = "bold"
        ),
        panel.grid=element_blank(),
        panel.border = element_blank(),  
        #legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
  )
fig_Dryfrequency


#Line chart
data_region_SDI <- data_region %>%
  filter(location_name %in% c("Global", "High SDI", "High-middle SDI", 
                             "Middle SDI", "Low-middle SDI", "Low SDI"))
library(ggplot2)
data_region_SDI$location_name <- factor(data_region_SDI$location_name,
                        levels = c("Global", "High SDI", "High-middle SDI",
                              "Middle SDI", "Low-middle SDI", "Low SDI"))
ggplot(data_region_SDI, aes(x = year, y = val, color = location_name,
                                              group = location_name)) +
  geom_line(size = 1) +                
  geom_point(size = 2) +          
  scale_color_manual(values = c("Global" = "#1f78b4", "High SDI" = "#33a02c", 
                                "High-middle SDI" = "#e31a1c", "Middle SDI" = "#ff7f00", 
                                "Low-middle SDI" = "#6a3d9a", "Low SDI" = "#b15928")) + 
  labs(title = "Congenital birth defects",
       x = "Year",
       y = "ASIR per 100,000",
       color = "location") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white"), 
    panel.background = element_rect(fill = "grey90"), 
    panel.grid.major = element_line(color = "white", size = 0.5),
    panel.grid.minor = element_line(color = "white", size = 0.3), 
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    axis.line = element_line(color = "black"),      
    axis.ticks = element_line(color = "black")     
  ) +
  guides(color = guide_legend(reverse = TRUE))  


#Random forest
library(randomForest)
set.seed(2025)
rf_model <- randomForest(Congenital_birth_defects ~ DryFrequency + DryDuration + DrySeverity + DryImpact_year_ +Temp + SDI_rank + Year,
                         data = data,
                         ntree = 500,  
                         mtry = 3,    
                         importance = TRUE, 
                         type = "regression")  
print(rf_model)
importance(rf_model)
drought_importance <- importance(rf_model)[c("DryFrequency", "DryDuration", "DrySeverity", "DryImpact_year_"), ]
drought_importance_df <- data.frame(
  Importance = drought_importance[, 1],  # 1(%IncMSE,Mean Decrease in Accuracy);2(IncNodePurity,Mean Decrease in Gini Impurity)
  Indicator = rownames(drought_importance)
)
ggplot(drought_importance_df, aes(x = reorder(Indicator, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() + 
  labs(
    title = "Mean Decrease in Accuracy",
    x = "Drought Indicators",
    y = "Importance Score") +
  theme_minimal()


#Forest plot
library(grid)
library(forestploter)
tm <- forest_theme(base_size = 10,
                   ci_pch =20,
                   ci_col = "#4575b4",#762a83
                   ci_lty = 1,
                   ci_lwd = 2,
                   ci_Theight = 0.2,
                   refline_lwd = 1.5,
                   refline_lty = "dashed",
                   refline_col = "red",
                   summary_fill  = "#4575b4",
                   summary_col = "#4575b4",
                   footnote_cex = 1.1,
                   footnote_fontface = "italic",
                   footnote_col = "blue")
p <- forest(data[,c(1,3:5)],
            est = data$RR,
            lower = data$LCI2,
            upper = data$UCI2,
            sizes = 0.6,
            ci_column = 3,
            ref_line = 1,
            xlim = c(1,1.0025),
            ticks_at = c(1,1.0025),
            theme = tm)
print(p)
