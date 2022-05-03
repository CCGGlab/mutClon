# Load overview donor data
###########################

# Note that these cannot be shared for privacy reasons
PM_data<- readRDS("raw/donors/anatomy_database/donor_data.rds")

# How many donors/year since overall registration started?
###########################################################

table(PM_data$year)
# 2016 2017 2018 2019 2020 2021 
# 108   82  100  133  102   84 

summary(as.numeric(table(PM_data$year)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 82.0    88.0   101.0   101.5   106.5   133.0 

summary(PM_data$Leeftijd)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   36.00   72.00   81.00   78.96   88.00  106.00      79 

# Distribution of PMI (only available for recent 2 years)?
##############################################################

PM_data_temp<- PM_data[!is.na(PM_data$PMI),]

summary(as.numeric(PM_data_temp$PMI))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.5167   19.3750   27.0833   56.9527   46.1458 2740.1667 

mean(as.numeric(PM_data_temp$PMI)>12&as.numeric(PM_data_temp$PMI)<=36) # 0.53
mean(as.numeric(PM_data_temp$PMI)<=72) # 0.9175824

# Plot
p_PMI<- ggplot(PM_data_temp, aes(x=PMI)) + 
  geom_histogram(aes(y = stat(count) / sum(count)), colour="white", fill="grey", binwidth = 6, center = 3)+
  geom_density(aes(y = 50 * stat(count) / sum(count)), alpha=.1, fill="#FF6666") +
  xlab("Post-mortem interval (h)") +
  ylab("") +
  scale_x_continuous(breaks = c(0,24,48,72), limits=c(0,72), expand = c(0, 0) ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0), limits=c(0,0.24)) +
  theme(axis.text = element_text(size=6),  
        axis.title = element_text(size=7),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# Age distribution of these donors?
#####################################

# Plot
p_age<- ggplot(PM_data_temp, aes(x=Leeftijd)) + 
  geom_histogram(aes(y = stat(count) / sum(count)), colour="white", fill="grey", binwidth = 5, center = 2.5)+
  geom_density(aes(y =45 * stat(count) / sum(count)), alpha=.1, fill="#FF6666") +
  xlab("Age (y)") +
  ylab("") +
  scale_x_continuous(breaks = seq(0, 110, 10), expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0), limits=c(0,0.24) ) +
  theme(axis.text = element_text(size=6),  
        axis.title = element_text(size=7),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# Save plots
##############
p<- grid.arrange(
  p_age,
  p_PMI,
  ncol=2
)
ggsave("results/figs/manuscript_donor_description.pdf",  p, width = 75, height = 40, units = "mm")

