#### Data Exploration - hypergraph project thesis

rm(list=ls())

##Load required R packages
library(viridis)
library(igraph)
library(ggplot2)
library(tidyr)
library(dplyr)

###################################################
#Define parameters 
#Family, pairs and singles
###40% JUVENILES TO ADULTS -> ~1:2 PAIR TO FAMILY RATIO
Fs <- c(20, 10, 40)
Ps <- c(10, 5, 20)
#singles aren't really interesting to vary, just going to use one value
Ss <- 2
runs <- 5
space_runs <- 5
#spatial threshold (for proximity networks)
thresh<-c(0.1, 0.05, 0.15)
#Starting infectious dose
d_start<-0.025
#Additional responsiveness for sub-group and family ################ WILL USE THIS 
d_sub<-0.025
d_fam<-0.1
#Calling threshold
thresh_call<-c(0.2, 0.1, 0.25)
#Leaving threshold (departures)
leave_thresh<-0.9
#timestep counter
time_count <- 0


#Create empty data frame to store results
dyad_model_stats <- data.frame(replicate_no = integer(), space_replicate = integer(),
                          F = integer(), P = integer(), S = integer(), t=integer(), tc=integer(),
                          dyad_mean_leave_size = numeric(), dyad_sd_ls = numeric(), dyad_cv_ls = numeric(),
                          dyad_mean_leave_times = numeric())

hyp_model_stats <- data.frame(replicate_no = integer(), space_replicate = integer(),
                              F = integer(), P = integer(), S = integer(), t=integer(), tc=integer(),
                               hyp_mean_leave_size = numeric(), hyp_sd_ls = numeric(), hyp_cv_ls = numeric(),
                               hyp_mean_leave_times = numeric())

c_ts<-numeric()
dyad_c_wl<-numeric()
hyp_c_wl<-numeric()

#Setting up model functions
#statistics calculations for the dyadic models
calculate_dyad_stats <- function(dyad_c_wl=dyad_c_wl, F=F, P=P, S=S, t=t, tc=tc, replicate_no=runs, space_replicate=space_runs) {
  #Calculate leave sizes
  dyad_leave_size <- diff(dyad_c_wl, na.rm = TRUE)
  dyad_leave_size <- dyad_leave_size[dyad_leave_size > 0]
  
  #Calculate when groups left
  dyad_leave_times <- which(diff(dyad_c_wl, na.rm = TRUE) > 0)
  
  #Calculate descriptive stats
  dyad_mean_leave_size <- mean(dyad_leave_size, na.rm = TRUE)
  dyad_sd_ls <- sd(dyad_leave_size, na.rm = TRUE)
  dyad_cv_ls <- dyad_sd_ls / dyad_mean_leave_size
  dyad_mean_leave_times <- mean(dyad_leave_times, na.rm = TRUE)
  
  #Return a data frame with the calculated statistics
  return(data.frame(replicate_no = replicate_no, space_replicate=space_replicate,
                    F = F, P = P, S = S, t=t, tc=tc, dyad_mean_leave_size = dyad_mean_leave_size, dyad_sd_ls = dyad_sd_ls,
                    dyad_cv_ls = dyad_cv_ls, dyad_mean_leave_times = dyad_mean_leave_times))
  
}

#statistics calculations for the hypergraph models
calculate_hyp_stats <- function(hyp_c_wl=hyp_c_wl, F=F, P=P, S=S, t=t, tc=tc, replicate_no=runs, space_replicate=space_runs) {
  #Calculate leave sizes
  hyp_leave_size <- diff(hyp_c_wl, na.rm = TRUE)
  hyp_leave_size <- hyp_leave_size[hyp_leave_size > 0]
  
  #Calculate when groups left
  hyp_leave_times <- which(diff(hyp_c_wl, na.rm = TRUE) > 0)
  
  #Calculate descriptive stats
  hyp_mean_leave_size <- mean(hyp_leave_size, na.rm = TRUE)
  hyp_sd_ls <- sd(hyp_leave_size, na.rm = TRUE)
  hyp_cv_ls <- hyp_sd_ls / hyp_mean_leave_size
  hyp_mean_leave_times <- mean(hyp_leave_times, na.rm = TRUE)
  
  #Return a data frame with the calculated statistics
  return(data.frame(replicate_no = replicate_no, space_replicate=space_replicate,
                    F = F, P = P, S = S, t=t, tc=tc, hyp_mean_leave_size = hyp_mean_leave_size, hyp_sd_ls = hyp_sd_ls,
    hyp_cv_ls = hyp_cv_ls, hyp_mean_leave_times = hyp_mean_leave_times))
}

#Dose response function
dr_curve<-function(x){
  return(1/(1+exp(-10*(x-0.5)))-1/(1+exp(-10*(-0.5))))
}

net_inf<-function(x){
  return(rbinom(length(x),1,prob=x))
}


#########################################
#Running the network code
for (t in thresh) {
  for (tc in thresh_call) {
    print(paste("Thresh:", t, "Call thresh:", tc))
    
for (run in 1:runs) {
  #Print the current run
  print(paste("Run:", run))
  
#Loop through each combination
for (F in Fs) {
  for (P in Ps) {
    for (S in Ss) {
      
      #Print current iteration for reference
      print(paste("F:", F, "P:", P, "S:", S))

for (space_run in 1:space_runs){
  #print(paste("Within space run:", space_run))
    
##Generate individuals in group
n_families<- F
n_pairs<- S
n_singles<- P

#Calculate family size
family_size<-3+rpois(n_families,1)

#Work out number of individuals
n_indivs<-sum(family_size)+n_pairs*2+n_singles

#assign IDs
id<-seq(1,n_indivs,1)

#assign Family IDs
family<-c(rep(1:n_families,family_size),rep((n_families+1):(n_families+n_pairs),each=2),
          (n_families+n_pairs+1):(n_families+n_pairs+n_singles))

#Assign Family status
status<-c(rep("fam",sum(family_size)),rep("pair",n_pairs*2),rep("sing",n_singles))

#Create individual dataframe
ind_data<-data.frame(id,family,status)

###################################################

##Assign locations in group

#Family coordinates
family_x<-runif(n_families+n_pairs+n_singles,0,(n_families*4+n_pairs*2+n_singles)/100) #might still need to fix this
family_y<-runif(n_families+n_pairs+n_singles,0,(n_families*4+n_pairs*2+n_singles)/100)

#Assign individual coordinates
ind_x<-ind_y<-rep(NA,nrow(ind_data))
for(i in 1:nrow(ind_data)){
  ind_x[i]<-rnorm(1,family_x[ind_data$family[i]],0.05)
  ind_y[i]<-rnorm(1,family_y[ind_data$family[i]],0.05)
}

###################################################

##Create family network and hypergraph
fam_net<-matrix(0,nr=n_indivs,nc=n_indivs)
fam_hyp<-matrix(0,nr=n_families+n_pairs+n_singles,nc=n_indivs)

for(i in 1:(n_indivs-1)){
  for(j in (i+1):n_indivs){
    if(ind_data$family[i]==ind_data$family[j]){fam_net[i,j]<-fam_net[j,i]<-1}
  }
}

#Create hypergraph incidence matrix
for(i in 1:(n_families+n_pairs+n_singles)){
  fam_hyp[i,which(ind_data$family==i)]<-1
}

###################################################

##Create spatial network and hypergraph

#Calculate distance matrix
dist_mat<-as.matrix(dist(cbind(ind_x,ind_y)))
diag(dist_mat)<-NA

#Create matrix of subgroup membership
#sp_mat<-dist_mat<thresh[t] 
sp_mat<-dist_mat<t ###NEEDS FIXING - think this one might work

sp_net<-apply(sp_mat,2,as.numeric)

###############################################
#Dyadic call transmission using dose response curve

call_transmit_net<-function(fam_net,sp_net,d_start,d_sub,d_fam,callers,left=NULL){
  
  dose_mat<-call_net[callers,]
  fam_mat<-dose_mat*fam_net[callers,]
  sg_mat<-dose_mat*sp_net2[callers,]
  
  dose_mat2<-dose_mat*d_start+fam_mat*d_fam+sg_mat*d_sub
  if(length(left)>0){dose_mat2[left,]<-0}
  dose_mat2[,callers]<-0
  
  transmit_mat<-t(apply(dose_mat2,1,dr_curve))
  
  nc_mat<-t(apply(transmit_mat,1,net_inf))
  
  nc<-sign(colSums(nc_mat,na.rm=T))
  
  new_callers<-which(nc>0)
  new_callers
  return(new_callers)
  
}

###############################################

#Find components of network
components<-decompose(graph_from_adjacency_matrix(sp_net,mode="undirected"))

#Use to create network and hypergraph
ind_components<-rep(NA,nrow(ind_data))
for(i in 1:length(components)){
  ind_components[as.numeric(V(components[[i]])$name)]<-i
}
sp_net2<-matrix(0,nr=n_indivs,nc=n_indivs) #dyadic
sp_hyp<-matrix(0,nr=length(components),nc=n_indivs) #hypergraph

for(i in 1:(n_indivs-1)){
  for(j in (i+1):n_indivs){
    if(ind_components[i]==ind_components[j]){sp_net2[i,j]<-sp_net2[j,i]<-1}
  }
}

#Create hypergraph
for(i in 1:length(components)){
  sp_hyp[i,which(ind_components==i)]<-1
}

###################################################

##Create an additional network for being able to hear calls
#sp_mat<-dist_mat<thresh_call ####NEEDS FIXING
sp_mat<-dist_mat<tc ###AGAIN THINK THIS MIGHT BE THE ONE
call_net<-apply(sp_mat,2,as.numeric)

###################################################
##Dosage transmission of calls with departures
#Seed calling individuals
callers<-sample(ind_data$id,2)
og_callers<-callers

left<-numeric()

call_net2<-call_net

#Model for 1000 timesteps
for(i in 1:1000){
  new_callers<-call_transmit_net(fam_net=fam_net,sp_net=sp_net2,d_start=d_start,d_sub=d_sub,d_fam=d_fam,callers=callers)
  if(length(new_callers)>0){callers<-sort(c(callers,new_callers))}
  c_ts[i]<-length(callers)
  callers2<-rep(0,nrow(ind_data))
  callers2[callers]<-1
  calling_net<-call_net2*callers2
  if(length(left)>0){
    call_net2[left,]<-call_net2[,left]<-0
    calling_net[left,]<-calling_net[,left]<-0
  }
  prop_calling<-colSums(calling_net,na.rm=T)/colSums(call_net2,na.rm=T)
  if(length(left)==0){
    left<-which(as.numeric(prop_calling>leave_thresh)>0)
  }
  if(length(left)>0){
    nl<-which(as.numeric(prop_calling>leave_thresh)>0)
    left<-sort(unique(c(left,nl))) ##ADDED THE UNIQUE
  }
  dyad_c_wl[i]<-length(left)
  if(length(callers)==nrow(ind_data)){break()}
  if(length(callers)<= nrow(ind_data)) {
    time_count <- time_count + 1
    if(time_count >= 100) {break()}
  } 
  else {time_count <- 0}
  #Resetting counter for next run
  time_count<-0
}

#run the calculation function
dyad_stats <- calculate_dyad_stats(dyad_c_wl=dyad_c_wl, F=F, P=P, S=S, t=t, tc=tc, 
                                   replicate_no=run, space_replicate=space_run)
dyad_model_stats <- rbind(dyad_model_stats, dyad_stats)



################################
# Repeating transmission/dose functions now applied to the hypergraph networks 

call_transmit_hyp<-function(fam_hyp,sp_hyp,call_net,d_start,d_sub,d_fam,callers,left=NULL){
  
  callers2<-rep(0,ncol(fam_hyp))
  callers2[callers]<-1
  
  fam_dose<-numeric()
  sg_dose<-numeric()
  base_dose<-numeric()
  
  for(i in 1:length(callers2)){
    if(callers2[i]==0){
      t_vec1<-fam_hyp[fam_hyp[,i]==1,]
      t_vec1[i]<-0
      if(length(left)>0){t_vec1[left]<-0}
      fam_dose[i]<-sum(callers2*t_vec1)*d_fam
      t_vec2<-sp_hyp[sp_hyp[,i]==1,]
      t_vec2[i]<-0
      if(length(left)>0){t_vec2[left]<-0}
      sg_dose[i]<-sum(callers2*t_vec2)*d_sub
      t_vec3<-call_net[i,]
      t_vec3[i]<-0
      if(length(left)>0){t_vec3[left]<-0}
      base_dose[i]<-sum(callers2*t_vec3)*d_start
    }
  }
  
  dose<-fam_dose+sg_dose+base_dose
  transmit<-sapply(dose,dr_curve)
  nc<-rbinom(length(transmit), 1, prob = transmit)
  
  new_callers<-which(nc>0)
  new_callers
  return(new_callers)
  
}
###################################################

##Test hypergraph model with departures
#Seed calling individuals

#callers<-sample(ind_data$id,2)
callers<-og_callers

left<-numeric()

call_net2<-call_net

#Model over 1000 timesteps
for(i in 1:1000){
  new_callers<-call_transmit_hyp(fam_hyp=fam_hyp,sp_hyp=sp_hyp,call_net=call_net,d_start=d_start,d_sub=d_sub,d_fam=d_fam,callers=callers,left=left)
  if(length(new_callers)>0){callers<-sort(c(callers,new_callers))}
  c_ts[i]<-length(callers)
  callers2<-rep(0,nrow(ind_data))
  callers2[callers]<-1
  calling_net<-call_net2*callers2
  if(length(left)>0){
    call_net2[left,]<-call_net2[,left]<-0
    calling_net[left,]<-calling_net[,left]<-0
  }
  prop_calling<-colSums(calling_net,na.rm=T)/colSums(call_net2,na.rm=T)
  if(length(left)==0){
    left<-which(as.numeric(prop_calling>leave_thresh)>0)
  }
  if(length(left)>0){
    nl<-which(as.numeric(prop_calling>leave_thresh)>0)
    left<-sort(unique(c(left,nl))) ##ADDED THE UNIQUE
  }
  hyp_c_wl[i]<-length(left)
  if(length(callers)==nrow(ind_data)){break()}
  if(length(callers)<= nrow(ind_data)) {
    time_count <- time_count + 1
    if(time_count >= 100) {break()}
  } 
  else {time_count <- 0}
  time_count<-0
}


#Run the calculation function
  hyp_stats <- calculate_hyp_stats(hyp_c_wl = hyp_c_wl, F=F, P=P, S=S, t=t, tc=tc,
                                   replicate_no=run, space_replicate=space_run)
  hyp_model_stats <- rbind(hyp_model_stats, hyp_stats)

}

     }
    }
   }
  }
 }
}
############################

#Replace NA values with 0s
dyad_model_stats[is.na(dyad_model_stats)] <- 0
hyp_model_stats[is.na(hyp_model_stats)] <- 0

combined_dyad_stats <- aggregate(cbind(dyad_mean_leave_size, dyad_sd_ls, dyad_cv_ls, dyad_mean_leave_times) ~ 
                                   F + P + S + t + tc + replicate_no + space_replicate,
                                 data = dyad_model_stats, 
                                 FUN = mean, na.rm = TRUE)

combined_hyp_stats <- aggregate(cbind(hyp_mean_leave_size, hyp_sd_ls, hyp_cv_ls, hyp_mean_leave_times) ~ 
                                  F + P + S + t + tc + replicate_no + space_replicate,
                                data = hyp_model_stats, 
                                FUN = mean, na.rm = TRUE)

#Merge dyadic and hypergraph statistics
all_model_stats <- merge(combined_dyad_stats, combined_hyp_stats, 
                         by = c("F", "P", "S", "t", "tc", "replicate_no", "space_replicate"), 
                         suffixes = c("_dyad", "_hyp"), all = TRUE)

#creating new column for dyad and hyp comparison
all_model_stats$comb_leave_size <- all_model_stats$hyp_mean_leave_size - all_model_stats$dyad_mean_leave_size


#merging F P and S into one column
all_model_stats$fam_combo <- paste0(all_model_stats$F, "_", 
                                    all_model_stats$P, "_",
                                    all_model_stats$S)

all_model_stats$replicate_unique <- paste0(all_model_stats$F,"_",all_model_stats$P,"_",all_model_stats$S,
                                           "_",all_model_stats$t,"_",all_model_stats$tc,"_",all_model_stats$replicate_no)
#all_model_stats <- all_model_stats[, c(1:3, ncol(all_model_stats), 4:(ncol(all_model_stats)-1))]


#saving the dataset
#write.csv(all_model_stats, file = "gby_data_nosingle_final.csv", row.names = FALSE)

#############
all_model_stats<-read.csv("gby_data_nosingle_final.csv")

library(lme4)

str(all_model_stats)

all_model_stats$t <- as.factor(all_model_stats$t)
all_model_stats$tc <- as.factor(all_model_stats$tc)

#hypergraph variance -> random effects model
hyp_variance_model <- lmer(hyp_mean_leave_size ~ 1 + (1|fam_combo) + (1|t) + (1|tc) + (1|replicate_unique), data = all_model_stats)
summary(hyp_variance_model)

#dyadic variance -> random effects model
dyad_variance_model <- lmer(dyad_mean_leave_size ~ 1 + (1|fam_combo) + (1|t) + (1|tc) + (1|replicate_unique), data = all_model_stats)
summary(dyad_variance_model)

#############
#Mixed effects models
#install.packages('MuMIn')
library(MuMIn)

hyp_mod2<-lmer(hyp_mean_leave_size ~ as.factor(fam_combo) + as.factor(t) + as.factor(tc) + (1|replicate_unique), data = all_model_stats)
summary(hyp_mod2)
#calculating the replicate variance
hyp_rep<-r.squaredGLMM(hyp_mod2)
re_var_hyp<-hyp_rep[2]-hyp_rep[1]

dyad_mod2<-lmer(dyad_mean_leave_size ~ as.factor(fam_combo) + as.factor(t) + as.factor(tc) + (1|replicate_unique), data = all_model_stats)
summary(dyad_mod2)
dyad_rep<-r.squaredGLMM(dyad_mod2)
re_var_dyad<-dyad_rep[2]-dyad_rep[1]

#PartR2 models
#install.packages("partR2")
library(partR2)

#hyp
hyp_r2<-partR2(hyp_mod2,partvars=c("as.factor(fam_combo)","as.factor(t)","as.factor(tc)"))
#str(hyp_r2)
hyp_r2$R2[,1:2]

#dyad
dyad_r2<-partR2(dyad_mod2,partvars=c("as.factor(fam_combo)","as.factor(t)","as.factor(tc)"))
#str(dyad_r2)
dyad_r2$R2[,1:2]


##############
#plotting the variance models

# Manually creating the data frame with the R2 values
#DYADIC
dyad_partR2_data <- data.frame(
  Predictor = c("Model", "as.factor(fam_combo)", "as.factor(t)", "as.factor(tc)", "as.factor(fam_combo)+as.factor(t)", 
                "as.factor(fam_combo)+as.factor(tc)", "as.factor(t)+as.factor(tc)", "as.factor(fam_combo)+as.factor(t)+as.factor(tc)"),
  R2 = c(0.583, 0.189, 0.00863, 0.385, 0.198, 0.574, 0.393, 0.583)
)

# Extract the necessary R2 values for the plot
dyad_var_data <- dyad_partR2_data %>% 
  filter(Predictor %in% c("as.factor(fam_combo)", "as.factor(t)", "as.factor(tc)")) %>%
  mutate(Predictor = factor(Predictor, levels = c("as.factor(fam_combo)", "as.factor(t)", "as.factor(tc)")))

# Calculate the residual variance
dyad_total_R2 <- dyad_partR2_data %>% filter(Predictor == "Model") %>% pull(R2)
dyad_residual_R2 <- 1 - dyad_total_R2

# Add the residual variance to the data frame
dyad_var_data <- rbind(dyad_var_data, data.frame(Predictor = "Residual", R2 = dyad_residual_R2))

#first the fixed effects (may need to change the 2:4)
dyad_var_data2<-dyad_r2$R2[2:4,1:2]
re_name<-'replicate_unique'
re_var2<-data.frame(re_name,re_var_dyad)
names(re_var2)<-names(dyad_var_data2)
dyad_var_data2<-rbind(dyad_var_data2,re_var2)
dyad_resid_hm<-1-sum(dyad_var_data2$estimate)
res_name<-'residual'
dyad_res_var2<-data.frame(res_name,dyad_resid_hm)
names(dyad_res_var2)<-names(dyad_var_data2)
dyad_var_data2<-rbind(dyad_var_data2,dyad_res_var2)
names(dyad_var_data2)<-names(dyad_var_data)

#########
#HYPERGRAPH
hyp_partR2_data <- data.frame(
  Predictor = c("Model", "as.factor(fam_combo)", "as.factor(t)", "as.factor(tc)", "as.factor(fam_combo)+as.factor(t)", 
                "as.factor(fam_combo)+as.factor(tc)", "as.factor(t)+as.factor(tc)", "as.factor(fam_combo)+as.factor(t)+as.factor(tc)"),
  R2 = c(0.366, 0.0931, 0.255, 0.0178, 0.348, 0.111, 0.273, 0.366)
)

# Extract the necessary R2 values for the plot
hyp_var_data <- hyp_partR2_data %>% 
  filter(Predictor %in% c("as.factor(fam_combo)", "as.factor(t)", "as.factor(tc)")) %>%
  mutate(Predictor = factor(Predictor, levels = c("as.factor(fam_combo)", "as.factor(t)", "as.factor(tc)")))

# Calculate the residual variance
hyp_total_R2 <- hyp_partR2_data %>% filter(Predictor == "Model") %>% pull(R2)
hyp_residual_R2 <- 1 - hyp_total_R2

# Add the residual variance to the data frame
hyp_var_data <- rbind(hyp_var_data, data.frame(Predictor = "Residual", R2 = hyp_residual_R2))

#first the fixed effects (may need to change the 2:4)
hyp_var_data2<-hyp_r2$R2[2:4,1:2]
re_name<-'replicate_unique'
re_var2<-data.frame(re_name,re_var_hyp)
names(re_var2)<-names(hyp_var_data2)
hyp_var_data2<-rbind(hyp_var_data2,re_var2)
hyp_resid_hm<-1-sum(hyp_var_data2$estimate)
res_name<-'residual'
hyp_res_var2<-data.frame(res_name,hyp_resid_hm)
names(hyp_res_var2)<-names(hyp_var_data2)
hyp_var_data2<-rbind(hyp_var_data2,hyp_res_var2)
names(hyp_var_data2)<-names(hyp_var_data)

library(gridExtra)
library(paletteer)
# Plot the data
dyad_plot<- ggplot(dyad_var_data2, aes(x = 1, y = R2, fill = Predictor)) +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_continuous(labels = NULL) +
  labs(x = "Dyadic", y = "Proportion of Variance Explained") +
  theme(axis.ticks.x = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 14),   
        axis.title = element_text(size = 16) 
  ) +
  scale_fill_paletteer_d("nationalparkcolors::Acadia") +
  guides(fill = guide_legend(title = "Predictor"))

hyp_plot<- ggplot(hyp_var_data2, aes(x = 1, y = R2, fill = Predictor)) +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_continuous(labels = NULL) +
  labs(x = "Hypergraph") +
  theme(axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 14),  
        axis.title = element_text(size = 16)
  ) +
  scale_fill_paletteer_d("nationalparkcolors::Acadia") +
  guides(fill = guide_legend(title = "Predictor"))

grid.arrange(dyad_plot, hyp_plot, ncol = 2, widths = c(0.55, 1))


##########################
library(tidyr)

#making longform data
long_data <- all_model_stats %>%
  gather(key = "Measure", value = "MeanLeaveSize", dyad_mean_leave_size, hyp_mean_leave_size) %>%
  mutate(Measure = recode(Measure, 
                          dyad_mean_leave_size = "Dyadic", 
                          hyp_mean_leave_size = "Hypergraph"))

#family combo
ggplot(long_data, aes(x = as.factor(fam_combo), y = MeanLeaveSize, fill = Measure)) +
  geom_boxplot() +
  labs(x = "Family Combo", y = "Mean Leave Size", fill = "Measure") +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12),
  ) +
  scale_fill_paletteer_d("nationalparkcolors::Acadia")

#proximity threshold 
ggplot(long_data, aes(x = as.factor(t), y = MeanLeaveSize, fill = Measure)) +
  geom_boxplot() +
  labs(x = "Proximity Threshold", y = "Mean Leave Size", fill = "Measure") +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12),
  ) +
  scale_fill_paletteer_d("nationalparkcolors::Acadia")

#calling threshold
ggplot(long_data, aes(x = as.factor(tc), y = MeanLeaveSize, fill = Measure)) +
  geom_boxplot() +
  labs(x = "Calling Threshold", y = "Mean Leave Size", fill = "Measure") +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12),
  ) +
  scale_fill_paletteer_d("nationalparkcolors::Acadia")


#################################
#HISTOGRAMS FOR THESIS

acadia_colors <- as.character(paletteer::paletteer_d("nationalparkcolors::Acadia", n = 3))

# Determine the range of data for both variables
max_leave_size <- max(all_model_stats$dyad_mean_leave_size, all_model_stats$hyp_mean_leave_size, na.rm = TRUE)

# Create breaks that cover the entire range
breaks_seq <- seq(0, max_leave_size, 0.5)

# Plot the histogram for dyad_mean_leave_size
hist(all_model_stats$dyad_mean_leave_size, 
     col = adjustcolor(acadia_colors[3], 0.60), 
     border = NA, 
     breaks = breaks_seq, 
     ylim = c(0, 400), 
     xlim = c(0, 60),
     xlab = "Mean Leave Size",
     ylab = "Frequency",
     cex.main = 1.5, 
     cex.lab = 1.2, 
     cex.axis = 1.1,
     yaxt = "n")  # Suppress y-axis

# Custom y-axis with rotated labels
axis(2, las = 2)

# Adding the histogram for hyp_mean_leave_size
hist(all_model_stats$hyp_mean_leave_size, 
     col = adjustcolor(acadia_colors[2], 0.60), 
     border = NA, 
     breaks = breaks_seq, 
     add = TRUE)

legend("topright", legend = c("Dyadic", "Hypergraph"), 
       fill = c(acadia_colors[3], acadia_colors[2]), border = NA, cex = 1.2)


##########################
### Descriptive statistics 
#install.packages("sjPlot")
library(sjPlot)


# Reshape data to long format
long_data_mls <- all_model_stats %>%
  select(replicate_unique, dyad_mean_leave_size, hyp_mean_leave_size) %>%
  pivot_longer(cols = c(dyad_mean_leave_size, hyp_mean_leave_size),
               names_to = "Measure",
               values_to = "MeanLeaveSize")

# Convert Measure to a factor
long_data_mls$Measure <- factor(long_data_mls$Measure, levels = c("dyad_mean_leave_size", "hyp_mean_leave_size"))

# Mixed-effects model comparing dyadic and hypergraph mean leave sizes
compare_model_mls <- lmer(MeanLeaveSize ~ Measure + (1|replicate_unique), data = long_data_mls)

# Summary of the model
summary(compare_model_mls)

# Visualize the mixed-effects model -> sjPlot
plot_model(compare_model_mls, type = "pred", terms = "Measure") +
  ggtitle("Comparison of Dyadic and Hypergraph Mean Leave Size") +
  xlab("Measure") +
  ylab("Predicted Mean Leave Size")

tab_model(compare_model_mls, show.re.var = TRUE, show.icc = TRUE, show.r2 = TRUE)
tab_model(compare_model_mls, file = "baseline_model_summary_mls.html")