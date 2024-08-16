#### Data Exploration - hypergraph project thesis

### Kinship models - increasing family subgroup responsiveness 
rm(list=ls())

##Load required R packages
library(viridis)
library(igraph)
library(ggplot2)
library(tidyr)

###################################################
#Define parameters 
#Family, pairs and singles
Fs <- 20
Ps <- 10
Ss <- 2
runs <- 100 #will adjust later
space_runs <- 5 #will adjust later
#spatial threshold (for proximity networks)
thresh<-c(0.1, 0.05, 0.15)
#Starting infectious dose
d_start<-0.025
#Additional responsiveness for sub-group and family ################ WILL USE THIS 
d_sub<-0.0125
d_fam<-0.2
#Calling threshold
thresh_call<-0.2
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
                    F = F, P = P, S = S, t=t, tc=tc,
                    dyad_mean_leave_size = dyad_mean_leave_size, dyad_sd_ls = dyad_sd_ls,
                    dyad_cv_ls = dyad_cv_ls, dyad_mean_leave_times = dyad_mean_leave_times))
  
}

#statistics calculations for the hypergraph models
calculate_hyp_stats <- function(hyp_c_wl=hyp_c_wl, F=F, P=P, S=S, t=t, tc=tc,
                                replicate_no=runs, space_replicate=space_runs) {
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
                    F = F, P = P, S = S, t=t, tc=tc,
                    hyp_mean_leave_size = hyp_mean_leave_size, hyp_sd_ls = hyp_sd_ls,
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
#multipliers = mean from random distribution
family_x<-runif(n_families+n_pairs+n_singles,0,(n_families*4+n_pairs*2+n_singles)/100) 
family_y<-runif(n_families+n_pairs+n_singles,0,(n_families*4+n_pairs*2+n_singles)/100)

#Assign individual coordinates
ind_x<-ind_y<-rep(NA,nrow(ind_data))
for(i in 1:nrow(ind_data)){
  ind_x[i]<-rnorm(1,family_x[ind_data$family[i]],0.05)
  ind_y[i]<-rnorm(1,family_y[ind_data$family[i]],0.05)
}

#Assign colours to families for a plot
fam_cols<-sample(viridis(n_families+n_pairs+n_singles),n_families+n_pairs+n_singles,replace=FALSE)

#Plot individual locations
plot(ind_x,ind_y,col=fam_cols[ind_data$family],pch=16,cex=2)

###################################################

##Create family network and hypergraph
fam_net<-matrix(0,nr=n_indivs,nc=n_indivs)
fam_hyp<-matrix(0,nr=n_families+n_pairs+n_singles,nc=n_indivs)

for(i in 1:(n_indivs-1)){
  for(j in (i+1):n_indivs){
    if(ind_data$family[i]==ind_data$family[j]){fam_net[i,j]<-fam_net[j,i]<-1}
  }
}

#Plot family network
#plot(graph_from_adjacency_matrix(fam_net))

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
sp_mat<-dist_mat<thresh
sp_net<-apply(sp_mat,2,as.numeric)

#Check this has worked by plotting netwok
#plot(graph_from_adjacency_matrix(sp_net,mode="undirected"),vertex.label=NA,vertex.size=8) #dyadic

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
sp_mat<-dist_mat<thresh_call
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
dyad_stats <- calculate_dyad_stats(dyad_c_wl=dyad_c_wl, F=F, P=P, S=S, replicate_no=run, space_replicate=space_run,
                                   t=t, tc=tc)
dyad_model_stats <- rbind(dyad_model_stats, dyad_stats)

#Plot results
#par(mfrow=c(1,1))
#plot(c_ts,type="l",col="#0000FF",ylim=c(0,nrow(ind_data)))
#lines(dyad_c_wl,col="#FFCC33")

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
hyp_stats <- calculate_hyp_stats(hyp_c_wl = hyp_c_wl, F=F, P=P, S=S,replicate_no=run, space_replicate=space_run,
                                 t=t, tc=tc)
hyp_model_stats <- rbind(hyp_model_stats, hyp_stats)
  
#Plot results
#par(mfrow=c(1,1))
#plot(c_ts,type="l",col="#0000FF",ylim=c(0,nrow(ind_data)))
#lines(hyp_c_wl,col="#FFCC33")

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
write.csv(all_model_stats, file = "kinship_model_replicated_3aug.csv", row.names = FALSE)
all_model_stats <- read.csv("kinship_model_replicated_3aug.csv")

#################################
#### MEAN LEAVE SIZE ANALYSIS 
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

#########################
#MEAN LEAVE TIMES PLOTS

#making longform data
long_data_mlt <- all_model_stats %>%
  gather(key = "Measure", value = "MeanLeaveTimes", dyad_mean_leave_times, hyp_mean_leave_times) %>%
  mutate(Measure = recode(Measure, 
                          dyad_mean_leave_times = "Dyadic", 
                          hyp_mean_leave_times = "Hypergraph"))

#family combo
ggplot(long_data_mlt, aes(x = as.factor(fam_combo), y = MeanLeaveTimes, fill = Measure)) +
  geom_boxplot() +
  labs(x = "Family Combo", y = "Mean Leave Times", fill = "Measure") +
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
ggplot(long_data_mlt, aes(x = as.factor(t), y = MeanLeaveTimes, fill = Measure)) +
  geom_boxplot() +
  labs(x = "Proximity Threshold", y = "Mean Leave Times", fill = "Measure") +
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
ggplot(long_data_mlt, aes(x = as.factor(tc), y = MeanLeaveTimes, fill = Measure)) +
  geom_boxplot() +
  labs(x = "Calling Threshold", y = "Mean Leave Times", fill = "Measure") +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12),
  ) +
  scale_fill_paletteer_d("nationalparkcolors::Acadia")

###################
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
     ylim = c(0, 600), 
     xlab = "Mean Leave Size",
     ylab = "Frequency")

# Adding the histogram for hyp_mean_leave_size
hist(all_model_stats$hyp_mean_leave_size, 
     col = adjustcolor(acadia_colors[2], 0.60), 
     border = NA, 
     breaks = breaks_seq, 
     add = TRUE)

################
# Determine the range of data for both variables
max_leave_times <- max(all_model_stats$dyad_mean_leave_times, all_model_stats$hyp_mean_leave_times, na.rm = TRUE)
min_leave_times <- min(all_model_stats$dyad_mean_leave_times, all_model_stats$hyp_mean_leave_times, na.rm = TRUE)

# Extend the range slightly to ensure coverage
breaks_seq_mlt <- seq(min_leave_times - 0.5, max_leave_times + 0.5, 0.5)

# Plot the histogram for dyad_mean_leave_size
hist(all_model_stats$dyad_mean_leave_size, 
     col = adjustcolor(acadia_colors[3], 0.60), 
     border = NA, 
     breaks = breaks_seq, 
     ylim = c(0, 600), 
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


######################
#mean leave times

breaks_seq_mlt2 <- seq(min_leave_times - 5, max_leave_times + 5, 10)

hist(all_model_stats$dyad_mean_leave_times, 
     col = adjustcolor(acadia_colors[3], 0.60), 
     border = NA, 
     breaks = breaks_seq_mlt2, 
     ylim = c(0, 150), 
     xlab = "Mean Leave Times (timesteps)",
     ylab = "Frequency",
     cex.main = 1.5, 
     cex.lab = 1.2, 
     cex.axis = 1.1,
     yaxt = "n")  # Suppress y-axis
# Custom y-axis with rotated labels
axis(2, las = 2)

# Adding the histogram for hyp_mean_leave_size
hist(all_model_stats$hyp_mean_leave_times, 
     col = adjustcolor(acadia_colors[2], 0.60), 
     border = NA, 
     breaks = breaks_seq_mlt2, 
     add = TRUE)

legend("topright", legend = c("Dyadic", "Hypergraph"), 
       fill = c(acadia_colors[3], acadia_colors[2]), border = NA, cex = 1.2)



### DESCRIPTIVE STATS

# Reshape data to long format
#install.packages("sjPlot")
library(sjPlot)
library(lme4)
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
tab_model(compare_model_mls, file = "social_model_summary_mls.html")

# Plot the model predictions
ggplot(long_data_mls, aes(x = Measure, y = MeanLeaveSize, fill = Measure)) +
  geom_boxplot() +
  labs(x = "Network Type", y = "Mean Leave Size", fill = "Measure") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),  
    axis.title = element_text(size = 16),
    legend.position = "none"
  ) +
  scale_fill_paletteer_d("nationalparkcolors::Acadia")

###########
#MLT
### Descriptive statistics 

# Reshape data to long format
long_data_mlt <- all_model_stats %>%
  select(replicate_unique, dyad_mean_leave_times, hyp_mean_leave_times) %>%
  pivot_longer(cols = c(dyad_mean_leave_times, hyp_mean_leave_times),
               names_to = "Measure",
               values_to = "MeanLeaveTimes")

# Convert Measure to a factor
long_data_mlt$Measure <- factor(long_data_mlt$Measure, levels = c("dyad_mean_leave_times", "hyp_mean_leave_times"))

# Mixed-effects model comparing dyadic and hypergraph mean leave sizes
compare_model_mlt <- lmer(MeanLeaveTimes ~ Measure + (1|replicate_unique), data = long_data_mlt)

# Summary of the model
summary(compare_model_mlt)

# Visualize the mixed-effects model -> sjPlot
plot_model(compare_model_mlt, type = "pred", terms = "Measure") +
  ggtitle("Comparison of Dyadic and Hypergraph Mean Leave Times") +
  xlab("Measure") +
  ylab("Predicted Mean Leave Times")

tab_model(compare_model_mlt, show.re.var = TRUE, show.icc = TRUE, show.r2 = TRUE)
tab_model(compare_model_mlt, file = "social_model_summary_mlt.html")

# Plot the model predictions
ggplot(long_data_mlt, aes(x = Measure, y = MeanLeaveTimes, fill = Measure)) +
  geom_boxplot() +
  labs(x = "Measure", y = "Mean Leave Times", fill = "Measure") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),  
    axis.title = element_text(size = 16),
    legend.position = "none"
  ) +
  scale_fill_paletteer_d("nationalparkcolors::Acadia")

