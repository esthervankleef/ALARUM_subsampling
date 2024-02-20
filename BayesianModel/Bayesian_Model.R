# replace "~/Desktop/Shared Code for Soon" by "folder-location"
#location_of_readme <- "~/Desktop/folder-location"
location_of_readme<- "C:/Users/evankleef/OneDrive - ITG/Documenten/GitHub/ITG/Other_projects/AMR/ALARUM/Analyses/ALARUM_subsampling/BayesianModel"

# load libraries 
library(tidyverse)
library(rstan)
library(bayesplot)
library(loo)
library(bayestestR)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
n_chains <- 4

# load data  
#dd <- read_csv(file=str_c(location_of_readme,'/data/Dataset_For_Bayesian_Model.csv') ) %>% 
dd <- read_csv(file=str_c(location_of_readme,'/Data_For_Model/Dataset_For_Bayesian_Model.csv') ) %>% 
                select(Antibiotic,Setting,
                       Rcgc_1_2,
                       Rcgc_2,
                       Rcgc2_zero_replaced,
                       Rtax_colkleebsal_Bracken,
                       Rtax_f_E_Bracken,
                       Rtax_o_E_Bracken,
                       Res_Inf_Ent,
                       Tot_Tested_Ent)

#  split column Antibiotic
dd1 <- dd %>% separate(Antibiotic, into=c("Antibiotic","count"),sep = " " ) %>% select(-count)
# rename cols
names(dd1)[3:8] <- c("ar12","ar2","ar1o2","b4div","entdiv","entbacdiv")

n_atbs <- dd1 %>% group_by(Antibiotic) %>% n_groups() # 16
#  extract infection data
df_transf2 <- dd1 %>% select(Setting,Antibiotic,Tot_Tested_Ent) %>% 
                spread(key=Antibiotic ,value=Tot_Tested_Ent)
df_transf3 <- dd1 %>% select(Setting,Antibiotic,Res_Inf_Ent) %>% 
                spread(key=Antibiotic ,value=Res_Inf_Ent)
#
setting_v <- df_transf2 %>% .$Setting
count_s <- df_transf2  %>% select(-Setting) %>% as.matrix() 
count_sr <- df_transf3  %>% select(-Setting) %>% as.matrix()
data_exists <- replace(count_sr, !is.na(count_sr), 1 ) # make helper that says where data is
data_exists <- replace(data_exists , is.na(data_exists), 0 ) # and where it isn't
count_s_hna <- replace(count_s,is.na(count_s),0) # hide NAs with 0
count_sr_hna <- replace(count_sr,is.na(count_sr),0) # hide NAs with 0

# taxonomic predictors (standadise)
df_transf <- dd1 %>% select(Setting,Antibiotic,b4div) %>% 
                spread(key=Antibiotic ,value=b4div)
z <- df_transf %>% select(-Setting) %>% as.matrix()
z_1a <- (z - mean(z)) / sd(z) # standardise
#
df_transf <- dd1 %>% select(Setting,Antibiotic,entdiv) %>% 
                spread(key=Antibiotic ,value=entdiv)
z <- df_transf %>% select(-Setting) %>% as.matrix()
z_1b <- (z - mean(z)) / sd(z) # standardise
#
df_transf <- dd1 %>% select(Setting,Antibiotic,entbacdiv) %>% 
                spread(key=Antibiotic ,value=entbacdiv)
z <- df_transf %>% select(-Setting) %>% as.matrix()
z_1c <- (z - mean(z)) / sd(z) # standardise

# CARD predictors (standadise)
df_transf <- dd1 %>% select(Setting,Antibiotic,ar12) %>% 
                spread(key=Antibiotic ,value=ar12)
z <- df_transf %>% select(-Setting) %>% as.matrix()
z_2a <- (z - mean(z)) / sd(z) # standardise
#
df_transf <- dd1 %>% select(Setting,Antibiotic,ar2) %>% 
                spread(key=Antibiotic ,value=ar2)
z <- df_transf %>% select(-Setting) %>% as.matrix()
z_2b <- (z - mean(z)) / sd(z) # standardise
#
df_transf <- dd1 %>% select(Setting,Antibiotic,ar1o2) %>% 
                spread(key=Antibiotic ,value=ar1o2)
z <- df_transf %>% select(-Setting) %>% as.matrix()
z_2c <- (z - mean(z)) / sd(z) # standardise

# impute counts where counts was 0 based on imputing from two other settings
count_s_for_pred <- count_s
for (i in 1:ncol(count_s)) {
                icol <-  count_s[,i]
                if (sum(is.na(icol)) > 0) {
                                icol_nonna <- icol[!is.na(icol)] 
                                count_s_for_pred[which(is.na(icol)),i] <- round(mean(icol_nonna))
                }
}


# run the Stan models -----------------------------------------------------


## model # c("ar12","ar2","ar1o2","b4div","entdiv","entbacdiv") = 2a 2b 2c , 1a 1b 1c
df_stan <- list(
                N=nrow(z),
                n=ncol(z),
                z1=z_1a, # e4 (b4div)
                z2=z_2a, # ar12 (coeff restricter to positive)
                count_s =count_s_hna,
                count_sr=count_sr_hna,
                data_exists=data_exists, 
                count_s_for_pred=count_s_for_pred
)
m00 <- stan(file=str_c(location_of_readme,"/stan/m_BAR_pred_sett_04_vary_d_miss.stan"),
            data = df_stan , 
            chains=5, iter=11000, warmup=1000, thin=10,cores = 1)
pop_1a2_e4 <- m00
loo_pop_1a2_e4 <- extract_log_lik(m00) %>% loo() 

## model # c("ar12","ar2","ar1o2","b4div","entdiv","entbacdiv") = 2a 2b 2c , 1a 1b 1c
df_stan <- list(
                N=nrow(z),
                n=ncol(z),
                z1=z_1b, # entdiv (bacterales)
                z2=z_2a, # ar12 (coeff restricter to positive)
                count_s =count_s_hna,
                count_sr=count_sr_hna,
                data_exists=data_exists, 
                count_s_for_pred=count_s_for_pred
)
m00 <- stan(file=str_c(location_of_readme,"/stan/m_BAR_pred_sett_04_vary_d_miss.stan"),
            data = df_stan , 
            chains=5, iter=11000, warmup=1000, thin=10)
pop_1a2_entles <- m00
loo_pop_1a2_entles <- extract_log_lik(m00) %>% loo() 

## model # c("ar12","ar2","ar1o2","b4div","entdiv","entbacdiv") = 2a 2b 2c , 1a 1b 1c
df_stan <- list(
                N=nrow(z),
                n=ncol(z),
                z1=z_1c, # entbacdiv (Enterobacteriaceae)
                z2=z_2a, # ar12 (coeff restricter to positive)
                count_s =count_s_hna,
                count_sr=count_sr_hna,
                data_exists=data_exists, 
                count_s_for_pred=count_s_for_pred
)
m00 <- stan(file=str_c(location_of_readme,"/stan/m_BAR_pred_sett_04_vary_d_miss.stan"),
            data = df_stan , 
            chains=5, iter=11000, warmup=1000, thin=10)
pop_1a2_entceae <- m00
loo_pop_1a2_entceae <- extract_log_lik(m00) %>% loo()

# without tax
m00 <- stan(file=str_c(location_of_readme,"/stan/m_BAR_pred_sett_04_vary_d_miss_notax.stan"),
            data = df_stan , 
            chains=5, iter=11000, warmup=1000, thin=10)
pop_1a2_notax<- m00
loo_pop_1a1_notax <- extract_log_lik(m00) %>% loo()

## model # c("ar12","ar2","ar1o2","b4div","entdiv","entbacdiv") = 2a 2b 2c , 1a 1b 1c
df_stan <- list(
                N=nrow(z),
                n=ncol(z),
                z1=z_1a, # b4div ()
                z2=z_2b, # ar2 (coeff restricter to positive)
                count_s =count_s_hna,
                count_sr=count_sr_hna,
                data_exists=data_exists, 
                count_s_for_pred=count_s_for_pred
)
m00 <- stan(file=str_c(location_of_readme,"/stan/m_BAR_pred_sett_04_vary_d_miss.stan"),
            data = df_stan , 
            chains=5, iter=11000, warmup=1000, thin=10)
pop_2_e4 <- m00
loo_pop_2_e4 <- extract_log_lik(m00) %>% loo()

## model # c("ar12","ar2","ar1o2","b4div","entdiv","entbacdiv") = 2a 2b 2c , 1a 1b 1c
df_stan <- list(
                N=nrow(z),
                n=ncol(z),
                z1=z_1b, # entdiv (bacterales)
                z2=z_2b, # ar12 (coeff restricter to positive)
                count_s =count_s_hna,
                count_sr=count_sr_hna,
                data_exists=data_exists, 
                count_s_for_pred=count_s_for_pred
)
m00 <- stan(file=str_c(location_of_readme,"/stan/m_BAR_pred_sett_04_vary_d_miss.stan"),
            data = df_stan , 
            chains=5, iter=11000, warmup=1000, thin=10)
pop_2_entles <- m00
loo_pop_2_entles <- extract_log_lik(m00) %>% loo()

## model # c("ar12","ar2","ar1o2","b4div","entdiv","entbacdiv") = 2a 2b 2c , 1a 1b 1c
df_stan <- list(
                N=nrow(z),
                n=ncol(z),
                z1=z_1c, # entbacdiv (Enterobacteriaceae)
                z2=z_2b, # ar12 (coeff restricter to positive)
                count_s =count_s_hna,
                count_sr=count_sr_hna,
                data_exists=data_exists, 
                count_s_for_pred=count_s_for_pred
)
m00 <- stan(file=str_c(location_of_readme,"/stan/m_BAR_pred_sett_04_vary_d_miss.stan"),
            data = df_stan , 
            chains=5, iter=11000, warmup=1000, thin=10)
pop_2_entceae<- m00
loo_pop_2_entceae <- extract_log_lik(m00) %>% loo()

# without tax
m00 <- stan(file=str_c(location_of_readme,"/stan/m_BAR_pred_sett_04_vary_d_miss_notax.stan"),
            data = df_stan , 
            chains=5, iter=11000, warmup=1000, thin=10)
pop_2_notax<- m00
loo_pop_2_notax <- extract_log_lik(m00) %>% loo()

# model without any covariates
m05 <- stan(file=str_c(location_of_readme,"/stan/m_BAR_pred_sett_05_miss.stan"),
            data = df_stan , 
            chains=5, iter=11000, warmup=1000, thin=10)
loo_m05 <- extract_log_lik(m05) %>% loo() 

# Model comparison --------------------------------------------------------


wts2 <- loo_model_weights(
                list(loo_m05,
                     loo_pop_2_entceae,
                     loo_pop_2_entles,
                     loo_pop_2_e4,
                     loo_pop_1a2_entceae,
                     loo_pop_1a2_entles,
                     loo_pop_1a2_e4,
                     loo_pop_2_notax,
                     loo_pop_1a1_notax),
                method = "pseudobma",
                optim_control = list(reltol=1e-10)
)

df_comp <- compare( loo_m05,
                    loo_pop_2_entceae,
                    loo_pop_2_entles,
                    loo_pop_2_e4,
                    loo_pop_1a2_entceae,
                    loo_pop_1a2_entles,
                    loo_pop_1a2_e4,
                    loo_pop_2_notax,
                    loo_pop_1a1_notax )
df_comp <- df_comp %>% as_tibble(rownames = "model") %>% 
                select(model,elpd_loo,elpd_diff) %>% 
                left_join( tibble(
                                model=c("loo_pop_1a2_e4",
                                        "loo_pop_1a2_entceae",
                                        "loo_pop_2_e4",
                                        "loo_pop_2_entceae",
                                        "loo_pop_1a2_entles",
                                        "loo_pop_2_entles",
                                        "loo_pop_1a1_notax",
                                        "loo_pop_2_notax",
                                        "loo_m05"), 
                                model_nname=c("ALL_e4",
                                              "ALL_entceae",
                                              "DEF_e4",
                                              "DEF_entceae",
                                              "ALL_entles",
                                              "DEF_entles",
                                              "ALL_notax",
                                              "DEF_notax",
                                              "baseline")
                ), by="model") %>% 
                mutate( model=model_nname ) %>% select(-model_nname) %>% 
                rename( Model=model,
                        loo_prediction=elpd_loo,
                        loo_diff_to_bestmodel=elpd_diff)

wts2 %>% enframe() %>% 
                mutate(name=c("baseline",
                              "DEF_entceae",
                              "DEF_entles",
                              "DEF_e4",
                              "ALL_entceae",
                              "ALL_entles",
                              "ALL_e4",
                              "DEF_notax",
                              "ALL_notax")) %>% 
                mutate(value=round(value,3)) %>% 
                rename(BMA_weight=value,
                       Model=name) %>% 
                left_join( df_comp, by="Model" ) %>% arrange(desc(BMA_weight)) %>% 
                mutate(loo_prediction=round(loo_prediction,2),
                       loo_diff_to_bestmodel=round(loo_diff_to_bestmodel,2)) 



# Plot model predictions with data ----------------------------------------

mrel <- 1
n_sett <- 3
n_genes <- n_atbs
# model with metagenomic information
fit_post <- rstan::extract(pop_1a2_e4)
post_pred <- fit_post$pred_count_sr
pred_lower <- matrix(NA,nrow=n_sett,ncol=n_genes)
pred_upper <- matrix(NA,nrow=n_sett,ncol=n_genes)
pred_mean <- matrix(NA,nrow=n_sett,ncol=n_genes)
for (i in 1:n_sett) {
                for (g in 1:n_genes){
                                post_PI <- bayestestR::hdi(post_pred[,i,g],ci=0.95 ) %>% unlist() 
                                pred_lower[i,g] <- post_PI[2]
                                pred_upper[i,g] <- post_PI[3]
                                pred_mean[i,g] <- mean(post_pred[,i,g])
                }
}

# wrrangle real data
df1 <- count_sr %>% as_tibble() %>% mutate(Setting=setting_v) %>% gather(key=Antibiotic,value=real,-Setting)
tot <- count_s_for_pred %>% as_tibble() %>% mutate(Setting=setting_v) %>% gather(key=Antibiotic,value=tot,-Setting) %>% .$tot
tot_withNA <- count_s %>% as_tibble() %>% mutate(Setting=setting_v) %>% gather(key=Antibiotic,value=tot,-Setting) %>% .$tot

# wrangle predictions
pred_mean_v <- pred_mean %>% as_tibble() %>% mutate(Setting=setting_v) %>% gather(key=Antibiotic,value=pred_mean,-Setting) %>% .$pred_mean
pred_lower_v <- pred_lower %>% as_tibble() %>% mutate(Setting=setting_v) %>% gather(key=Antibiotic,value=pred_lower,-Setting) %>% .$pred_lower
pred_upper_v <- pred_upper %>% as_tibble() %>% mutate(Setting=setting_v) %>% gather(key=Antibiotic,value=pred_upper,-Setting) %>% .$pred_upper

df2 <- df1 %>% mutate( pred_lower=pred_lower_v,
                       pred_upper=pred_upper_v,
                       pred_mean=pred_mean_v,
                       tot=tot,
                       tot_withNA=tot_withNA,
                       close_brack=")",
                       open_brack=" (") %>% 
                mutate(Setting=replace(Setting,Setting=="CAPOP","Cambodia"),
                       Setting=replace(Setting,Setting=="KEPOP","Kenya"),
                       Setting=replace(Setting,Setting=="UKPOP","UK"))

sett_df <- tibble( old_sett=c("CAMBODIA_POP_POOL",
                              "KENYA_POP_POOL",
                              "UK_POP_POOL"),
                   new_sett=c("Cambodia",
                              "Kenya",
                              "UK"))

(p1 <- df2 %>% arrange(Setting) %>% unite(col = Sett_Anti,Setting,Antibiotic,sep="-",remove=F) %>% 
                                unite(col = Antibiotic_count,Antibiotic,open_brack,tot_withNA,close_brack,sep="",remove=F) %>% 
                                # for color of bars
                                mutate( col_bars=ifelse(is.na(tot_withNA), "no_clin_cases" , Setting) ) %>% 
                                ggplot( aes(y=Antibiotic_count) ) + 
                                geom_point( aes(x=pred_mean, col=col_bars ),alpha=0.7,size=4,pch="|" ,show.legend = F) +
                                geom_segment( aes(yend=Antibiotic_count,x=pred_lower,xend=pred_upper,
                                                  col=col_bars),alpha=0.50,size=2.3,show.legend = F ) + # 1model
                                scale_color_manual(values=c('CAMBODIA_POP_POOL'='#efc750',
                                                            'KENYA_POP_POOL'='#a6cdd9',
                                                            'UK_POP_POOL'='#b7b079',
                                                            'no_clin_cases'='darkgrey') ) +
                                geom_point( aes(x=real), size=2) + 
                                facet_wrap(~Setting,ncol=1,scales="free") +
                                labs(x="Resistance Proportion",y="",title="")) 


