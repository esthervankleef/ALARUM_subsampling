 data {
   // give numbers and sizes
   int<lower=1> N; // number of Settings
   int<lower=1> n; // number of antibiotics
   // data 
   //real z[N,n]; // the score for each setting and antibiotic
   int<lower=0> count_sr[N,n]; // the count of sick-resistant patients in setting and antibiotic
   int<lower=0> count_s[N,n]; // the count of sick patients in setting and antibiotic
   int<lower=0,upper=1> data_exists[N,n];
   int<lower=0> count_s_for_pred[N,n]; // for prediction we cannot have missing number of sick
 }
 
 parameters {
   // Define parameters to estimate
   //real b_sett[N]; // coefficient for each antibiotic
   // real b_z;//
   real b_res[n]; // mean of the varying coefficients
   //real b_z_mu; // mean for the varying coefficients
   //real<lower=0> b_z_sigma; // sigma of the varying coefficients
   //real<lower=0> b_sett_sigma; // sigma of the varying coefficients
 }
 
 transformed parameters  {
   // make parameter transformations
   real logit_p[N,n];
   // i: setting g: gene
   for (i in 1:N){ 
                   for (g in 1:n) {
                                 logit_p[i,g] = b_res[g] ; // learns world g-profile only
                   }
   }
 }
 
 model {
   // Priors (no need to specify if non-informative)
   //b_z_mu ~ normal(0,5); // don't know where it will be
   //b_z_sigma ~ cauchy(0,1);
   // b_sett_sigma ~ cauchy(0,1); // restrict this because only 3 settings
   // through antibiotics
   // b_z ~ normal( 0,1 ); // same coefficient
   for (g in 1:n) b_res[g] ~ normal(0,1); // independent priors, how big of a difference do we allow? 

   // Likelihoods
   for (i in 1:N){
                   for (g in 1:n) {
                                 if (data_exists[i,g]==1) count_sr[i,g] ~ binomial_logit( count_s[i,g],logit_p[i,g] ) ;  
                   }
   }
}

 generated quantities {
   // define predicted vector
     int pred_count_sr[N,n];
     real log_lik[N,n];
     for (i in 1:N){
                   for (g in 1:n) {
                                 pred_count_sr[i,g] = binomial_rng( count_s_for_pred[i,g],inv_logit(logit_p[i,g]) ) ;
                                 log_lik[i,g] = binomial_logit_lpmf( count_sr[i,g] | count_s[i,g] , logit_p[i,g] ); 
                   }
   }
 }
