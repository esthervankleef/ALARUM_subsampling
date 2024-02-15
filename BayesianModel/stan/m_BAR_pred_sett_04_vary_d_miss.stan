 data {
   // give numbers and sizes
   int<lower=1> N; // number of Settings
   int<lower=1> n; // number of antibiotics
   // data 
   real z1[N,n]; // the score for each setting and antibiotic
   real z2[N,n]; // the score for each setting and antibiotic
   int<lower=0> count_sr[N,n]; // the count of sick-resistant patients in setting and antibiotic
   int<lower=0> count_s[N,n]; // the count of sick patients in setting and antibiotic
   int<lower=0, upper=1> data_exists[N,n];
   int<lower=0> count_s_for_pred[N,n]; // for prediction we cannot have missing number of sick
 }
 
 parameters {
   // Define parameters to estimate
   //real b_sett[N]; // coefficient for each antibiotic
   real b_z1[n];//
   real<lower=0> b_z2[n];//
   real b_res[n]; // mean of the varying coefficients
   real b_z_mu1; // mean for the varying coefficients
   real<lower=0> b_z_mu2; // mean for the varying coefficients

   real<lower=0> b_z_sigma1; // sigma of the varying coefficients
   real<lower=0> b_z_sigma2; // sigma of the varying coefficients

   //real<lower=0> b_sett_sigma; // sigma of the varying coefficients
 }
 
 transformed parameters  {
   // make parameter transformations
   real logit_p[N,n];
   // i: setting g: gene
   for (i in 1:N){ 
                   for (g in 1:n) {
                                 logit_p[i,g] = b_res[g] + b_z1[g]*z1[i,g] + b_z2[g]*z2[i,g]; // learns g-profile from all settings  + additional info about low and high
                   }
   }
 }
 
 model {
   // Priors (no need to specify if non-informative)
   b_z_mu1 ~ normal(0,1); // don't know where it will be
   b_z_mu2 ~ normal(0,1); // don't know where it will be
   b_z_sigma1 ~ normal(0,1);
   b_z_sigma2 ~ normal(0,1);
   // b_sett_sigma ~ cauchy(0,1); // restrict this because only 3 settings
   // through antibiotics
   // b_z ~ normal( 0,1 ); // same coefficient
   b_z1 ~ normal(b_z_mu1 , b_z_sigma1); // varying effect 
   b_z2 ~ normal(b_z_mu2 , b_z_sigma2); // varying effect 
   for (g in 1:n) {
                   b_res[g] ~ normal(0,1); // independent priors, how big of a difference do we allow?
                   } 

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
     real my_inv_logit ;
     for (i in 1:N){
                   for (g in 1:n) {
                                 my_inv_logit = 1 / (1 + exp(-logit_p[i,g])) ; // manually do the inverse 
                                 pred_count_sr[i,g] = binomial_rng( count_s_for_pred[i,g] , my_inv_logit ) ;
                                 log_lik[i,g] = binomial_logit_lpmf( count_sr[i,g] | count_s[i,g] , logit_p[i,g] );
                   }
   }
 }
