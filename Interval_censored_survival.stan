data {
  
  int<lower=0> N;
  int<lower=0> N_clones;
  int<lower=0> N_strains;
  int<lower=0> N_clonestraincombinations;
  array[N] int<lower=1,upper=4> Interval;
  vector[N] Daph_season;
  vector[N] MC_season;
  vector[N] Synchronicity;
  array[N] int Daph_year;
  array[N] int Daph_clone;
  array[N] int MC_strain;
  array[N] int Clonestraincombination;
  int<lower=0,upper=1> GxG;
  int N_preds;
  vector[N_preds] pred_grid;
  array[2] real beta_prior;
  array[2] real re_sd_prior;
  
}

parameters {
  
  // Fixed effects
  array[2] real beta_0;
  array[2] real beta_1;
  array[2] real beta_2;
  array[2] real beta_3;
  
  // Random effects
  vector[N_clones] clone_std;
  vector[N_strains] strain_std;
  array[GxG]vector[N_clonestraincombinations] clonestraincombination_std;
  real<lower=0> clone_sd;
  real<lower=0> strain_sd;
  array[GxG]real<lower=0> clonestraincombination_sd;
  
}

transformed parameters {
  
  vector[N_clones] clone_effect = clone_std * clone_sd;
  vector[N_strains] strain_effect = strain_std * strain_sd;
  array[GxG] vector[N_clonestraincombinations] clonestraincombination_effect;
  vector[N] mu;
  
  // Without genotype x genotype interaction random effect
  if (GxG == 0) {
    mu = to_vector(beta_0[Daph_year]) + to_vector(beta_1[Daph_year]) .* Daph_season + to_vector(beta_2[Daph_year]) .* MC_season + to_vector(beta_3[Daph_year]) .* Synchronicity + clone_effect[Daph_clone] + strain_effect[MC_strain];
  }
  // With genotype x genotype interaction random effect
  if (GxG == 1) {
    clonestraincombination_effect[1] = clonestraincombination_std[1] * clonestraincombination_sd[1];
    mu = to_vector(beta_0[Daph_year]) + to_vector(beta_1[Daph_year]) .* Daph_season + to_vector(beta_2[Daph_year]) .* MC_season + to_vector(beta_3[Daph_year]) .* Synchronicity + clone_effect[Daph_clone] + strain_effect[MC_strain] + clonestraincombination_effect[1,Clonestraincombination];
  }

}

model {
  
  // Priors on fixed effects
  for (i in 1:2) {
    target += normal_lpdf(beta_0[i] | beta_prior[1], beta_prior[2]);
    target += normal_lpdf(beta_1[i] | beta_prior[1], beta_prior[2]);
    target += normal_lpdf(beta_2[i] | beta_prior[1], beta_prior[2]);
    target += normal_lpdf(beta_3[i] | beta_prior[1], beta_prior[2]);
  }
  
  // Priors on random effects
  target += std_normal_lpdf(clone_std);
  target += std_normal_lpdf(strain_std);
  if (GxG == 1) target += std_normal_lpdf(clonestraincombination_std[1]);
  target += normal_lpdf(clone_sd | re_sd_prior[1], re_sd_prior[2]);
  target += normal_lpdf(strain_sd | re_sd_prior[1], re_sd_prior[2]);
  if (GxG == 1) target += normal_lpdf(clonestraincombination_sd[1] | re_sd_prior[1], re_sd_prior[2]);
  
  // Likelihood
  for (i in 1:N) {
    if (Interval[i] == 1) {target += exponential_lcdf(0.5 | exp(-mu[i]));}
    if (Interval[i] == 2) {target += log_diff_exp(exponential_lcdf(0.9 | exp(-mu[i])), exponential_lcdf(0.5 | exp(-mu[i])));}
    if (Interval[i] == 3) {target += log_diff_exp(exponential_lcdf(1.2 | exp(-mu[i])), exponential_lcdf(0.9 | exp(-mu[i])));}
    if (Interval[i] == 4) {target += exponential_lccdf(1.2 | exp(-mu[i]));}
  }

}

generated quantities {
  
  // Posterior predictive samples
  array[N] real Survival_times;
  // Estimated survival per treatment through time
  array[3,2] vector[N_preds] S_averaged;
  
  Survival_times = exponential_rng(exp(-mu));
  
  // Early Daphnia
  S_averaged[1,1] = exp(-exp((beta_0[1] + beta_0[2])/2 - 0.5*0.5*(beta_1[1] + beta_1[2]))*pred_grid);
  // Late Daphnia
  S_averaged[1,2] = exp(-exp((beta_0[1] + beta_0[2])/2 + 0.5*0.5*(beta_1[1] + beta_1[2]))*pred_grid);
  // Early Microcystis
  S_averaged[2,1] = exp(-exp((beta_0[1] + beta_0[2])/2 - 0.5*0.5*(beta_2[1] + beta_2[2]))*pred_grid);
  // Late Microcystis
  S_averaged[2,2] = exp(-exp((beta_0[1] + beta_0[2])/2 + 0.5*0.5*(beta_2[1] + beta_2[2]))*pred_grid);
  // Contemporal
  S_averaged[3,1] = exp(-exp((beta_0[1] + beta_0[2])/2 - 0.5*0.5*(beta_3[1] + beta_3[2]))*pred_grid);
  // Allotemporal
  S_averaged[3,2] = exp(-exp((beta_0[1] + beta_0[2])/2 + 0.5*0.5*(beta_3[1] + beta_3[2]))*pred_grid);
  
}
