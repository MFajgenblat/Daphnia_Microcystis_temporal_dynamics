data {
  
  int<lower=0> N;
  int<lower=0> N_clones;
  int<lower=0> N_strains;
  int<lower=0> N_clonestraincombinations;
  array[N] int Survived_individuals;
  array[N] int Total_individuals;
  vector[N] Daph_season;
  vector[N] MC_season;
  vector[N] Synchronicity;
  array[N] int Daph_year;
  array[N] int Daph_clone;
  array[N] int MC_strain;
  array[N] int Clonestraincombination;
  int GxG;
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
  target += binomial_logit_lpmf(Survived_individuals | Total_individuals, mu);

}

generated quantities {
  
  // Posterior predictive samples
  array[N] int Survived_individuals_rep;
  // Estimated survival per treatment
  array[2,4] real average_survival;
  
  Survived_individuals_rep = binomial_rng(Total_individuals, inv_logit(mu));
  
  for (i in 1:2) {
    // Early Daphnia, early Microcystis
    average_survival[i,1] = inv_logit(beta_0[i] + beta_1[i]*0 + beta_2[i]*0 + beta_3[i]*1);
    // Early Daphnia, late Microcystis
    average_survival[i,2] = inv_logit(beta_0[i] + beta_1[i]*0 + beta_2[i]*1 + beta_3[i]*0);
    // Late Daphnia, early Microcystis
    average_survival[i,3] = inv_logit(beta_0[i] + beta_1[i]*1 + beta_2[i]*0 + beta_3[i]*0);
    // Late Daphnia, late Microcystis
    average_survival[i,4] = inv_logit(beta_0[i] + beta_1[i]*1 + beta_2[i]*1 + beta_3[i]*1);
  }


}
