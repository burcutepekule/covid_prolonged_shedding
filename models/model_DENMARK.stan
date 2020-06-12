functions {
  real switch_lock(real t, real t_lock, real r_lock, real shift_lock, real m_lock) {
    return(r_lock+(1-r_lock)/(1+exp(m_lock*(t-t_lock-shift_lock))));
  }
  real switch_relax(real t, real t_relax, real r_lock, real r_end, real m_relax) {
    return(r_lock+1./(1/(r_end-r_lock)+exp(-m_relax*(t-t_relax-(6*(1/m_relax)*(1-m_relax)+1/m_relax)))));
  }
  real[] SEIR(real t,
  real[] y,
  real[] theta,
  real[] x_r,
  int[]  x_i 
  ) {
    real tswitch = x_r[1]; // lockdown time
    real trelax  = x_r[2];
    real r_end   = x_r[3];
    real m_relax = x_r[4];
    real gamma_c = x_r[5];
    int  r_c     = x_i[1]; //the only fixed integer is the reduction in inf due to being chronic
    real p_tswitch;
    real p_relax;
    real dydt[9];
    real init[2]; // inital intercept for S and E
    
    real pi; // number of cases at t0
    real R0; //
    real tau;
    real gamma_s;
    real gamma_H;
    real gamma_ICU;
    real eps_H;
    real eps_H2ICU;
    real eps_x_ICU;
    real r_d_s; // detection rate of symp
    real r_d_c; // detection rate of chronics
    real r_lock; // reduction in transmission rate after lockdown
    real shift_lock;
    real m_lock;
    real coeffR;
    
    // Free parameters
    pi        = theta[1];
    R0        = theta[2];
    tau       = theta[3];
    gamma_s   = theta[4];
    gamma_H   = theta[5];
    gamma_ICU = theta[6];
    eps_H     = theta[7];
    eps_H2ICU = theta[8];
    eps_x_ICU = theta[9];
    r_d_s     = theta[10];
    r_d_c     = theta[11];
    r_lock    = theta[12];
    shift_lock= theta[13];
    m_lock    = theta[14];
    
    p_tswitch = switch_lock(t,tswitch,r_lock,shift_lock,m_lock);
    
    if(m_relax>0){ //m_relax=0 no relaxation
    p_relax   = switch_relax(t,trelax,r_lock,r_end,m_relax);
    if(p_tswitch>p_relax){
      coeffR = p_tswitch;
    }
    else{
      coeffR = p_relax;
    }
    }else{
      coeffR = p_tswitch;
    }
    
    // Initial conditions
    init[1] = 1-pi; // -> sensitives, the actual initial contidition (this is for speed, check below)
    init[2] = pi; // -> exposed

    dydt[1] = - coeffR*(y[1]+init[1])*(R0*gamma_s*y[3]+R0*gamma_c*(1-0.01*r_c)*y[4]); // -lockdown*S*(beta_s*Is+beta_c*Ic)
    // E
    dydt[2] = + coeffR*(y[1]+init[1])*(R0*gamma_s*y[3]+R0*gamma_c*(1-0.01*r_c)*y[4]) - tau*(y[2]+init[2]); // +lockdown*S*(beta_s*Is+beta_c*Ic) - tau*E
    // Is
    dydt[3] = + tau*(y[2]+init[2]) - gamma_s*y[3];
    // Ic
    dydt[4] = + (1-eps_H)*gamma_s*y[3] - gamma_c*y[4];
    // H
    dydt[5] = +eps_H*gamma_s*y[3] - gamma_H*y[5];
    // ICU
    dydt[6] = +gamma_H*eps_H2ICU*y[5] - gamma_ICU*y[6];
    // R
    dydt[7] = +gamma_H*(1-eps_H2ICU)*y[5]+gamma_ICU*(1-eps_x_ICU)*y[6]+gamma_c*y[4]; 
    // X      
    dydt[8] = +gamma_ICU*eps_x_ICU*y[6];
    // C_p_d
    dydt[9]  = +r_d_s*gamma_s*y[3]+(1-(r_d_s-eps_H))*r_d_c*gamma_c*y[4];
    return(dydt);
  }
}

data {
  // data
  int pop_t; // total population
  real tswitch; 
  real trelax;
  real r_end;
  real m_relax;
  real gamma_c;
  int D; // number of days with reported incidence
  int k_daily_cases[D];
  int k_daily_deaths[D];
  
  // priors
  real p_pi[2];
  real p_R0[2];
  real p_tau;
  real p_gamma_s;
  real p_gamma_H;
  real p_gamma_ICU;
  real p_eps_H[2];
  real p_eps_H2ICU[2];
  real p_eps_x_ICU[2];
  real p_r_d_s[2];
  real p_r_d_c[2];
  real p_r_lock[2];
  real p_phi;
  
  // Simulation
  real t0; //starting time
  int t_data; //time of first data
  int S; // total simulation time
  int E; // extra simulation time (S=D+E)
  real ts[D]; // time bins
  real ts_pred[S];
  
  // controls
  int r_c;
}

transformed data {
  real x_r[5]; 
  int x_i[1]; // this is for r_c
  real init[9] = rep_array(1e-9,9); // initial values -> this is for speed, keep it like that
  x_r[1] = tswitch;
  x_r[2] = trelax;
  x_r[3] = r_end;
  x_r[4] = m_relax;
  x_r[5] = gamma_c;
  x_i[1] = r_c;
}

parameters{
  real<lower=0, upper=1> pi; // number of cases at t0
  real<lower=0> R0;
  real<lower=0, upper=1> tau; 
  real<lower=0, upper=1> gamma_s;  
  real<lower=0, upper=1> gamma_H; 
  real<lower=0, upper=1> gamma_ICU;
  real<lower=0, upper=1> eps_H;
  real<lower=0, upper=1> eps_H2ICU;
  real<lower=0, upper=1> eps_x_ICU;
  real<lower=0, upper=1> r_d_s;
  real<lower=0, upper=1> r_d_c;
  real<lower=0, upper=1> r_lock;
  real<lower=0, upper=1> m_lock_raw; // slope of quarantine implementation
  real<lower=0> shift_lock; // shift of quarantine implementation
  real<lower=0> phi[2]; // dispersion parameters
}
transformed parameters {
  real m_lock = m_lock_raw+0.5;
  real theta[14]; // vector of parameters
  real y[D,9]; // raw ODE output
  vector[D] output_cumC;
  vector[D] output_cumX;
  // outcomes
  vector[D] output_k_daily_deaths; 
  vector[D] output_k_daily_cases;
  
  theta = {pi,R0,tau,gamma_s,gamma_H,gamma_ICU,eps_H,eps_H2ICU,eps_x_ICU,r_d_s,r_d_c,r_lock,shift_lock,m_lock};
  // run ODE solver
  y = integrate_ode_bdf(
    SEIR, // ODE function
    init, // initial states
    t0, // t0
    ts, // evaluation dates (ts)
    theta, // parameters
    x_r, // real data
    x_i, // integer data
    1.0E-10, 1.0E-10, 1.0E3); // tolerances and mam_lockmum steps
    // extract and format ODE results (1.0E-9 correction to avoid negative values due to unprecise estimates of zeros as tolerance is 1.0E-10)
    for(i in 1:D) {
      output_cumC[i]= (y[i,9]+1.0E-9)*pop_t;
      output_cumX[i]= y[i,8]*pop_t;
      // fitted outputs
      output_k_daily_cases[i] = i==1 ? output_cumC[i] : 1.0E-9*pop_t + output_cumC[i] - output_cumC[i-1]; 
      output_k_daily_deaths[i] = i==1 ? output_cumX[i] : 1.0E-9*pop_t + output_cumX[i] - output_cumX[i-1]; 
    }
}

model {
  // priors
  pi ~ beta(p_pi[1],p_pi[2]);
  R0 ~ normal(p_R0[1],p_R0[2]);
  tau~ exponential(p_tau);
  gamma_s ~ exponential(p_gamma_s);
  gamma_H ~ exponential(p_gamma_H);
  gamma_ICU ~ exponential(p_gamma_ICU);
  eps_H ~ normal(p_eps_H[1],p_eps_H[2]);
  eps_H2ICU ~ normal(p_eps_H2ICU[1],p_eps_H2ICU[2]);
  eps_x_ICU ~ normal(p_eps_x_ICU[1],p_eps_x_ICU[2]);
  r_d_s ~ normal(p_r_d_s[1],p_r_d_s[2]);
  r_d_c ~ normal(p_r_d_c[1],p_r_d_c[2]);
  r_lock ~ beta(p_r_lock[1],p_r_lock[2]);
  phi ~ exponential(p_phi);
  m_lock_raw ~ beta(1,1); 
  shift_lock ~ exponential(1/5.0);
  // likelihood
  for(i in 1:D) {
    target += neg_binomial_2_lpmf( k_daily_cases[i] | output_k_daily_cases[i], phi[1]);
    target += neg_binomial_2_lpmf( k_daily_deaths[i] | output_k_daily_deaths[i], phi[2]);
  }
}

generated quantities{
  
  real y_pred[S,9]; // raw ODE output
  
  vector[S] comp_S;
  vector[S] comp_E;
  vector[S] comp_Is;
  vector[S] comp_Ic;
  vector[S] comp_H;
  vector[S] comp_ICU;
  vector[S] comp_R;
  vector[S] comp_X;
  vector[S] comp_C;
  vector[S] comp_diffX;
  vector[S] comp_diffC;
  
  int fitted_k_daily_cases[D];
  int fitted_k_daily_deaths[D];
  
  int predicted_k_daily_cases[E];
  int predicted_k_daily_deaths[E];
  
  real predicted_daily_cases[S];
  real predicted_cum_cases[S];
  real predicted_current_hospit[S];
  real predicted_current_icu[S];
  real predicted_daily_deaths[S];
  real predicted_cum_deaths[S];
  
  y_pred = integrate_ode_bdf(
    SEIR, // ODE function
    init, // initial states
    t0, // t0
    ts_pred, // evaluation dates (ts)
    theta, // parameters
    x_r, // real data
    x_i, // integer data
    1.0E-10, 1.0E-10, 1.0E3);
    
    for(i in 1:S) {
      
      comp_S[i] = (y_pred[i,1] + 1 - pi)*pop_t;
      comp_E[i] = (y_pred[i,2] + pi)*pop_t;
      comp_Is[i] = y_pred[i,3]*pop_t;
      comp_Ic[i] = y_pred[i,4]*pop_t;
      comp_H[i]  = y_pred[i,5]*pop_t;
      comp_ICU[i]= y_pred[i,6]*pop_t;
      comp_R[i]  = y_pred[i,7]*pop_t;
      comp_X[i]  = y_pred[i,8]*pop_t;
      comp_C[i]  = (y_pred[i,9]+1.0E-9)*pop_t;
      // lagged differences of cumulative compartments (daily deaths and daily cases)
      comp_diffX[i] = i==1 ? comp_X[i] : 1.0E-9*pop_t + comp_X[i] - comp_X[i-1];
      comp_diffC[i] = i==1 ? comp_C[i] : 1.0E-9*pop_t + comp_C[i] - comp_C[i-1];
    }
    
    for(i in 1:D) {
      fitted_k_daily_cases[i]  = neg_binomial_2_rng( output_k_daily_cases[i], phi[1]);
      fitted_k_daily_deaths[i] = neg_binomial_2_rng( output_k_daily_deaths[i], phi[2]);
    }
    for(i in 1:E) {
      predicted_k_daily_cases[i]   = neg_binomial_2_rng( comp_diffC[D+i], phi[1]);
      predicted_k_daily_deaths[i] = neg_binomial_2_rng( comp_diffX[D+i], phi[2]);
    }
    for(i in 1:S) {
      predicted_current_hospit[i] = comp_H[i];
      predicted_current_icu[i] = comp_ICU[i];
      predicted_daily_deaths[i] = comp_diffX[i];
      predicted_cum_deaths[i] = comp_X[i];
      predicted_daily_cases[i] = comp_diffC[i];
      predicted_cum_cases[i] = comp_C[i];
    }
}
