data {
  int<lower=1> N;                                   // Total time points
  int<lower=1> M;                                   // Number of dimensions
  matrix[N, M] Y;                                   // Event matrix
  matrix[N, M] A;                                   // Alarm event matrix
  int day_or_night[N];                              // Indicator for day (0) or night (1)
  real seasonal_hist[N];                            // Seasonal historical factor
  int total_number_events;                          // total number of events in any dimension (=length(total_events))
  int total_number_alarm_events;                    // total number of events alarm in any dimension (=length(total_events))
  int total_events[total_number_events];            // Times where at least one event occurs
  int total_alarm_events[total_number_alarm_events];// Times where at least one alarm occurs
  real alpha_hyp;                                   // Hyperparameter for alpha
  real K_hyp;                                       // Hyperparameter for K
  int mu_size_per;                                  // size of each mu period (e.g., 365 for yearly baseline)
  int rounded;                                      // rounded = ceil(N / mu_size_per)
  int mu_index_per_event[total_number_events];      // index of mu vector at each event time
  matrix[N, M] Y_day;                               // Day-time event matrix
  matrix[N, M] A_day;                               // Day-time alarm event matrix 
  matrix[N, M] Y_night;                             // Night-time event matrix    
  matrix[N, M] A_night;                             // Night-time alarm event matrix
}

parameters {
  real<lower=0, upper=1> mu[rounded, M];            // Baseline intensity parameter
  real<lower=1e-5, upper=0.999> K[M,M];          // Influence matrix (day)
  real<lower=1e-5, upper=0.999> alpha[M,M];      // Influence matrix (alarm)
  real<lower=0, upper=1> beta_d;                    // Decay parameter (day)
  real<lower=0, upper=1> beta_s;                    // Decay parameter (self-excitation)
  real<lower=0, upper=1> K_n;                       // Influence matrix (night)
  real<lower=0, upper=1> alpha_n;                   // Influence matrix (alarm, night)
  real<lower=0, upper=1> beta_d_alarm;                    // Decay parameter (day) after alarm
  real<lower=0, upper=1> beta_s_alarm;                    // Decay parameter (self-excitation) after alarm
}

transformed parameters {
  matrix[mu_size_per*rounded,M] mu_total;           // repeating mu over the required interval. We actually repeat it more than we need to for ease

  for (m in 1:M){
    for (t in 1:rounded){
      mu_total[((t-1)*mu_size_per+1):(t*mu_size_per),m]=rep_vector(mu[t,m],mu_size_per);
    }
  }
  
}

model {
  matrix[total_number_events, M] R_diff_day = rep_matrix(0.0, total_number_events, M);
  matrix[total_number_events, M] R_diff_night = rep_matrix(0.0, total_number_events, M);
  matrix[total_number_events, M] R_alarm_diff_day = rep_matrix(0.0, total_number_events, M);
  matrix[total_number_events, M] R_alarm_diff_night = rep_matrix(0.0, total_number_events, M);
  matrix[total_number_events, M] R_same_day = rep_matrix(0.0, total_number_events, M);
  matrix[total_number_events, M] R_same_night = rep_matrix(0.0, total_number_events, M);
  matrix[total_number_events, M] R_alarm_same_day = rep_matrix(0.0, total_number_events, M);
  matrix[total_number_events, M] R_alarm_same_night = rep_matrix(0.0, total_number_events, M);
  
  real sum1 = 0;
  real sum2 = 0;
  real sum3 = 0;
  real sum4 = 0;
  real sum5 = 0;
  real sum6 = 0;

  real mu_reg = 0;
  real K_reg = 0;
  real alpha_reg = 0;

  // Initialize R matrices
  for (l in 1:M) {
    for (j in 2:total_number_events) {
      R_diff_day[j, l] = (1 - beta_d)^(total_events[j] - total_events[j - 1]) * R_diff_day[j - 1, l]
                          + Y_day[total_events[j - 1], l] * beta_d * (1 - beta_d)^(total_events[j] - total_events[j - 1] - 1);
      R_diff_night[j, l] = (1 - beta_d)^(total_events[j] - total_events[j - 1]) * R_diff_night[j - 1, l]
                           + Y_night[total_events[j - 1], l] * beta_d * (1 - beta_d)^(total_events[j] - total_events[j - 1] - 1);
      R_alarm_diff_day[j, l] = (1 - beta_d_alarm)^(total_events[j] - total_events[j - 1]) * R_alarm_diff_day[j - 1, l]
                               + A_day[total_events[j - 1], l] * beta_d_alarm * (1 - beta_d_alarm)^(total_events[j] - total_events[j - 1] - 1);
      R_alarm_diff_night[j, l] = (1 - beta_d_alarm)^(total_events[j] - total_events[j - 1]) * R_alarm_diff_night[j - 1, l]
                                 + A_night[total_events[j - 1], l] * beta_d_alarm * (1 - beta_d_alarm)^(total_events[j] - total_events[j - 1] - 1);
    }
  }

  for (l in 1:M) {
    for (j in 2:total_number_events) {
      R_same_day[j, l] = (1 - beta_s)^(total_events[j] - total_events[j - 1]) * R_same_day[j - 1, l]
                          + Y_day[total_events[j - 1], l] * beta_s * (1 - beta_s)^(total_events[j] - total_events[j - 1] - 1);
      R_same_night[j, l] = (1 - beta_s)^(total_events[j] - total_events[j - 1]) * R_same_night[j - 1, l]
                           + Y_night[total_events[j - 1], l] * beta_s * (1 - beta_s)^(total_events[j] - total_events[j - 1] - 1);
      R_alarm_same_day[j, l] = (1 - beta_s_alarm)^(total_events[j] - total_events[j - 1]) * R_alarm_same_day[j - 1, l]
                               + A_day[total_events[j - 1], l] * beta_s_alarm * (1 - beta_s_alarm)^(total_events[j] - total_events[j - 1] - 1);
      R_alarm_same_night[j, l] = (1 - beta_s_alarm)^(total_events[j] - total_events[j - 1]) * R_alarm_same_night[j - 1, l]
                                 + A_night[total_events[j - 1], l] * beta_s_alarm * (1 - beta_s_alarm)^(total_events[j] - total_events[j - 1] - 1);
    }
  }

  // sum1
  vector[M] lam;
  for (m in 1:M) {
    for (j in 1:total_number_events) {
      
      lam[m] = mu[mu_index_per_event[j], m]*seasonal_hist[total_events[j]];

      for (l in 1:M) {
        if (l != m) {
          lam[m] += K[l, m] * R_diff_day[j, l] + K_n * K[l, m] * R_diff_night[j, l] 
                 + alpha[l, m] * R_alarm_diff_day[j, l] + alpha_n * alpha[l, m] * R_alarm_diff_night[j, l];
        }
      }
      lam[m] += K[m, m] * R_same_day[j, m] + K_n * K[m, m] * R_same_night[j, m]
             + alpha[m, m] * R_alarm_same_day[j, m] + alpha_n * alpha[m, m] * R_alarm_same_night[j, m];
      sum1 += Y[total_events[j], m] * log(lam[m]);
    }
  }

  // sum2: Day event influence terms
  for (l in 1:M) {
    for (m in 1:M) {
      if (l != m) {
        real inside_sum = 0;
        for (t in total_events) {
          inside_sum += Y_day[t, l] * (1 - (1 - beta_d)^(N - t));
        }
        sum2 += K[l, m] * inside_sum;
      }
    }
  }

  for (m in 1:M) {
    real inside_sum = 0;
    for (t in total_events) {
      inside_sum += Y_day[t, m] * (1 - (1 - beta_s)^(N - t));
    }
    sum2 += K[m, m] * inside_sum;
  }

  // sum3: Night event influence terms
  for (l in 1:M) {
    for (m in 1:M) {
      if (l != m) {
        real inside_sum = 0;
        for (t in total_events) {
          inside_sum += Y_night[t, l] * (1 - (1 - beta_d)^(N - t));
        }
        sum3 += K_n * K[l, m] * inside_sum;
      }
    }
  }

  for (m in 1:M) {
    real inside_sum = 0;
    for (t in total_events) {
      inside_sum += Y_night[t, m] * (1 - (1 - beta_s)^(N - t));
    }
    sum3 += K_n * K[m, m] * inside_sum;
  }

  // sum4: Alarm influence terms (day)
  for (l in 1:M) {
    for (m in 1:M) {
      if (l != m) {
        real inside_sum = 0;
        for (t in total_alarm_events) {
          inside_sum += A_day[t, l] * (1 - (1 - beta_d_alarm)^(N - t));
        }
        sum4 += alpha[l, m] * inside_sum;
      }
    }
  }

  for (m in 1:M) {
    real inside_sum = 0;
    for (t in total_alarm_events) {
      inside_sum += A_day[t, m] * (1 - (1 - beta_s_alarm)^(N - t));
    }
    sum4 += alpha[m, m] * inside_sum;
  }

  // sum5: Alarm influence terms (night)
  for (l in 1:M) {
    for (m in 1:M) {
      if (l != m) {
        real inside_sum = 0;
        for (t in total_alarm_events) {
          inside_sum += A_night[t, l] * (1 - (1 - beta_d_alarm)^(N - t));
        }
        sum5 += alpha_n * alpha[l, m] * inside_sum;
      }
    }
  }

  for (m in 1:M) {
    real inside_sum = 0;
    for (t in total_alarm_events) {
      inside_sum += A_night[t, m] * (1 - (1 - beta_s_alarm)^(N - t));
    }
    sum5 += alpha_n * alpha[m, m] * inside_sum;
  }



   for (m in 1:M){
    for (n in 1:N){
      sum6=sum6+mu_total[n,m]*seasonal_hist[n];
    }
  }
  
  

  // Regularization terms
  for (l in 1:M){
    for (m in 1:M){
      alpha_reg=alpha_reg+alpha[l,m];
      if (l != m){
        K_reg=K_reg+K[l,m];
      }
    }
  }

  // Full log-likelihood
  target += (sum1 - sum2 - sum3 - sum4 - sum5 - sum6)
             - (alpha_hyp * alpha_reg + K_hyp *K_reg);

}
