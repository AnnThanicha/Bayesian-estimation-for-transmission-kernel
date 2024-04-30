data {
  int<lower=0> ID_day;
  int<lower=0> N;
  int<lower=0> K;
  vector[K] distance_kernel; 
  array[ID_day] int start; 
  array[ID_day] int stop;
  array[ID_day] int<lower=0> event; 
}

parameters {
  real<lower=0> k0;
  real<lower=0> r0;
  real<lower=0> alpha;
}

model {
  k0 ~ normal(0.005, 0.004) T[0, ]; 
  r0 ~ normal(0.19, 0.08) T[0, ];
  alpha ~ normal(1.56, 0.14) T[0, ];
  
  vector[K] kernel = ((1+((distance_kernel/r0)^alpha))^-1)*k0;
  vector[ID_day] prob;
  
  for (i in 1:ID_day) {
    
    if (start[i] == stop[i]){
      prob[i] = kernel[start[i]];
    }
    else{prob[i] = sum(kernel[start[i]:stop[i]]);
    }
   
  }

  event[1:ID_day] ~ bernoulli( 1 - exp(-1*(prob[1:ID_day])));
}
