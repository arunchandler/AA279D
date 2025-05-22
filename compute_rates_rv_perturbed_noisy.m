function rv_dot_noisy = compute_rates_rv_perturbed_noisy(t, rv)
  global sigma
  rv_dot = compute_rates_rv_perturbed(t, rv);
  a_true = rv_dot(4:6);
  a_noisy = a_true + sigma*randn(3,1);
  rv_dot_noisy = [rv_dot(1:3); a_noisy];
end