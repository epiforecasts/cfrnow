/* Two-outcome mixture-cure log-likelihood (death + timed recovery).

   Holes as in cure_lpmf_death.stan, plus the recovery delay (its own family):
     <<recovery_pars>>     recovery-delay parameter declarations (r-prefixed)
     <<recovery_id>>       primarycensored recovery distribution id
     <<recovery_reparam>>  recovery parameters in primarycensored's native order
     <<recovery_upper>>    upper truncation of the recovery delay

   outcome: 1 = death, 2 = timed recovery, 3 = resolved, else = censored. */
real cfrnow_<<family>>_lpmf(data int y, <<death_pars>>, real prob,
                            <<recovery_pars>>, data real outcome,
                            data real pwindow, data real swindow,
                            array[] real primary_params) {
  if (outcome == 1) {
    return log(prob) + primarycensored_lpmf(
        y | <<death_id>>, {<<death_reparam>>}, pwindow, y + swindow, 0.0,
        <<death_upper>>, <<primary_id>>, primary_params);
  } else if (outcome == 2) {
    return log1m(prob) + primarycensored_lpmf(
        y | <<recovery_id>>, {<<recovery_reparam>>}, pwindow, y + swindow, 0.0,
        <<recovery_upper>>, <<primary_id>>, primary_params);
  } else if (outcome == 3) {
    return log1m(prob);
  } else {
    real fbar_d = primarycensored_cdf(
        y | <<death_id>>, {<<death_reparam>>}, pwindow, 0.0,
        <<death_upper>>, <<primary_id>>, primary_params);
    real fbar_r = primarycensored_cdf(
        y | <<recovery_id>>, {<<recovery_reparam>>}, pwindow, 0.0,
        <<recovery_upper>>, <<primary_id>>, primary_params);
    return log_sum_exp(log(prob) + log1m(fbar_d), log1m(prob) + log1m(fbar_r));
  }
}
