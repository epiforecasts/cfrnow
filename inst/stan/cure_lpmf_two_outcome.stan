/* Two-outcome mixture-cure log-likelihood (death + timed recovery).

   Holes as in cure_lpmf_death.stan, plus the recovery delay (its own family):
     <<recovery_pars>>     recovery-delay parameter declarations (r-prefixed)
     <<recovery_id>>       primarycensored recovery distribution id
     <<recovery_reparam>>  recovery parameters in primarycensored's native order
     <<recovery_upper>>    upper truncation of the recovery delay
     <<loss_pars>>         loss-to-follow-up parameter declarations, empty when
                           no loss is estimated
     <<death_kept>>/<<recovery_kept>>
                           "log1m(<loss>) + " for a recorded outcome
     <<unresolved_terms>>  the mixture an unresolved case sums over

   outcome: 1 = death, 2 = timed recovery, 3 = resolved, else = censored. */
real cfrnow_<<family>>_lpmf(data int y, <<death_pars>>, real prob,
                            <<recovery_pars>>, <<loss_pars>>data real outcome,
                            data real pwindow, data real swindow,
                            array[] real primary_params) {
  if (outcome == 1) {
    return <<death_kept>>log(prob) + primarycensored_lpmf(
        y | <<death_id>>, {<<death_reparam>>}, pwindow, y + swindow, 0.0,
        <<death_upper>>, <<primary_id>>, primary_params);
  } else if (outcome == 2) {
    return <<recovery_kept>>log1m(prob) + primarycensored_lpmf(
        y | <<recovery_id>>, {<<recovery_reparam>>}, pwindow, y + swindow, 0.0,
        <<recovery_upper>>, <<primary_id>>, primary_params);
  } else if (outcome == 3) {
    return <<recovery_kept>>log1m(prob);
  } else {
    // On the log scale: either cdf can round to (just above) 1 for a long
    // follow-up, where log1m() would reject the draw. The floor keeps such a
    // case very unlikely instead of impossible, so the sampler can move.
    real log_surv_d = log1m_exp(fmin(primarycensored_lcdf(
        y | <<death_id>>, {<<death_reparam>>}, pwindow, 0.0,
        <<death_upper>>, <<primary_id>>, primary_params), -1e-12));
    real log_surv_r = log1m_exp(fmin(primarycensored_lcdf(
        y | <<recovery_id>>, {<<recovery_reparam>>}, pwindow, 0.0,
        <<recovery_upper>>, <<primary_id>>, primary_params), -1e-12));
    return log_sum_exp({<<unresolved_terms>>});
  }
}
