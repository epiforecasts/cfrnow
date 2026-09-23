/* Death-only mixture-cure log-likelihood for one case.

   Holes (<<name>>) are filled in by epidist_stancode():
     <<family>>         delay family name (lognormal / gamma)
     <<death_pars>>     delay parameter declarations, e.g. "real mu, real sigma"
     <<death_id>>       primarycensored delay distribution id
     <<death_reparam>>  delay parameters in primarycensored's native order
     <<death_upper>>    upper truncation of the delay (positive_infinity() when
                        the delay has no max)
     <<primary_id>>     primarycensored primary (uniform) distribution id
     <<loss_pars>>      "real loss, " when the fit allows loss to follow-up,
                        otherwise empty
     <<observed>>       "log1m(loss) + " for a case whose outcome was recorded
     <<lost_open>>/<<lost_close>>
                        wrap the unresolved case's term in the loss mixture

   outcome: 1 = observed death, 3 = resolved non-death, else = censored. */
real cfrnow_<<family>>_lpmf(data int y, <<death_pars>>, real prob,
                            <<loss_pars>>data real outcome, data real pwindow,
                            data real swindow, array[] real primary_params) {
  if (outcome == 1) {
    return <<observed>>log(prob) + primarycensored_lpmf(
        y | <<death_id>>, {<<death_reparam>>}, pwindow, y + swindow, 0.0,
        <<death_upper>>, <<primary_id>>, primary_params);
  } else if (outcome == 3) {
    return <<observed>>log1m(prob);
  } else {
    // On the log scale: the cdf can round to (just above) 1 for a long
    // follow-up, where log1m() would reject the draw. The floor keeps such a
    // case very unlikely instead of impossible, so the sampler can move.
    real log_fbar = primarycensored_lcdf(
        y | <<death_id>>, {<<death_reparam>>}, pwindow, 0.0,
        <<death_upper>>, <<primary_id>>, primary_params);
    return <<lost_open>>log1m_exp(fmin(log(prob) + log_fbar, -1e-12))<<lost_close>>;
  }
}
