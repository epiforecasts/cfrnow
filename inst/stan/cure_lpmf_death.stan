/* Death-only mixture-cure log-likelihood for one case.

   Holes (<<name>>) are filled in by epidist_stancode():
     <<family>>         delay family name (lognormal / gamma)
     <<death_pars>>     delay parameter declarations, e.g. "real mu, real sigma"
     <<death_id>>       primarycensored delay distribution id
     <<death_reparam>>  delay parameters in primarycensored's native order
     <<primary_id>>     primarycensored primary (uniform) distribution id

   outcome: 1 = observed death, 3 = resolved non-death, else = censored. */
real cfrnow_<<family>>_lpmf(data int y, <<death_pars>>, real prob,
                            data real outcome, data real pwindow,
                            data real swindow, array[] real primary_params) {
  if (outcome == 1) {
    return log(prob) + primarycensored_lpmf(
        y | <<death_id>>, {<<death_reparam>>}, pwindow, y + swindow, 0.0,
        positive_infinity(), <<primary_id>>, primary_params);
  } else if (outcome == 3) {
    return log1m(prob);
  } else {
    real fbar = primarycensored_cdf(
        y | <<death_id>>, {<<death_reparam>>}, pwindow, 0.0,
        positive_infinity(), <<primary_id>>, primary_params);
    return log1m(prob * fbar);
  }
}
