## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* The spell checker flags "CFR". It abbreviates "case fatality ratio", spelled out at its first use in the Description.

## Notes for the reviewers

* `cmdstanr` is a suggested package available from the Stan r-universe, listed under `Additional_repositories`. Tests and vignette chunks that fit models run only when `cmdstanr` and CmdStan are installed.
* The examples that fit a model are in `\dontrun{}`: each fit first compiles a Stan model, for several minutes. The other examples run.
