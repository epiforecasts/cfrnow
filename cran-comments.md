## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Notes for the reviewers

* `cmdstanr` is a suggested package available from the Stan r-universe, listed under `Additional_repositories`. Tests and vignette chunks that fit models run only when `cmdstanr` and CmdStan are installed, and are skipped otherwise.
* Examples that fit a model are wrapped in `\dontrun{}`, because each fit compiles a Stan model, which takes several minutes. Examples that do not fit a model run as usual.
