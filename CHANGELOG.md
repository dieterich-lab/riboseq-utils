# Change Log
All notable changes to the riboseq utilities will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/), 
and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]
### Changed
- Option `use_chi_square` added to function `get_predicted_orfs`, so selection based on chisq is
    only made if specified in config file.

### Fixed
- Bayes factors (or chisq) filtering applied on top of based-filtered ORFs in `get_predicted_orfs` 
    to avoid Python RuntimeWarning when `bayes_factor_var` is not defined.

## [0.2.6] - 2017-12-08
### Added
- File name for b-tea violin plots
- Consistent filtering for all b-tea fields
- Small changes for Rp-Bp reports

### Updated
- Remove deprecated function calls from `extract_metagene_profiles`, 
    `estimate_metagene_profile_bayes_factors` and `select_periodic_offsets`.

### Added
- Flag added to differentiate between *exons* file (list of exons for any given transcript),
    and *orfs* file (technically, a list of exons per orf). 
    See [Issue #1](https://github.com/dieterich-lab/riboseq-utils/issues/1) for more details.

## [0.2.5] - 2017-10-26
### Updated
- Version specifications for prereqs

## [0.2.4] - 2017-07-17
### Updated
- The `bootstrap-ribo-analysis` script to use `misc.pandas_utils`. Please see
    [Issue #4](https://github.com/dieterich-lab/riboseq-utils/issues/4) for
    more details.

- The `bootstrap-ribo-analysis` script to include replicate and condition
    names and to output the entries in a user-friendly order. Please see
    [Issue #5](https://github.com/dieterich-lab/riboseq-utils/issues/5) for
    more details.

### Fixed
- Removed duplicate definition of `get_te_values` filename

## [0.2.3] - 2017-06-13
### Updated
- All references to `misc.bio_utils` replaced with `bio_utils`
- All references to `misc.bio` replaced with `bio_utils.bio`
- `requirements.txt` to point to the github repositories. Please see
    [Issue #2](https://github.com/dieterich-lab/riboseq-utils/issues/2) for
    more details.

## [0.2.2] - 2017-05-21
### Updated
- The reverse sample maps are now always `_return_key_dicts`.

## [0.2.1] - 2017-05-09
### Fixed
- A reference to `bio.bed6_field_names`

## [0.2.0] - 2017-03-21
### Changed
- The filenames for the motif reports have been completely changed from the
  previous format. In particular, any additional analysis of existing motifs
  will need to use version `0.1.2` of this library.


## [0.1.2] and previous versions

The initial versions have not been documented in the change long. They include
all of the utilities required to run the version of Rp-Bp used in the paper, as
well as things such as filenames which are compatible with early versions of the
B-tea pipeline.
