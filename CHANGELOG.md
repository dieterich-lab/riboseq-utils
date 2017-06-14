# Change Log
All notable changes to the riboseq utilities will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/), 
and this project adheres to [Semantic Versioning](http://semver.org/).

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
