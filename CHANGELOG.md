# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [1.3.1] - 2023-09-16
- modified gridss module (reduced time spent generating metrics plots). Added optional parameter to gridss tasks (additionalParameters)

## [1.3.0] - 2023-08-01
- svprep job scattered by chromosome, RAM allocated based on chromosome size and additional performance-affecting
  parameters are exposed for tuning

## [1.2.0] - 2023-05-24
- Added step for generating names for tumor and normal, basically the same way mutect2 uses

## [1.1.0] - 2023-04-08
- Added additional pre-processing step

## [1.0.1] - 2023-01-09
- Added new gridss aliases (gridss_matched and gridss_matched_by_tumor_group)

## [1.0.0] - 2022-12-13
- modified README.md
- added parameter blacklist
- removed tumor/normal name parameters, added outputFilePrefix
- modified RT 
- Split the workflow in parts to improve speed/reliability

## [Unreleased]

## [0.0.1] - 2022-10-27
### Added
- CHANGELOG.md
- vidarrbuild.json
- vidarrtest-regression.json
- Jenkinsfile

## [0.0.1] - 2022-10-31
### Added
- tests/compare.sh
- tests/calculate.sh

### Changed
- renamed vidarrtest-regression.json to vidarrtest-regression.json.in

