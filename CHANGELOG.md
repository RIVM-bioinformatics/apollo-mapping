# Changelog

## [0.4.7](https://github.com/RIVM-bioinformatics/apollo-mapping/compare/v0.4.6...v0.4.7) (2025-09-03)


### Bug Fixes

* changing ref genome to latest clade 1 genome ([cb97df0](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/cb97df0d66eb9aa65b653110e9cbbc0998f4cb34))

## [0.4.6](https://github.com/RIVM-bioinformatics/apollo-mapping/compare/v0.4.5...v0.4.6) (2025-05-06)


### Bug Fixes

* Add setuptools v.79.0.1 dependency to yaml ([7b118b2](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/7b118b2f97ec4886320fb88b58087a169d46ed6e))

## [0.4.5](https://github.com/RIVM-bioinformatics/apollo-mapping/compare/v0.4.4...v0.4.5) (2025-04-16)


### Bug Fixes

* update juno library ([7b756dc](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/7b756dc6cc362df3a02fe25283ce349cc84d3994))

## [0.4.4](https://github.com/RIVM-bioinformatics/apollo-mapping/compare/v0.4.3...v0.4.4) (2025-03-27)


### Bug Fixes

* added new reference genome for aspergillus fumigatus ([68afc76](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/68afc766d4642e602fccf13b576744904db491e6))

## [0.4.3](https://github.com/RIVM-bioinformatics/apollo-mapping/compare/v0.4.2...v0.4.3) (2024-06-03)


### Bug Fixes

* disable filter_status in mqc report ([327fd5a](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/327fd5ae34dbabd26089d723cb8cc5ee29e2cf45))

## [0.4.2](https://github.com/RIVM-bioinformatics/apollo-mapping/compare/v0.4.1...v0.4.2) (2023-11-08)


### Bug Fixes

* update master env path ([2314e0a](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/2314e0a9e626d6fe238d562377ae9d6360ca108a))

## [0.4.1](https://github.com/RIVM-bioinformatics/apollo-mapping/compare/v0.4.0...v0.4.1) (2023-10-26)


### Performance Improvements

* set new ref for cauris ([0ada052](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/0ada05238e220b9ae5ecf7d4eed31c86e83daf22))

## [0.4.0](https://github.com/RIVM-bioinformatics/apollo-mapping/compare/v0.3.0...v0.4.0) (2023-10-19)


### Features

* allow running on cauris project ([9eb3e75](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/9eb3e7549bf974a923ac2a55219bf89c66f575c6))

## [0.3.0](https://github.com/RIVM-bioinformatics/apollo-mapping/compare/v0.2.2...v0.3.0) (2023-10-18)


### Miscellaneous Chores

* release 0.3.0 ([f13ed94](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/f13ed946a673e6000f8b1a4c17b6c106395d0b9f))

## [0.2.2](https://github.com/RIVM-bioinformatics/apollo-mapping/compare/v0.2.1...v0.2.2) (2023-07-18)


### Dependencies

* remove anaconda and defaults channel and add nodefaults channel to environments ([0279943](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/0279943ca9631b9020cb4479840e1036cd58d835))

## [0.2.1](https://github.com/RIVM-bioinformatics/apollo-mapping/compare/v0.2.0...v0.2.1) (2023-07-04)


### Bug Fixes

* add HPC resources to last rules ([cd5d1f2](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/cd5d1f2b950d0507ec128d24e8148696024103dc))


### Performance Improvements

* set longer job time ([d8e95fb](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/d8e95fbce8f0919bf6afd4f12b9a881407b317c9))

## [0.2.0](https://github.com/RIVM-bioinformatics/apollo-mapping/compare/v0.1.0...v0.2.0) (2023-06-27)


### Features

* add extra picard qc ([c4fdde3](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/c4fdde38f853e8fb444561959582674e0fb55535))
* customise mqc report ([d9140f6](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/d9140f6c35fad93c6cc9bbbc31224dcf7b5052f1))


### Bug Fixes

* finalise qc variant calling ([6111f81](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/6111f81ce1f530dac6373e07f552096120fd4c7c))
* switch conda activation ([528183a](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/528183a7bac9c09460aa571dc70faf66d357b536))


### Dependencies

* upgrade smk to prevent pip updating smk ([bf081f9](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/bf081f9dde2d4ddbcf6313f957d31533a90f67b8))

## 0.1.0 (2023-05-08)


### Features

* add bbtools pileup and bbtools summary ([5641319](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/564131956c540281faa639239bd7a9c743a9676a))
* add multiqc and some bug fixes: add zstd and zstandard dependency to master environment, run on default with -w 180 for variant calling ([e721806](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/e721806cbcff6d16e69a4688fe144f0060165c19))
* add tables with variant qc ([2da6045](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/2da6045ac3839a904b6927502383d346182b848f))
* add variant evaluation functionality of gatk and modify multiqc file ([da0f4c9](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/da0f4c92651d4070c4282f03baf608a010c22b8d))
* Give default time limit, replace required reference with required species ([cdc1ba0](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/cdc1ba02a4fb1ba12937e0742c00b2be6a9bba44))


### Bug Fixes

* args typo ([22a326e](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/22a326e6efc87fb4fbdd3c5917ffd0898f9c090b))


### Dependencies

* remove mamba env ([7a674c0](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/7a674c0a89120588e54b394a5242816fb70f608b))


### Documentation

* add license ([1dc90b8](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/1dc90b8d688fd26862a3fb417f00e824fa71de8a))
* add readme ([315ffe2](https://github.com/RIVM-bioinformatics/apollo-mapping/commit/315ffe2d92028784722be9af9b39d4b64b71cdf5))
