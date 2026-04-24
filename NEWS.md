# metabinR 2.0.0

## Breaking changes

* The three binning entry points — `abundance_based_binning()`,
  `composition_based_binning()`, `hierarchical_binning()` — now return a
  `MetabinResult` S4 object instead of a plain `data.frame`. The v1.x
  tabular layout is available via `as.data.frame(result)`.
* The Java backend no longer exposes a command-line interface; the shaded
  JAR dropped the `commons-cli` dependency. R is the only supported entry
  point.
* Minimum Java runtime is now 17 (previously 8).

## New features

* New `MetabinResult` class with accessors `assignments()`, `nClusters()`,
  `parameters()`, `algorithm()`, plus `show()` and `as.data.frame()`
  methods.
* Binning functions accept, in addition to file paths,
  `Biostrings::DNAStringSet`, `Biostrings::QualityScaledDNAStringSet`, and
  `ShortRead::ShortReadQ` objects. Non-file inputs are staged to a
  tempfile before the Java backend runs.
* `numOfThreads` defaults to `BiocParallel::bpworkers()`.
* `metabinR_jvm_options()` and the `metabinR.jvm.flags` option let users
  customise JVM flags at load time without losing the default G1GC +
  string-deduplication tuning.
* Parameter validation is centralised via `checkmate`; errors surface
  through `cli::cli_abort()` with classed conditions
  (`metabinR_error_*`).

## Internal

* Java sources build via Maven (`java/metabinR/pom.xml`,
  `tools/build-jar.sh`); the shaded JAR is reproducible
  (`project.build.outputTimestamp` pinned).
* Typed JNI entry points on `MTxAB`, `MTxCB`, `MTxABxCB` return
  tab-separated assignments directly to R, removing the previous
  round-trip through on-disk CLI output.

# metabinR 1.5.1 (2024-04-07)

* Preparing for next Bioconductor Release.

# metabinR 1.5.0 (2023-10-29)

* Bump x.y.z version to odd y following creation of RELEASE_3_18 branch.

# metabinR 1.2.0 (2023-04-21)

* Bump x.y.z version to even y prior to creation of RELEASE_3_17 branch.

# metabinR 1.1.0 (2022-11-01)

* Bump x.y.z version to odd y following creation of RELEASE_3_16 branch.

# metabinR 1.0.0 (2022-11-01)

* Bioconductor 3.16 Release. New package **metabinR**, Abundance and
  Compositional Based Binning of Metagenomes.

# metabinR 0.99.3 (2022-10-30)

* Update NEWS.

# metabinR 0.99.2 (2022-10-26)

* Remove citation message on package attach.

# metabinR 0.99.1 (2022-10-11)

* Precheck changes on vignette.

# metabinR 0.99.0 (2022-10-05)

* Submitted to Bioconductor.
