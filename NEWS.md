# carnation

## 1.1.2

- New volcano plot added
- various plot updates
  - MA plot
    - added aspect ratio  and alpha controls
    - new 'refresh' button for reduced reactivity
  - Gene plot
    - add option to show unnormalized counts. This is helpful when the data in the `dds` element is pre-normalized, e.g. in proteomics
  - scatter plot
    - raw p-values added as a plotting option
    - plot/table selection panels now open upon selection
  - aspect ratio control for PCA & scatter plots
- analysis table
  - updated to show all available datasets
  - obtain file info (size, mtime) and show in the table
  - improved selection behavior

- misc fixes
  - correctly reset gene scratchpad when loading data so previously added genes don't show up again.
  - the gene plot function now sets the 'text' aesthetic to `factor()` to stop warnings about statistical transformation
  - `fromList.with.names` now uses base R column indexing instead of `dplyr::relocate`.

## v0.99.10

- add `edgeR`/`limma` support for input files and objects.
- allow addition of degPatterns/pattern analysis results via UI
- improve object validation
- add helper functions to edit carnation default settings

## v0.99.9

- carnation now supports `edgeR` and `limma` output in differential
  expression analysis.
- The scatter plot now allows gene selection directly from the plot
  and viewed in the table. Genes selected in the table can be added
  to the gene scratchpad and tracked across the app.
- The main loading page now shows the currently loaded dataset.
  This prevents accidentally reloading or replacing the current
  data.
- Pattern analysis can now be added to a carnation object via
  the `Load data` module, either in TSV or Rds format.

## v0.99.8

- Add packages to `Suggests::` to address build issues.

## v0.99.7

- Add unit tests
- Address Bioconductor review comments

## v0.99.6

- Made `R CMD check` more efficient by switching multiple module server
  examples to `examplesIf interactive()`.

## v0.99.3

- removed `clusterProfiler` & `DEGreport` from Suggests
- added examples to exported functions

## v0.99.1

- Updated vignette with carnation app screenshots
- Added vignette dependencies to Suggests

## v0.99.0

- carnation submitted to bioconductor for review
