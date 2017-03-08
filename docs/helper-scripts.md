# Helper scripts

These scripts are useful for both the [Rp-Bp](https://github.com/dieterich-lab/rp-bp)
and [B-tea](https://github.com/dieterich-lab/b-tea) pipelines.

## Creating YAML config from sample sheet

The YAML config file used by the pipelines is designed to resemble sample
sheets commonly used, for example, for submitting sequencing data to online
repositories. However, we find such sample sheets tend to be in formats like
Excel and CSV. The `bootstrap-ribo-analysis` script converts such sample sheets
into the YAML config files used by the pipelines.

This script does not include things like the genome reference files; they must
be manually added later. However, it does include them in the output YAML file
with dummy values. Some of the fields included in the output may not be 
necessary depending on the type of analysis. For example, `rnaseq_data` is not
required for Rp-Bp. These values can be deleted or ignored.

Additionally, this script creates symlinks from the original sequencing files to
human-readable files.

```
bootstrap-ribo-analysis <sample_sheet> <out> [--sample-sheet-file-type <file_type>] [--sheet-name <sheet_name>] [--riboseq-sample-types <type_1> ...] [--rnaseq-sample-types <type_1>] [--overwrite]
```

### Command line options

* `sample_sheet`. The sample sheet. See below for more details about expected
  column names and behavior.

* `out`. The YAML config file

* `--sample-sheet-file-type`. By default, the file type of the sample sheet is
  guessed based on the file extension. However, the format can be explicitly
  given using this option. The valid types and associated file extensions are
  as follows.
  
    * `excel`. xls, xlsx

    * `hdf5`. hdf, hdf5, h5, he5

    * `csv`. All other file extensions

* `--sheet-name`. If the sample sheet is an Excel file, then this sheet from 
  the workbook will be treated as the sample sheet. If the sample sheet is from
  an HDF5 file, then this key is used to retrieve it. This is ignored if the
  sample sheet is a CSV file.

* `--riboseq-sample-types`. The `sample_type`s (see below) which indicate a
  sample is a riboseq sample.

* `--rnaseq-sample-types`. The `sample_type`s (see below) which indicate a
  sample is an RNA-seq sample.

* `--overwrite`. If this flag is given, existing files will be deleted and
  rewritten. **N.B.** Be very careful with this flag!

### Columns in sample sheet

The script looks for the following columns in the sample sheet.

#### Required

* `sample_filename`. The absolute path to the sample (fastq) file. 
  Example: `/prj/riechert-riboseq/2016_12_16_riechert/Ribo-seq/CA0EGACXX_EvRe_Nov16_16s007718-1-1_Ibberson_lane51_sequence.txt.gz`

* `condition`. An identifier for the condition of the sample. Example: `sedentary-2wks-wt`

* `sample_type`. The type of this sample. This should match the 
  `--riboseq-sample-types` or `--rnaseq-sample-types`, as appropriate. The
  sample sheet can also contain information about other sample types (ChIP-seq,
  for example), but this script ignores those.

* `location`. The location for the symlinks, concatenated samples, etc. The
  files are placed directly in this directory. Exampe: `/prj/riechert-riboseq/RPF/raw-data/`

#### Optional

* `cell_type`. The type of cell from which the sample comes. Example: `cm`

* `replicate`. A name for this replicate. Example: `mouse-403`

* `lane`. If the reads for a single replicate are spread across multiple lanes,
  then these are often distributed as multiple files. As described below, if 
  multiple lanes are present for a single replicate, they will be concatenated
  as part of the bootstrap process.

 The format of the symlinks and sample names is as follows.
```
<condition>.<sample_type>[.cell-type-<cell_type>][.rep-<replicate>][.lane-<lane>]
```
The optional parts are skipped if they are not present in the sample sheet.
The script first concatenates samples with the same condition, sample type, 
cell type and replicate identifiers (but with different lanes). The 
`biological_replicates` group samples with the same condition, sample type 
and cell type.

The sample sheet can also contain additional columns, but they are ignored.