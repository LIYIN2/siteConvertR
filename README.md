# siteConvertR

`siteConvertR` is an R package for converting VCF marker coordinates between old and new genome assemblies.

It can convert:

* chromosome names
* marker IDs
* marker positions
* VCF `##contig` header lines

The package supports two conversion directions:

* old genome coordinates to new genome coordinates
* new genome coordinates back to old genome coordinates

This package is designed for genotype data stored in VCF format, especially for marker datasets where old and new marker IDs are connected by a correspondence table.

---

## Features

`siteConvertR` can:

* convert VCF sites from old marker IDs to new marker IDs
* convert VCF sites from new marker IDs back to old marker IDs
* update `CHROM`, `POS`, and `ID` columns in VCF files
* update VCF header lines such as `##contig=<ID=...>`
* remove unmapped markers if requested
* keep unmapped markers if requested
* sort variants by chromosome and position before writing
* support both `.vcf` and `.vcf.gz` input files
* write compressed VCF output

---

## Installation

Install the package from GitHub:

```r
install.packages("devtools")

devtools::install_github("LIYIN2/siteConvertR")
```

Load the package:

```r
library(siteConvertR)
```

---

## Input files

### 1. VCF file

The input file should be a standard VCF file.

Example old-genome VCF records:

```text
#CHROM  POS     ID          REF ALT QUAL FILTER INFO FORMAT sample1 sample2
LG1     1388    LG1_1388    T   .   .    .      PR   GT     ./.     0/0
LG1     24682   LG1_24682   N   .   .    .      PR   GT     ./.     ./.
LG1     40734   LG1_40734   C   A   .    .      PR   GT     1/1     0/1
```

Example new-genome VCF records:

```text
#CHROM  POS      ID            REF ALT QUAL FILTER INFO FORMAT sample1 sample2
chr1    741675   chr1_741675   T   .   .    .      PR   GT     ./.     0/0
chr1    741692   chr1_741692   N   .   .    .      PR   GT     ./.     ./.
```

---

### 2. Mapping table

The marker correspondence table must contain four columns:

```text
old_ID    old_chr    new_chr    new_ID
```

Example:

```text
old_ID          old_chr   new_chr   new_ID
LG10_10004727   LG10      chr10     chr10_741675
LG10_10004744   LG10      chr10     chr10_741692
LG10_10004785   LG10      chr10     chr10_741733
```

Column meaning:

| Column    | Description                       |
| --------- | --------------------------------- |
| `old_ID`  | Marker ID in the old genome       |
| `old_chr` | Chromosome name in the old genome |
| `new_chr` | Chromosome name in the new genome |
| `new_ID`  | Marker ID in the new genome       |

The marker ID should contain the genomic position after an underscore, for example:

```text
LG10_10004727
chr10_741675
```

The package extracts the position from the part after `_`.

---

## Basic usage

### Old genome to new genome

Use:

```r
convert_vcf_sites(
  vcf_file = "old.vcf",
  mapping_file = "mapping.tsv",
  output_file = "new.vcf.gz",
  direction = "old_to_new"
)
```

This converts:

```text
LG10_10004727
```

to:

```text
chr10_741675
```

and updates:

```text
CHROM = chr10
POS   = 741675
ID    = chr10_741675
```

---

### New genome to old genome

Use:

```r
convert_vcf_sites(
  vcf_file = "new.vcf.gz",
  mapping_file = "mapping.tsv",
  output_file = "old.vcf.gz",
  direction = "new_to_old"
)
```

This converts:

```text
chr10_741675
```

back to:

```text
LG10_10004727
```

and updates:

```text
CHROM = LG10
POS   = 10004727
ID    = LG10_10004727
```

---

## Main function

```r
convert_vcf_sites(
  vcf_file,
  output_file = "converted.vcf.gz",
  direction = c("old_to_new", "new_to_old"),
  mapping_file = NULL,
  keep_unmapped = FALSE,
  update_contig_header = TRUE,
  sort_by_position = TRUE
)
```

---

## Arguments

### `vcf_file`

Path to the input VCF file.

The input can be:

```text
input.vcf
input.vcf.gz
```

Example:

```r
vcf_file = "tmp.vcf"
```

---

### `output_file`

Path to the output VCF file.

Recommended output name:

```text
converted.vcf.gz
```

Example:

```r
output_file = "new.vcf.gz"
```

The output file is compressed. Use `zcat` or `zless` to view it.

---

### `direction`

Conversion direction.

Available options:

```r
direction = "old_to_new"
```

Convert old genome sites to new genome sites.

```r
direction = "new_to_old"
```

Convert new genome sites back to old genome sites.

---

### `mapping_file`

Path to the marker correspondence table.

Example:

```r
mapping_file = "mapping.tsv"
```

If the package contains an internal mapping table, `mapping_file` can be omitted:

```r
convert_vcf_sites(
  vcf_file = "old.vcf",
  output_file = "new.vcf.gz",
  direction = "old_to_new"
)
```

If no internal mapping table is available, the user must provide `mapping_file`.

---

### `keep_unmapped`

Whether to keep markers not found in the mapping table.

Default:

```r
keep_unmapped = FALSE
```

This means unmapped markers are removed.

Recommended setting:

```r
keep_unmapped = FALSE
```

This produces a clean output VCF containing only converted markers.

If set to `TRUE`:

```r
keep_unmapped = TRUE
```

unmapped markers are kept in the output VCF. This may produce a mixed-coordinate VCF, where converted markers use the target genome coordinates and unmapped markers still use the original coordinates.

---

### `update_contig_header`

Whether to update VCF `##contig` header lines.

Default:

```r
update_contig_header = TRUE
```

Recommended setting:

```r
update_contig_header = TRUE
```

This changes header lines such as:

```text
##contig=<ID=LG1,length=34872297>
```

to:

```text
##contig=<ID=chr1,length=...>
```

when converting from old genome to new genome.

When converting from new genome back to old genome, it changes:

```text
##contig=<ID=chr1,length=...>
```

back to:

```text
##contig=<ID=LG1,length=...>
```

The contig length is estimated from the maximum marker position in the current VCF file. It may not be the true reference genome chromosome length.

---

### `sort_by_position`

Whether to sort variants by chromosome and position before writing.

Default:

```r
sort_by_position = TRUE
```

Recommended setting:

```r
sort_by_position = TRUE
```

This helps avoid PLINK errors such as:

```text
Error: .bim file has a split chromosome.
```

---

## Example workflow

### Step 1. Convert old VCF to new VCF

```r
library(siteConvertR)

convert_vcf_sites(
  vcf_file = "tmp.vcf",
  mapping_file = "mapping.tsv",
  output_file = "new.vcf.gz",
  direction = "old_to_new",
  keep_unmapped = FALSE,
  update_contig_header = TRUE,
  sort_by_position = TRUE
)
```

---

### Step 2. Check output VCF

In Linux terminal:

```bash
zcat new.vcf.gz | head -40
```

Check contig header:

```bash
zcat new.vcf.gz | grep "^##contig" | head
```

Check variant records:

```bash
zcat new.vcf.gz | grep -v "^##" | head
```

Check chromosome names:

```bash
zcat new.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq
```

For old-to-new conversion, the chromosome names should look like:

```text
chr1
chr2
chr3
chr4
...
```

---

### Step 3. Convert new VCF back to old VCF

```r
convert_vcf_sites(
  vcf_file = "new.vcf.gz",
  mapping_file = "mapping.tsv",
  output_file = "old.vcf.gz",
  direction = "new_to_old",
  keep_unmapped = FALSE,
  update_contig_header = TRUE,
  sort_by_position = TRUE
)
```

Check the result:

```bash
zcat old.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq
```

The chromosome names should look like:

```text
LG1
LG2
LG3
LG4
...
```

---

## PLINK usage

After conversion, the VCF can be used with PLINK.

Recommended workflow:

```bash
plink --vcf new.vcf.gz \
  --make-bed \
  --allow-extra-chr \
  --out new_tmp
```

Then recode if needed:

```bash
plink --bfile new_tmp \
  --recode \
  --allow-extra-chr \
  --out new_ped
```

Avoid doing VCF to PED conversion directly in one step for large or complex files:

```bash
plink --vcf new.vcf.gz --recode --allow-extra-chr --out tmp
```

If PLINK reports:

```text
Error: .bim file has a split chromosome.
```

try sorting the VCF first or use the two-step PLINK workflow above.

---

## Large VCF files

For very large VCF files, such as millions or tens of millions of variants, reading the whole VCF into R may be slow or memory-intensive.

For large datasets, a streaming conversion function may be preferred. The idea is to read and write the VCF line by line instead of loading the entire VCF into memory.

Recommended strategy for very large files:

```bash
# convert VCF using streaming method
# then sort and compress with bcftools/bgzip

bcftools sort converted.vcf.gz -Oz -o converted.sorted.vcf.gz
tabix -p vcf converted.sorted.vcf.gz
```

Install bcftools and tabix with conda:

```bash
conda install -c bioconda bcftools htslib
```

---

## Notes about mapping privacy

If the mapping table is private, do not upload `mapping.tsv` to a public GitHub repository.

If the mapping table is included inside the R package as internal data, users may not see `mapping.tsv` directly, but the data can still be extracted by someone with access to the package.

For private mapping data, recommended options are:

1. Keep the GitHub repository private.
2. Do not include the mapping table in the public package.
3. Require users to provide their own `mapping_file`.

Example using an external private mapping file:

```r
convert_vcf_sites(
  vcf_file = "input.vcf.gz",
  mapping_file = "/path/to/private/mapping.tsv",
  output_file = "output.vcf.gz",
  direction = "old_to_new"
)
```

---

## Output

The output file is a gzipped VCF file.

View it with:

```bash
zcat output.vcf.gz | head
```

or:

```bash
zless output.vcf.gz
```

To create an uncompressed VCF:

```bash
gunzip -c output.vcf.gz > output.vcf
```

---

## Troubleshooting

### Error: mapping table is missing required columns

Make sure your mapping table contains exactly these columns:

```text
old_ID
old_chr
new_chr
new_ID
```

---

### Error: duplicated old_ID or duplicated new_ID

The mapping table should contain one-to-one marker correspondence.

For old-to-new conversion, `old_ID` should not contain duplicated values.

For new-to-old conversion, `new_ID` should not contain duplicated values.

---

### Output VCF looks like unreadable characters

The output is compressed.

Use:

```bash
zcat output.vcf.gz | head
```

Do not use normal `cat` or `less` on compressed files.

---

### Header still contains old contigs

Use:

```r
update_contig_header = TRUE
```

Then check:

```bash
zcat output.vcf.gz | grep "^##contig" | head
```

---

### PLINK reports split chromosome

Use:

```r
sort_by_position = TRUE
```

or sort with bcftools:

```bash
bcftools sort output.vcf.gz -Oz -o output.sorted.vcf.gz
```

Then run PLINK:

```bash
plink --vcf output.sorted.vcf.gz --make-bed --allow-extra-chr --out tmp
```

---

## Citation

If you use this package in a project, please cite the GitHub repository:

```text
siteConvertR: Convert VCF sites between old and new genome coordinates.
https://github.com/LIYIN2/siteConvertR
```

---

## License

This package is released under the MIT License.
