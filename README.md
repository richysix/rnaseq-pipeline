# RNA-seq analysis pipeline

Our HPC cluster uses Univa Grid Engine to submit jobs and the Modules package with Apptainer containers to control which software is available. The jobs are submitted as shell scripts but I've copied just the commands.

```
# set a variable for the analysis directory
basedir=`pwd`
```

## QC

We use [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to run QC on the initial fastq files and then MultiQC to collect that into a report on all samples together.

### FASTQC

Fastq files are downloaded into a directory called `fastq`  

```
cd $basedir/fastq
fastqc --quiet --threads 12 --noextract *.fq.gz
cd $basedir
```

This uses 12 threads to speed it up, so you will need to make sure that your job has enough cpus available or change the `--threads` option.

### MultiQC

We use [MultiQC](https://multiqc.info/) to aggregate the FASTQC info on indiviual fastq files into one report.

```
find fastq -name *fastqc.zip | sort -V > fastq/multiqc-input.txt
multiqc --file-list fastq/multiqc-input.txt -m fastqc -o fastq -n multiqc_report.html
```

The report is saved as an html file (fastq/multiqc_report.html)

## Mapping

### Download and index reference

We use [STAR](https://github.com/alexdobin/STAR) to map the reads to the reference genome. For information on options etc. see the [manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)

Download the current reference and index it with STAR
```
mkdir -p reference/grcz11
cd reference/
# set release number
version=109
# primary genome assembly
wget ftp://ftp.ensembl.org/pub/release-${version}/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.gz
# gene annotation
wget ftp://ftp.ensembl.org/pub/release-${version}/gtf/danio_rerio/Danio_rerio.GRCz11.${version}.gtf.gz
# unzip
gunzip Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.gz
gunzip Danio_rerio.GRCz11.${version}.gtf.gz

# index genome with STAR
mkdir grcz11 \
STAR \
--outFileNamePrefix grcz11. \
--runThreadN 4 \
--runMode genomeGenerate \
--genomeDir grcz11 \
--genomeFastaFiles Danio_rerio.GRCz11.dna_sm.primary_assembly.fa \
--sjdbGTFfile Danio_rerio.GRCz11.${version}.gtf \
--sjdbOverhang 149
```

The STAR command uses 4 cpus so make sure your job has enough allocated. The `--sjdbOverhang` option should be set to `readlength - 1`.

### Mapping (round 1)

Then STAR is run to map the fastq files for each sample individually. Fastq files are supplied using the `--readFilesIn` option, with read1 file(s), a space and then read2 file(s).  
If a sample has multiple fastq files these can be supplied as a comma-separated list of files.  
e.g. fastq/sample1-1_1.fa.gz,fastq/sample1-2_1.fa.gz fastq/sample1-1_2.fa.gz,fastq/sample1-2_2.fa.gz

If the fastq files are compressed you need to supply the `--readFilesCommand` option to tell STAR how to uncompress the files.

```
cd $basedir
mkdir star1
cd star1
sample='sample_name'
mkdir $sample

STAR \
--runThreadN 1 \
--genomeDir ../reference/grcz11 \
--readFilesIn path/to/fastq1 path/to/fastq2 \
--readFilesCommand zcat \
--outFileNamePrefix $sample/ \
--quantMode GeneCounts \
--outSAMtype BAM SortedByCoordinate
```

### Mapping (round 2)

STAR outputs a file of discovered splice junctions that aren't in the supplied annotation. Then we run STAR again and include those files with the `--sjdbFileChrStartEnd` option

```
cd $basedir
mkdir star2
cd star2
mkdir $sample
STAR \
--runThreadN 1 \
--genomeDir ../reference/grcz11 \
--readFilesIn path/to/fastq1 path/to/fastq2 \
--readFilesCommand zcat \
--outFileNamePrefix $sample/ \
--quantMode GeneCounts \
--outSAMtype BAM SortedByCoordinate \
--sjdbFileChrStartEnd `find ../star1 | grep SJ.out.tab$ | sort | tr '\n' ' '`
```

STAR outputs a file (ReadsPerGene.out.tab) of counts for each gene in the annotation with four columns which correspond to different strandedness options. These files are used as input to the DESeq script.

### Differential Expression

Get Ensembl annotation. This can be downloaded from BioMart. The DESeq2 script expects the annotation to contain the following columns. It shouldn't have a header line.

- Gene (Ensembl Gene ID)
- Chr (Chromosome)
- Start (Gene start coordinate)
- End (Gene end coordinate)
- Strand (Gene orientation 1/-1)
- Biotype (Gene type)
- Name (Gene name)
- Description (Gene description)

Alternatively, if you need a specific version of the annotation, we have a script that uses the Ensembl perl API.

Run the DESeq2 script. The scripts expects a tab-separated file of sample info with either two or three columns and no header line. For two columns, column 1 should be the sample name and column 2 the experimental condition. e.g.  

```
sample1    wt-unexposed
sample2    wt-unexposed
sample3    wt-unexposed
sample4    mut-unexposed
sample5    mut-unexposed
sample6    mut-unexposed
sample7    wt-exposed
sample8    wt-exposed
sample9    wt-exposed
sample10    mut-exposed
sample11    mut-exposed
sample12    mut-exposed
```

It has an option to specify the name of the directory with the count files in and expects the ReadsPerGene.out.tab files to be in directories with the same sample names as those in the samples file.
The script does one comparison, so you'll need to run the script for each comparison that you want to do. e.g.  

```
dir=deseq2-noMn-mut-vs-noMn-wt
exp=mut-unexposed
con=wt-unexposed
Rscript deseq2.R -s $dir/samples.tsv -e $exp -c $con \
-a $basedir/annotation/annotation.txt -d $basedir/star2 -o $dir
```

The script outputs the following results files:

- `all.tsv`: Contains log2 fold changes, pvalues, adjusted pvalues, and normalised counts for all genes
- `sig.tsv`: Same as all.tsv, but just the genes that are DE
- `size-factors.tsv`: the size factors used the normalise the counts based on library size
- `qc.pdf`: QC plots (PCA plot, MA plots)
- `pca.pdf`: plots showing all the principal components that explain more than 1% of the variance
- `PCs.tsv`: results of the PCA. Also for each principal component there is a file of the weighting of each gene for that principal component. (PC-1.tsv etc.)
- `counts.pdf`: count plots for all the DE genes
