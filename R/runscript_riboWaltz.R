#!/usr/bin/env Rscript

# library(GenomicFeatures,quietly=TRUE,warn.conflicts=FALSE);
# library(riboWaltz,quietly=TRUE,warn.conflicts=FALSE);
# library(ggplot2,quietly=TRUE,warn.conflicts=FALSE);
# library(ggrepel,quietly=TRUE,warn.conflicts=FALSE);
# library(cowplot,quietly=TRUE,warn.conflicts=FALSE);
# library(Biostrings,quietly=TRUE,warn.conflicts=FALSE);

suppressMessages(require(GenomicFeatures));
suppressMessages(require(ggplot2));
suppressMessages(require(ggrepel));
suppressMessages(require(cowplot));
suppressMessages(require(Biostrings));
suppressMessages(require(riboWaltz));

args <- commandArgs(trailingOnly=TRUE);

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
	stop('No cmd line arguments supplied to program!\n', call.=FALSE);
} else if (length(args)==7) {
	bamfile = args[1];
	project <- args[2];
	projectfolder <- args[3];
	path_to_bam <- args[4];
	cellline <- args[5];
	organism <- args[6];
	gtf_source <- args[7];
} else if (length(args)==6) {
	bamfile = args[1];
	project <- args[2];
	projectfolder <- args[3];
	path_to_bam <- args[4];
	cellline <- args[5];
	organism <- args[6];
	gtf_source <- 'Appris_gencode24';
} else if (length(args)==5) {
	bamfile = args[1];
	project <- args[2];
	projectfolder <- args[3];
	path_to_bam <- args[4];
	cellline <- args[5];
	organism <- 'Homo sapiens';
	gtf_source <- 'Appris_gencode24';
} else if (length(args)==4) {
	bamfile = args[1];
	project <- args[2];
	projectfolder <- args[3];
	path_to_bam <- args[4];
	cellline <- 'unspecified';
	organism <- 'Homo sapiens';
	gtf_source <- 'Appris_gencode24';
} else {
	stop('Not all required cmd line arguments supplied to program!\n', call.=FALSE);
}

sample <- gsub('.bam','',bamfile);
if (!dir.exists(paste('riboWaltz_results/',sample,sep=''))) {
	dir.create(paste('riboWaltz_results/',sample,sep=''));
}
# Import patches
cat('Importing patched functions...\n')
source(paste(projectfolder,'patch_riboWaltz.R',sep='/'));

path_to_fasta <- paste(projectfolder,'/appris_principal_isoforms_hg38_gencode_v24_',cellline,'_',project,'.fasta',sep='');
gtf_file <- paste(projectfolder,'/appris_principal_isoforms_hg38_gencode_v24_',cellline,'_',project,'.gtf',sep='');

# Import annotation
cat('Importing annotation...\n');
if (!file.exists(gsub('.gtf','_riboWaltz.RData',gtf_file))) {
	annotation <- create_annotation(gtfpath=gtf_file,dataSource=gtf_source,organism=organism);
	rownames(annotation) <- annotation[,'transcript'];
	save(file=gsub('.gtf','_riboWaltz.RData',gtf_file),annotation);
} else {
	load(file=gsub('.gtf','_riboWaltz.RData',gtf_file));
}

# Import transcript sequences
cat('Importing transcript sequences...\n');
sequences <- Biostrings::readDNAStringSet(path_to_fasta, format = 'fasta', use.names = TRUE);

# Import reads of sample
cat('Importing reads...\n');
df <- as.data.frame(GenomicAlignments::readGAlignments(paste(path_to_bam,bamfile, sep = '/')));
df <- df[, c('seqnames', 'start', 'end', 'width', 'strand')];
colnames(df) <- c('transcript', 'end5', 'end3', 'length', 'strand');
cat(sprintf('reads (total): %f M\n', (nrow(df)/1e+06)));
df <- subset(df, as.character(transcript) %in% rownames(annotation));
cat(sprintf('reads (kept): %f M\n', (nrow(df)/1e+06)));
df <- subset(df, strand == '+');
cat(sprintf('positive strand: %s %%\n', format(round((nrow(df)/nrow(df)) * 100, 2), nsmall = 2)));
cat(sprintf('negative strand: %s %%\n\n', format(round(((nrow(df) - nrow(df))/nrow(df)) * 100, 2), nsmall = 2)));
df$start_pos <- annotation[as.character(df$transcript), 'l_utr5'] + 1;
df$stop_pos <- annotation[as.character(df$transcript), 'l_utr5'] + annotation[as.character(df$transcript), 'l_cds'];
df$start_pos <- ifelse(df$start_pos == 1 & df$stop_pos == 0, 0, df$start_pos);
df <- df[,!(names(df) %in% 'strand')]
reads_list <- list();
reads_list[[sample]] <- df;
rm(df);
suppressMessages(gc());

# QC plot
cat('Plotting read length histogram...');
png(file=paste('riboWaltz_results/',sample,'/readLengthDist_',sample,'.png',sep=''));
length_dist_zoom <- rlength_distr(data=reads_list, sample=sample, cl=99);
length_dist_zoom[['plot']];
suppressMessages(dev.off());

# Determine p-site offsets
cat('Determining p-site offsets...\n');
psite_offset <- psite_onesample(data=reads_list,flanking=12,extremity='auto',plot=TRUE,plotdir='riboWaltz_results');
suppressMessages(dev.off());

# Adding p-site information to reads
reads_psite_list <- psite_info_onesample(data=reads_list,offset=psite_offset,fastapath=path_to_fasta);
rm(reads_list);
suppressMessages(gc());

# More QC plots
cat('Creating QC plots (meta profile etc.)...');
png(file=paste('riboWaltz_results/',sample,'/metaProfile_',sample,'.png',sep=''));
metaprofile <- metaprofile_psite(data=reads_psite_list, annotation=annotation, sample = sample, utr5l = 20, cdsl = 40, utr3l = 20)
metaprofile[['plot']];
suppressMessages(dev.off());

png(file=paste('riboWaltz_results/',sample,'/framesStratified_',sample,'.png',sep=''));
frames_stratified <- frame_psite_length(data=reads_psite_list, sample=sample,region='all', cl=90)
frames_stratified[['plot']]
suppressMessages(dev.off());

png(file=paste('riboWaltz_results/',sample,'/frames_',sample,'.png',sep=''));
frames <- frame_psite(data=reads_psite_list, sample=sample, region='all')
frames[['plot']]
suppressMessages(dev.off());

png(file=paste('riboWaltz_results/',sample,'/metaHeatmap_',sample,'.png',sep=''));
comparison_df <- list()
comparison_df[['subsample_31nt']] <- subset(reads_psite_list[[sample]], length == 31);
comparison_df[['whole_sample']] <- reads_psite_list[[sample]];
names_list <- list('Only_31' = c('subsample_31nt'),'All' = c('whole_sample'));
metaheatmap <- metaheatmap_psite(comparison_df, annotation=annotation, sample = names_list,utr5l = 20, cdsl = 40, utr3l = 20, log=F);
metaheatmap[['plot']];
suppressMessages(dev.off());

# Calculating counts per codon
cat('Calculating read counts per codon...\n');
codon_counts <- riboWaltz:::codon_coverage(data=reads_psite_list,annotation=annotation,sample=sample,psite=TRUE);
write.table(file=paste('riboWaltz_results/codon_counts_',sample,'.tsv',sep=''),codon_counts,sep='\t',quote=FALSE);

# Calculating number of p-sites per CDS
cat('Calculating read counts per CDS...\n');
cds_counts <- psite_per_cds(data=reads_psite_list, annotation=annotation);
write.table(file=paste('riboWaltz_results/cds_counts_',sample,'.tsv',sep=''),cds_counts,sep='\t',quote=FALSE);

# Create codon usage barplot
#cat('Plotting codon usage (barplot)...\n');
#codon_usage_barplot <- codon_usage_psite(data=reads_psite_list, annotation=annotation, sample = sample,fastapath = path_to_fasta);

# Create codon usage scatterplot
#cat('Plotting codon usage (scatterplot)...\n');
#codon_usage_scatter <- codon_usage_psite(data=reads_psite_list, annotation=annotation, sample = sample,fastapath = path_to_fasta, codon_usage = cub_mouse);

cat('------------------------------------------------------\nDone!\n------------------------------------------------------\n');

