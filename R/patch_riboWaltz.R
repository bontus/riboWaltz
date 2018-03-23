#####
# fixed psite() that can handle single samples as well
psite_onesample = function (data, flanking = 6, samplename = 'sample1', start = TRUE, extremity = "auto", plot = FALSE, plotdir = NULL, plotformat = "png", cl = 99) {
	offset <- NULL
	if (class(data) == 'list') {
		names <- names(data)
	} else {
		names <- samplename
	}
	for (n in names) {
		cat(sprintf("processing %s\n", n))
		if (class(data) == 'list') {
			df <- data[[n]]
		} else {
			df <- data
		}		
		lev <- sort(unique(df$length))
		if (start == T | start == TRUE) {
			base <- 0
			df$site.dist.end5 <- df$end5 - df$start_pos
			df$site.dist.end3 <- df$end3 - df$start_pos
		}
		else {
			base <- -5
			df$site.dist.end5 <- df$end5 - df$stop_pos - base
			df$site.dist.end3 <- df$end3 - df$stop_pos - base
		}
		site.sub <- df[which(df$site.dist.end5 <= -flanking & 
			df$site.dist.end3 >= flanking - 1), ]
		minlen <- min(site.sub$length)
		maxlen <- max(site.sub$length)
		t <- table(factor(site.sub$length, levels = lev))
		offset.temp <- data.frame(length = as.numeric(as.character(names(t))), 
			percentage = (as.vector(t)/sum(as.vector(t))) * 100)
		rownames(offset.temp) <- offset.temp$length
		offset.temp$around_site <- ifelse(offset.temp$percentage == 
			0, "F", "T")
		offset.temp$offset_from_5 <- as.numeric(as.vector(by(site.sub$site.dist.end5, 
			factor(site.sub$length, levels = lev), function(x) names(which.max(table(x))))))
		offset.temp$offset_from_3 <- as.numeric(as.vector(by(site.sub$site.dist.end3, 
			factor(site.sub$length, levels = lev), function(x) names(which.max(table(x))))))
		best.offset.from3.tab <- by(offset.temp$percentage, offset.temp$offset_from_3, 
			function(x) sum(x))
		best.offset.from5.tab <- by(offset.temp$percentage, offset.temp$offset_from_5, 
			function(x) sum(x))
		best.offset.from3 <- names(which.max(best.offset.from3.tab))
		best.offset.from5 <- names(which.max(best.offset.from5.tab))
		if ((extremity == "auto" & best.offset.from3.tab[best.offset.from3] > 
			best.offset.from5.tab[best.offset.from5] & as.numeric(best.offset.from3) <= 
			minlen - 2) | (extremity == "auto" & best.offset.from3.tab[best.offset.from3] <= 
			best.offset.from5.tab[best.offset.from5] & as.numeric(best.offset.from5) > 
			minlen - 1) | extremity == "3end") {
			best.offset <- as.numeric(best.offset.from3)
			line_plot <- "from3"
			cat(sprintf("best offset: %i nts from the 3' end\n", 
				best.offset))
			adj_offset_from_3 <- as.numeric(as.character(do.call(rbind, 
				list(by(site.sub, factor(site.sub$length, levels = lev), 
				  function(x) {
					t <- table(factor(x$site.dist.end3, levels = seq(min(x$site.dist.end3) - 
					  2, max(x$site.dist.end3))))
					t[1:2] <- t[3] + 1
					locmax <- as.numeric(as.character(names(t[which(diff(sign(diff(t))) == 
					  -2)]))) + 1
					adjoff <- locmax[which.min(abs(locmax - best.offset))]
					ifelse(length(adjoff) != 0, as.numeric(as.character(adjoff)), 
					  best.offset)
				  })))))
			adj_offset_from_3[is.na(adj_offset_from_3)] <- best.offset
			offset.temp$adj_offset_from_5 <- -adj_offset_from_3 + 
				offset.temp$length - 1
			offset.temp$adj_offset_from_3 <- adj_offset_from_3
		}
		else {
			if ((extremity == "auto" & best.offset.from3.tab[best.offset.from3] <= 
				best.offset.from5.tab[best.offset.from5] & as.numeric(best.offset.from5) <= 
				minlen - 1) | (extremity == "auto" & best.offset.from3.tab[best.offset.from3] > 
				best.offset.from5.tab[best.offset.from5] & as.numeric(best.offset.from3) > 
				minlen - 2) | extremity == "5end") {
				best.offset <- as.numeric(best.offset.from5)
				line_plot <- "from5"
				cat(sprintf("best offset: %i nts from the 5' end\n", 
				  -best.offset))
				adj_offset_from_5 <- as.numeric(as.character(do.call(rbind, 
				  list(by(site.sub, factor(site.sub$length, levels = lev), 
					function(x) {
					  t <- table(factor(x$site.dist.end5, levels = seq(min(x$site.dist.end5) - 
						2, max(x$site.dist.end5) + 1)))
					  t[1:2] <- t[3] + 1
					  locmax <- as.numeric(as.character(names(t[which(diff(sign(diff(t))) == 
						-2)]))) + 1
					  adjoff <- locmax[which.min(abs(locmax - 
						best.offset))]
					  ifelse(length(adjoff) != 0, as.numeric(as.character(adjoff)), 
						best.offset)
					})))))
				adj_offset_from_5[is.na(adj_offset_from_5)] <- best.offset
				offset.temp$adj_offset_from_5 <- abs(adj_offset_from_5)
				offset.temp$adj_offset_from_3 <- abs(offset.temp$adj_offset_from_5 - 
				  offset.temp$length + 1)
			}
		}
		t <- table(factor(df$length, levels = lev))
		offset.temp$offset_from_5 <- -offset.temp$offset_from_5
		offset.temp$total_percentage <- as.numeric(format(round((as.vector(t)/sum(as.vector(t))) * 
			100, 3), nsmall = 4))
		offset.temp$site_percentage <- as.numeric(format(round(offset.temp$percentage, 
			3), nsmall = 4))
		offset.temp$sample <- n
		offset.temp <- offset.temp[, c("length", "total_percentage", 
			"site_percentage", "around_site", "offset_from_5", 
			"offset_from_3", "adj_offset_from_5", "adj_offset_from_3", 
			"sample")]
		if (start == TRUE) {
			colnames(offset.temp) <- c("length", "total_percentage", 
				"start_percentage", "around_start", "offset_from_5", 
				"offset_from_3", "adj_offset_from_5", "adj_offset_from_3", 
				"sample")
		}
		else {
			colnames(offset.temp) <- c("length", "total_percentage", 
				"stop_percentage", "around_stop", "offset_from_5", 
				"offset_from_3", "adj_offset_from_5", "adj_offset_from_3", 
				"sample")
		}
		if (plot == T || plot == TRUE) {
			options(warn = -1)
			if (length(plotdir) == 0) {
				dir <- getwd()
				plotdir <- paste(dir, "/offset_plot", sep = "")
			}
			if (!dir.exists(plotdir)) {
				dir.create(plotdir)
			}
			minlen <- ceiling(quantile(site.sub$length, (1 - 
				cl/100)/2))
			maxlen <- ceiling(quantile(site.sub$length, 1 - (1 - 
				cl/100)/2))
			for (len in minlen:maxlen) {
				progress <- ceiling(((len + 1 - minlen)/(maxlen - 
				  minlen + 1)) * 25)
				cat(sprintf("\rplotting   %s\r", paste(paste(rep(c(" ", 
				  "<<", "-"), c(25 - progress, 1, progress)), 
				  collapse = ""), " ", as.character(progress * 
				  4), "% ", paste(rep(c("-", ">>", " "), c(progress, 
				  1, 25 - progress)), collapse = ""), sep = "")))
				site.temp <- df[which(df$site.dist.end5 %in% 
				  seq(-len + 1, 0) & df$length == len), ]
				site.tab5 <- as.data.frame(table(factor(site.temp$site.dist.end5, 
				  levels = (-len + 1):(len))))
				site.temp <- df[which(df$site.dist.end3 %in% 
				  seq(0, len - 2) & df$length == len), ]
				site.tab3 <- as.data.frame(table(factor(site.temp$site.dist.end3, 
				  levels = (-len):(len - 2))))
				colnames(site.tab5) <- colnames(site.tab3) <- c("distance", 
				  "reads")
				site.tab5$distance <- as.numeric(as.character(site.tab5$distance))
				site.tab3$distance <- as.numeric(as.character(site.tab3$distance))
				site.tab5$extremity = "5'end"
				site.tab3$extremity = "3'end"
				final.tab <- rbind(site.tab5[site.tab5$distance <= 
				  0, ], site.tab3[site.tab3$distance >= 0, ])
				final.tab$extremity <- factor(final.tab$extremity, 
				  levels = c("5'end", "3'end"))
				p <- ggplot(final.tab, aes(distance, reads, color = extremity)) + 
				  geom_line() + geom_vline(xintercept = seq(-round(len/3) * 
				  3, round(len/3) * 3, 3), linetype = 2, color = "gray90") + 
				  geom_vline(xintercept = 0, color = "gray50") + 
				  geom_vline(xintercept = -offset.temp[as.character(len), 
					"offset_from_5"], color = "#D55E00", linetype = 2, 
					size = 1.1) + geom_vline(xintercept = offset.temp[as.character(len), 
				  "offset_from_3"], color = "#56B4E9", linetype = 2, 
				  size = 1.1) + geom_vline(xintercept = -offset.temp[as.character(len), 
				  "adj_offset_from_5"], color = "#D55E00", size = 1.1) + 
				  geom_vline(xintercept = offset.temp[as.character(len), 
					"adj_offset_from_3"], color = "#56B4E9", 
					size = 1.1) + annotate("rect", ymin = -Inf, 
				  ymax = Inf, xmin = flanking - len, xmax = -flanking, 
				  fill = "#D55E00", alpha = 0.1) + annotate("rect", 
				  ymin = -Inf, ymax = Inf, xmin = flanking - 
					1, xmax = len - flanking - 1, fill = "#56B4E9", 
				  alpha = 0.1) + labs(x = "Distance from start (nt)", 
				  y = "Number of read extremities", title = paste(n, 
					" - length=", len, " nts", sep = ""), color = "Extremity") + 
				  theme_bw(base_size = 20) + scale_fill_discrete("") + 
				  theme(panel.grid.major.x = element_blank(), 
					panel.grid.minor.x = element_blank(), strip.placement = "outside") + 
				  theme(plot.title = element_text(hjust = 0.5))
				if (line_plot == "from3") {
				  p <- p + geom_vline(xintercept = best.offset, 
					color = "black", linetype = 3, size = 1.1) + 
					geom_vline(xintercept = best.offset - len + 
					  1, color = "black", linetype = 3, size = 1.1)
				}
				else {
				  p <- p + geom_vline(xintercept = best.offset, 
					color = "black", linetype = 3, size = 1.1) + 
					geom_vline(xintercept = best.offset + len - 
					  1, color = "black", linetype = 3, size = 1.1)
				}
				p <- p + scale_x_continuous(breaks = seq(-floor(len/5) * 
				  5, floor(len/5) * 5, 5), labels = as.character(seq(-floor(len/5) * 
				  5, floor(len/5) * 5, 5) + base))
				subplotdir <- paste(plotdir, n, sep = "/")
				dir.create(subplotdir)
				save_plot(paste(subplotdir, "/", len, ".", plotformat, 
				  sep = ""), p, base_height = 5, base_aspect_ratio = 3)
			}
			cat(sprintf("\rplotting   %s\n", paste(paste(rep(c(" ", 
				"<<", "-"), c(25 - progress, 1, progress)), collapse = ""), 
				" ", as.character(progress * 4), "% ", paste(rep(c("-", 
				  ">>", " "), c(progress, 1, 25 - progress)), 
				  collapse = ""), sep = "")))
			options(warn = 0)
		}
		offset <- rbind(offset, offset.temp)
		gc();
	}
	return(offset)
}

#####
# fixed psite_info() that can candle single samples as well
psite_info_onesample = function (data, offset, samplename = 'sample1', fastapath = NULL, bsgenome_dp = NULL, txdb = NULL, granges = FALSE) {
	if (length(fastapath) != 0 & length(bsgenome_dp) != 0) {
		warning("fastapath and bsgenome_dp are both specified. Only fastapath will be considered\n")
		bsgenome_dp = NULL
	}
	if (length(bsgenome_dp) != 0 & length(txdb) == 0) {
		cat("\n")
		stop("\nERROR: txdb is not specified \n\n")
	}
	if (length(fastapath) != 0 | length(bsgenome_dp) != 0) {
		if (length(fastapath) != 0) {
			cat("importing FASTA file\n\n")
			sequences <- Biostrings::readDNAStringSet(fastapath, 
				  format = "fasta", use.names = TRUE)
		}
		else {
			cat("preparing sequences based on bsgenome and txdb\n\n")
			if (length(bsgenome_dp) != 0) {
				if (bsgenome_dp %in% installed.genomes()) {
					library(bsgenome_dp, character.only = TRUE)
				}
				else {
					source("http://www.bioconductor.org/biocLite.R")
					biocLite(bsgenome_dp, suppressUpdates = TRUE)
					library(bsgenome_dp, character.only = TRUE)
				}
			}
			sequences <- extractTranscriptSeqs(get(bsgenome_dp), 
				txdb, use.names = T)
		}
	}
	if (class(data) == 'list') {
		names <- names(data)
	} else {
		names <- samplename
	}
	for (n in names) {
		cat(sprintf("processing %s\n", n))
		df <- data[[n]]
		cat("adding p-site position\n")
		df <- dplyr::left_join(df, offset[which(offset$sample == n), c("length", "adj_offset_from_3")], 
			by = "length")
		colnames(df)[colnames(df) == "adj_offset_from_3"] <- "psite"
		df$psite <- df$end3 - df$psite
		df <- df[, c("transcript", "end5", "psite", "end3", "length", 
			"start_pos", "stop_pos")]
		df$psite_from_start <- ifelse(df$stop_pos == 0, 0, df$psite - 
			df$start_pos)
		df$psite_from_stop <- ifelse(df$stop_pos == 0, 0, df$psite - 
			df$stop_pos)
		cat("adding region\n")
		df$psite_region <- ifelse(df$stop_pos == 0, NA, ifelse(df$psite_from_start >= 
			0 & df$psite_from_stop <= 0, "cds", ifelse(df$psite_from_start < 
			0 & df$psite_from_stop < 0, "5utr", "3utr")))
		if (length(fastapath) != 0 | length(bsgenome_dp) != 0) {
			df$psite_codon <- as.character(subseq(sequences[as.character(df$transcript)], 
				start = df$psite, end = df$psite + 2))
		}
		if (granges == T || granges == TRUE) {
			df <- GenomicRanges::makeGRangesFromDataFrame(df, 
				keep.extra.columns = TRUE, ignore.strand = TRUE, 
				seqnames.field = c("transcript"), start.field = "end5", 
				end.field = "end3", strand.field = "strand", 
				starts.in.df.are.0based = FALSE)
			strand(df) <- "+"
		}
		data[[n]] <- df
		gc();
	}
	if (granges == T || granges == TRUE) {
		data <- GenomicRanges::GRangesList(data)
	}
	return(data)
}

