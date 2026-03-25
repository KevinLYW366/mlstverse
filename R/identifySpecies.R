#' read SAM format file
#'
#' @param filename
#'
#' @return
#' @export
#'
#' @examples
readSAM <- function(filename) {
  sam <- suppressMessages(suppressWarnings(readr::read_delim(filename, delim="\t", col_names=F, quote='', progress=F, comment="@")))
  colnames(sam) <- c("qname", # Query template NAME
                     "flag",  # bitwise FLAG
                     "rname", # References sequence NAME
                     "pos",   # 1 - based leftmost mapping POSition
                     "mapq",  # MAPping Quality
                     "cigar", # CIGAR String
                     "rnext", # Ref. name of the mate/next read
                     "pnext", # Position of the mate/next read
                     "tlen",  # observed Template LENgth
                     rep("", ncol(sam) - 9))
  return(sam[,1:11])
}

#' create an input query for MLST score calculation
#'
#' @param sam
#'
#' @return
#' @export
#'
#' @examples
buildDepthSummary <- function(depth, len.loci) {
  if (nrow(depth) == 0) {
    return(data.frame(locus_tag=character(),
                      depth=numeric(),
                      coverage=numeric(),
                      coverRatio=numeric(),
                      len=numeric(),
                      mapped=integer(),
                      locus=character(),
                      stringsAsFactors=FALSE))
  }

  seqnames <- as.character(depth$seqnames)
  seq_ids <- unique(seqnames)
  split_depth <- split(depth$count, seqnames)
  split_depth <- split_depth[seq_ids]
  mapped <- lengths(split_depth)
  total_depth <- vapply(split_depth, sum, numeric(1))
  seq_len <- as.numeric(len.loci[seq_ids])

  data.frame(locus_tag=sub("^[^_]*_", "", seq_ids),
             depth=total_depth,
             coverage=total_depth/seq_len,
             coverRatio=as.numeric(mapped/seq_len),
             len=seq_len,
             mapped=as.integer(mapped),
             locus=sub("_.*$", "", seq_ids),
             stringsAsFactors=FALSE)
}

createQuery <- function(loci, depth, len.loci, threads=1) {
  continueSF <- TRUE
  if (!snowfall::sfIsRunning()) {
    suppressMessages(snowfall::sfInit(parallel=T, cpus=threads, useRscript=T))
    continueSF <- FALSE
  }
  depth_summary <- buildDepthSummary(depth, len.loci)
  split_summary <- split(depth_summary[c("locus_tag", "depth", "coverage", "coverRatio", "len", "mapped")],
                         depth_summary$locus)
  empty_query <- data.frame(locus_tag=character(),
                            depth=numeric(),
                            coverage=numeric(),
                            coverRatio=numeric(),
                            len=numeric(),
                            mapped=integer(),
                            stringsAsFactors=FALSE)
  query <- snowfall::sfLapply(loci, function(l, split_summary, empty_query) {
    q <- split_summary[[l]]
    if (is.null(q)) {
      return(empty_query)
    }
    rownames(q) <- NULL
    q
  }, split_summary, empty_query)
  if (!continueSF) {
    snowfall::sfStop()
  }
  names(query) <- loci
  return(query)
}

#' calculate MLST score
#'
#' @param threads
#' @param query
#' @param mlstdb
#'
#' @return
#' @export
#'
#' @examples
buildQueryLookup <- function(query) {
  lapply(query, function(q) {
    if (nrow(q) == 0) {
      return(list(cover_ratio=numeric(),
                  mapped=numeric(),
                  coverage=numeric(),
                  locus_tag=character()))
    }

    tags <- as.character(q$locus_tag)
    cover_ratio <- q$coverRatio
    names(cover_ratio) <- tags
    mapped <- q$mapped
    names(mapped) <- tags
    coverage <- q$coverage
    names(coverage) <- tags

    list(cover_ratio=cover_ratio,
         mapped=mapped,
         coverage=coverage,
         locus_tag=tags)
  })
}

calcMLSTScore <- function(query,
                          loci,
                          mlstdb=mlstverse.NTM.db,
                          method="default",
                          normalize=TRUE,
                          threads=1) {
  continueSF <- TRUE
  if (!snowfall::sfIsRunning()) {
    suppressMessages(snowfall::sfInit(parallel=T, cpus=threads, useRscript=T))
    continueSF <- FALSE
  }
  query_lookup <- buildQueryLookup(query)
  score_limit <- rowSums(!vapply(mlstdb[, loci, drop=FALSE], function(col) {
    vapply(col, function(x) {-1 %in% x}, logical(1))
  }, logical(nrow(mlstdb))))

  g <- function(db_entry, query_lookup, score_limit, method) {
    db_entry <- lapply(db_entry, "[[", 1)
    f <- function(l, db_entry, query_lookup, method) {
      if (-1 %in% db_entry[[l]] | length(query_lookup[[l]]$locus_tag) == 0) {
        return(0)
      }
      found <- db_entry[[l]] %in% query_lookup[[l]]$locus_tag
      if (any(found)) {
        if (method=="default") {
          matched <- db_entry[[l]][found]
          return(mean(query_lookup[[l]]$cover_ratio[matched]) / length(db_entry[[l]]))
        } else if (method=="sensitive") {
          return(sum(found) / length(db_entry[[l]]))
        }
      } else {
        return(0)
      }
    }
    tmp <- sapply(names(db_entry), f, db_entry, query_lookup, method)
    if (score_limit > 0) {
      if (normalize) {
        return(sum(tmp, na.rm=T) / score_limit)
      } else {
        return(sum(tmp, na.rm=T))
      }
    } else {
      return(0)
    }
  }
  scores <- snowfall::sfApply(mlstdb[,loci], 1, g, query_lookup, score_limit, method)
  if (!continueSF) {
    snowfall::sfStop()
  }
  return(scores)
}

#' get read count
#'
#' @param entry
#' @param query
#' @param loci
#' @param method
#' @param fill
#'
getCounts <- function(entry, query, loci, method="coverage", fill=TRUE) {
  x <- c()
  for (l in loci) {
    isFound <- FALSE
    if (nrow(query[[l]]) > 0) {
      tmp <- subset(query[[l]], locus_tag %in% as.character(entry[[l]]))
      if (nrow(tmp) > 0) {
        isFound <- TRUE
        if (method=="coverage") {
          x <- c(x, mean(tmp$coverage)/nrow(tmp))
        } else if (method=="mapped") {
          x <- c(x, mean(tmp$mapped)/nrow(tmp))
        }
      }
    }
    if (!isFound) {
      x <- c(x, NA)
    }
  }
  if (fill) {
    i <- is.na(x) & !sapply(entry[loci], function(x) {-1%in%x})
    x[i] <- 0
  }
  return(x)
}

getCountsFast <- function(entry, query_lookup, loci, method="coverage", fill=TRUE) {
  x <- numeric(length(loci))
  is_missing <- logical(length(loci))

  for (idx in seq_along(loci)) {
    l <- loci[idx]
    lookup <- query_lookup[[l]]
    entry_tags <- as.character(entry[[l]])
    found <- entry_tags %in% lookup$locus_tag

    if (any(found)) {
      matched <- entry_tags[found]
      values <- if (method == "coverage") lookup$coverage[matched] else lookup$mapped[matched]
      x[idx] <- mean(values) / length(matched)
    } else {
      x[idx] <- NA_real_
      is_missing[idx] <- TRUE
    }
  }

  if (fill) {
    valid_loci <- !vapply(entry[loci], function(x) {-1 %in% x}, logical(1))
    x[is_missing & valid_loci] <- 0
  }
  x
}


obtainConsensus <- function(bam, query, locus) {

}


#' Title
#'
#' @param filenames
#' @param mlstdb database for MLST (default: mlstverse.NTM.db)
#' @param min_depth only use genes larger than the minimum reads depth (default: 0)
#' @param min_ratio only use genes larger than the ratio to maximum reads depth (default: 0.1)
#' @param th.pvalue filter by threshold value in Kolmogorov–Smirnov test (default: 0.05)
#' @param th.score filter by threshold value in MLST score (default: 0.1)
#' @param threads number of threads (default: 1)
#' @param normalize boolean, if TRUE, normalize coverage (default: TRUE)
#' @param samfile return value of readSAM() (default: NULL)
#'
#' @return
#' @export
#'
#' @examples
mlstverse <- function(filenames,
                      mlstdb=mlstverse.Mycobacterium.db::mlstverse.Mycobacterium.db,
                      min_depth=0,
                      min_ratio=0.1,
                      th.pvalue=0.05,
                      th.score=0.1,
                      threads=1,
                      normalize=TRUE,
                      samfile=NULL,
                      method="default") {
  suppressWarnings(suppressMessages(snowfall::sfInit(parallel=T, cpus=threads, useRscript=T)))
  loci <- colnames(mlstdb)[grep("Locus_[0-9]+", colnames(mlstdb))]

  score <- list()
  query <- list()
  for (filename in filenames) {
    cat(paste("Processing", filename, "\n"))
    cat(paste("  Loading bam file...\n"))
    if (is.null(samfile)) {
      pileup.param <- Rsamtools::PileupParam(max_depth=10000,
                              min_base_quality=0,
                              min_mapq=0,
                              min_nucleotide_depth=0,
                              min_minor_allele_depth=0,
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=TRUE,
                              include_deletions=TRUE,
                              include_insertions=TRUE,
                              left_bins=NULL,
                              query_bins=NULL,
                              cycle_bins=NULL)
      depth <- Rsamtools::pileup(filename, pileupParam=pileup.param)
      len.loci <- Rsamtools::scanBamHeader(filename)[[1]]$targets
    } else {
      sam <- samfile
    }
    cat(paste("  Generating query...\n"))
    query[[filename]] <-
      lapply(createQuery(loci, depth, len.loci),
             function(q) {
               if (nrow(q) > 0) {
                 min_d <- max(max(q$depth)*min_ratio, min_depth)
                 subset(q, min_d < depth)
                 } else {
                   return(q)
                   }
               })
    cat(paste("  Calculating MLST score...\n"))
    # sfApply() can return a 1-column matrix for some inputs; flatten it so
    # downstream logical indexing against mlstdb always uses a plain vector.
    results <- as.numeric(calcMLSTScore(query[[filename]], loci, mlstdb=mlstdb, threads=threads, method=method, normalize=normalize))
    query_lookup <- buildQueryLookup(query[[filename]])

    if (normalize) {
      i <- results > th.score
    } else {
      i <- results > max(results)*th.score
    }
    speciesNames <- unique(apply(
      mlstdb[i, c("genus", "species")], 1, paste, collapse="|"))
    if (length(speciesNames) == 0) {
      score[[filename]] <- dplyr::data_frame(genus=character(),
                                             species=character(),
                                             strain=character(),
                                             score=numeric(),
                                             p.value=numeric(),
                                             mean=numeric(),
                                             var=numeric())
      next
    }

    tmp <- list()

    tmp$mean <- c()
    tmp$var <- c()
    tmp$dist <- c()
    tmp$score <- c()
    tmp$pvalue <- c()
    tmp$strains <- c()
    cat("  Starting Kolmogorov–Smirnov test...\n")
    for (speciesName in speciesNames) {
      cat(paste("    testing for", speciesName, "\n"))
      x <- strsplit(speciesName, "\\|")[[1]]
      if (is.na(x[2])) {
        j <- which(i & mlstdb$genus==x[1] & mlstdb$species=="")
      } else {
        j <- which(i & mlstdb$genus==x[1] & mlstdb$species==x[2])
      }
      db_entries <- mlstdb[j, loci, drop=FALSE]
      q.dist <- vapply(seq_len(nrow(db_entries)), function(idx) {
        getCountsFast(db_entries[idx, , drop=FALSE], query_lookup, loci, fill=TRUE)
      }, numeric(length(loci)))
      if (!is.matrix(q.dist)) {
        q.dist <- matrix(q.dist, nrow=length(loci), dimnames=list(NULL, NULL))
      }

      q.mean <- apply(q.dist, 2, mean, na.rm=T)
      q.var <- apply(q.dist, 2, var, na.rm=T)
      snowfall::sfExport("q.dist", "q.mean", "q.var")

      tmp$ks.test <- snowfall::sfLapply(
        1:ncol(q.dist),
        function(i) {
          m <- q.mean[i]
          v <- q.var[i]
          fit <- try(MASS::fitdistr(round(na.omit(q.dist[,i])),
                          densfun="negative binomial",
                          start=list(mu=m, size=m^2/(m+v))), silent=FALSE)
          if (class(fit) == "try-error") {
            return(NA)
          } else {
            ks.test(x=q.dist[,i],
                    y="pnbinom",
                    size=fit$estimate[2],
                    prob=fit$estimate[2]/fit$estimate[1])
          }
        }
      )
      
      tmp$ks.pvalue <- sapply(tmp$ks.test, function(x) {
        if(all(is.na(x))) {
          return(NA)
        } else {
          return(x[["p.value"]])
        }
      })
      table.score <- data.frame(index=j, score=results[j], pvalue=tmp$ks.pvalue)
      if (th.pvalue == 0 | all(is.na(table.score$pvalue))) {
        table.score.p <- table.score
      } else if (any(table.score$pvalue > th.pvalue, na.rm=T)) {
        table.score.p <- subset(table.score, pvalue > th.pvalue)
      } else {
        table.score.p <- subset(table.score, pvalue==max(table.score$pvalue, na.rm=T))
      }
      table.score.p.s <- subset(table.score.p, score == max(table.score.p$score))
      l <- table.score.p.s$index[order(table.score.p.s$pvalue, decreasing=T)[1]]
      if (is.na(l)) {
        tmp$mean <- c(tmp$mean, NA)
        tmp$var <- c(tmp$var, NA)
        tmp$dist <- c(tmp$dist, list(NA))
        tmp$pvalue <- c(tmp$pvalue, NA)
        tmp$score <- c(tmp$score, results[l])
        tmp$strains <- c(tmp$strains, gsub("strain=", "", lapply(strsplit(mlstdb$notes[l], ","), "[", 3)))
      } else {
        tmp$mean <- c(tmp$mean, q.mean[j == l])
        tmp$var <- c(tmp$var, q.var[j == l])
        tmp$dist <- c(tmp$dist, list(q.dist[, j == l]))
        tmp$pvalue <- c(tmp$pvalue, tmp$ks.pvalue[j == l])
        tmp$score <- c(tmp$score, results[l])
        tmp$strains <- c(tmp$strains, gsub("strain=", "", lapply(strsplit(mlstdb$notes[l], ","), "[", 3)))
      }
    }
    if (any(is.na(tmp$strains))) {
      tmp$strains[is.na(tmp$strains)] <- ""
    }
    score[[filename]] <-
      dplyr::data_frame(genus=gsub("\\|.*", "", speciesNames),
                 species=gsub(".*\\|", "", speciesNames),
                 strain=tmp$strains,
                 score=tmp$score,
                 p.value=tmp$pvalue,
                 mean=tmp$mean,
                 var=tmp$var)
    if (nrow(score[[filename]]) > 0) {
      score[[filename]] <- score[[filename]][order(score[[filename]]$p.value, decreasing=T),]
      score[[filename]] <- score[[filename]][order(score[[filename]]$score, decreasing=T),]
    }
  }
  snowfall::sfStop()
  return(list(query=query, score=score))
}
