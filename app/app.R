# =============================================================================
# MS-DIAL -> Agilent dMRM method builder (positive mode, identities from negative)
#
# Run:  shiny::runApp("app")        Packages: shiny, DT
#
# Input : MS-DIAL alignment result export(s) (tab-delimited .txt), positive
#         (required) and negative (optional, used only to identify acyl chains).
# Output: Agilent MassHunter dynamic-MRM transition list (.csv) + annotated table.
#
# Source: https://github.com/sarahehancock/MSDIAL-dMRM-builder
#
# Copyright (C) 2026 Sarah Hancock
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, version 3. It is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the LICENSE file for details.
# =============================================================================

library(shiny)
library(DT)

REPO_URL <- "https://github.com/sarahehancock/MSDIAL-dMRM-builder"

options(shiny.maxRequestSize = 500 * 1024^2)

# ---- masses -----------------------------------------------------------------
mH <- 1.00782503223; mC <- 12; mN <- 14.00307400443
mO <- 15.99491461957; mP <- 30.97376199842; mE <- 0.000548579909
PROTON <- mH - mE
H2O <- 2 * mH + mO
NH3 <- mN + 3 * mH
ADDUCT <- c("[M+H]+" = PROTON, "[M+NH4]+" = NH3 + PROTON,
            "[M+Na]+" = 22.98976928 - mE, "[M+H-H2O]+" = PROTON - H2O)

PCHO   <- 184.07332   # phosphocholine C5H15NO4P+
CHOL   <- 104.10699   # choline
NL_PE  <- 141.01909   # phosphoethanolamine
NL_PS  <- 185.00892   # phosphoserine
NL_PI  <- 260.02972   # inositol phosphate
NL_PG  <- 172.01367   # glycerophosphate
NL_PA  <- 97.97690    # H3PO4
HEX    <- 162.05282
C3H6O2 <- 74.03678
CHOLESTADIENE <- 27 * mC + 45 * mH - mE   # 369.3516

fa_mass   <- function(c, d) c * mC + (2 * c - 2 * d) * mH + 2 * mO
acylium   <- function(c, d) fa_mass(c, d) - mO - mH - mE
lcb_ion   <- function(c, d, ox = 2) c * mC + (2 * c + 3 - 2 * d) * mH + mN + ox * mO - 2 * H2O + PROTON
dg_like   <- function(c, d) 3 * mC + 8 * mH + 3 * mO + c * mC + (2 * c - 2 * d) * mH + 4 * mO - 2 * H2O - mO - mH - mE
plas_sn1  <- function(c, dO) (c + 2) * mC + (2 * c + 7 - 2 * (dO - 1)) * mH + mN + 4 * mO + mP - mE

# plausible fatty acyl chains considered when mining MS2 spectra
FA_TAB <- do.call(rbind, lapply(8:26, function(c) {
  maxd <- if (c < 14) 0 else if (c < 16) 1 else if (c < 18) 2 else if (c < 20) 4 else if (c < 22) 5 else 6
  data.frame(c = c, d = 0:maxd)
}))
FA_TAB$tok <- paste0(FA_TAB$c, ":", FA_TAB$d)
FA_TAB$m   <- fa_mass(FA_TAB$c, FA_TAB$d)

NCHAIN <- c(TG = 3, DG = 2, PC = 2, PE = 2, PS = 2, PI = 2, PG = 2, PA = 2, CL = 2,
            Cer = 2, SM = 2, HexCer = 2, Hex2Cer = 2)
SLASH_CLASSES <- c("Cer", "SM", "HexCer", "Hex2Cer")

default_ce <- function() data.frame(
  Class = c("PC", "PC-O", "LPC", "SM", "PE", "PE-O", "LPE", "PS", "PI", "PG", "PA",
            "DG", "TG", "TG-O", "CE", "Cer", "HexCer", "CL", "NAE", "CAR", "other"),
  CE_headgroup = c(28, 28, 28, 28, 20, 20, 18, 20, 20, 18, 18, 15, 25, 25, 12, 28, 32, 30, 15, 28, 25),
  CE_chain     = c(32, 32, 28, 28, 28, 22, 18, 28, 28, 25, 25, 20, 25, 25, 12, 28, 32, 30, 15, 28, 25),
  stringsAsFactors = FALSE)

# ---- reading MS-DIAL --------------------------------------------------------
read_msdial <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  hdr <- grep("^\"?Alignment ID", lines)[1]
  if (is.na(hdr)) stop("No 'Alignment ID' header row found - is this an MS-DIAL alignment result export?")
  sep <- if (grepl("\t", lines[hdr])) "\t" else ","
  df <- read.table(text = lines[hdr:length(lines)], sep = sep, header = TRUE,
                   quote = if (sep == ",") "\"" else "", comment.char = "",
                   check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE)
  need <- c("Average Rt(min)", "Average Mz", "Metabolite name", "Adduct type", "MS/MS spectrum")
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("Missing column(s): ", paste(miss, collapse = ", "))
  # intensity columns = everything after the MS/MS column that has a file type (Sample / QC / Blank / Standard)
  k <- match("MS/MS spectrum", names(df)); idx <- if (k < ncol(df)) (k + 1):ncol(df) else integer(0); type <- rep("Sample", length(idx))
  ft <- grep("File type", lines[seq_len(hdr - 1)])
  if (length(ft) && length(idx)) {
    v <- strsplit(lines[ft[1]], sep, fixed = TRUE)[[1]]; v <- c(v, rep(NA, max(0, ncol(df) - length(v))))[idx]
    ok <- !is.na(v) & nzchar(v) & v != "NA" & v != "File type"; idx <- idx[ok]; type <- v[ok]
  }
  attr(df, "samples") <- data.frame(idx = idx, name = names(df)[idx], type = type, stringsAsFactors = FALSE)
  df
}

parse_spec <- function(s) {
  if (is.na(s) || !nzchar(s) || s == "null") return(matrix(numeric(0), ncol = 2))
  p <- strsplit(strsplit(trimws(s), "[ ;]+")[[1]], ":", fixed = TRUE)
  p <- p[lengths(p) == 2]
  m <- matrix(suppressWarnings(as.numeric(unlist(p))), ncol = 2, byrow = TRUE)
  m[stats::complete.cases(m), , drop = FALSE]
}

pk <- function(spec, mz, tol) {
  if (!nrow(spec) || !length(mz)) return(numeric(length(mz)))
  vapply(mz, function(m) {
    i <- abs(spec[, 1] - m) <= tol
    if (isTRUE(any(i))) max(spec[i, 2]) else 0
  }, 0)
}

# ---- lipid name parsing -----------------------------------------------------
tok_parse <- function(t) {
  m <- regmatches(t, regexec("^(O-|P-)?(\\d+):(\\d+)(;O(\\d*))?$", t))[[1]]
  if (!length(m)) return(NULL)
  ox <- if (nzchar(m[5])) { if (nzchar(m[6])) as.integer(m[6]) else 1L } else 0L
  d <- as.integer(m[4]); if (m[2] == "P-") d <- d + 1L
  list(ether = nzchar(m[2]), c = as.integer(m[3]), d = d, ox = ox)
}
tok_str <- function(x) paste0(if (x$ether) "O-" else "", x$c, ":", x$d,
                              if (x$ox > 0) paste0(";O", if (x$ox > 1) x$ox else "") else "")

parse_part <- function(p) {
  p <- trimws(p)
  p <- sub("^([A-Za-z][A-Za-z0-9]*)\\((.*)\\)$", "\\1 \\2", p)          # DG(15:0/15:0/0:0)
  m <- regmatches(p, regexec("^([A-Za-z][A-Za-z0-9]*)\\s+(.+)$", p))[[1]]
  if (!length(m)) return(NULL)
  toks <- strsplit(m[3], "[_/]")[[1]]
  toks <- toks[toks != "0:0"]
  pr <- lapply(toks, tok_parse)
  if (!length(pr) || any(vapply(pr, is.null, TRUE))) return(NULL)
  list(base = m[2], chains = pr)
}

species_key <- function(chains, base) {
  eth <- vapply(chains, function(x) x$ether, TRUE)
  if (!(base %in% SLASH_CLASSES)) {
    o <- order(!eth, vapply(chains, function(x) x$c, 1), vapply(chains, function(x) x$d, 1))
    chains <- chains[o]
  }
  vapply(chains, tok_str, "")
}

species_label <- function(base, toks, plasmalogen = FALSE) {
  lab <- toks
  if (plasmalogen) lab <- vapply(toks, function(t) {
    x <- tok_parse(t)
    if (x$ether && x$d >= 1) paste0("P-", x$c, ":", x$d - 1) else t
  }, "")
  paste(base, paste(lab, collapse = if (base %in% SLASH_CLASSES) "/" else "_"))
}

parse_lipid <- function(name) {
  out <- list(ok = FALSE, name = name, cls = "other", base = NA, key = NA, sumname = name,
              species = NULL, C = NA, D = NA, ox = 0, ether = FALSE)
  out$istd <- grepl("\\[IS\\]|\\(IS\\)|ISTD", name)
  nm <- trimws(gsub("\\[IS\\]|\\(IS\\)", "", name))
  nm <- sub("^(w/o MS2:|low score:)\\s*", "", nm)
  out$labelled <- grepl("(^|[ (\\-])d\\d+([ )\\-]|$)|13C", nm)
  nm <- trimws(gsub("\\(d\\d+\\)", "", nm)); nm <- sub("^d\\d+[- ]", "", nm)
  parts <- strsplit(nm, "|", fixed = TRUE)[[1]]
  out$sumname <- trimws(parts[1])
  p1 <- parse_part(parts[1]); if (is.null(p1)) return(out)
  base <- p1$base; if (base %in% c("Ac", "ACar")) base <- "CAR"
  eth <- any(vapply(p1$chains, function(x) x$ether, TRUE))
  out$ok <- TRUE; out$base <- base; out$ether <- eth
  out$cls <- paste0(base, if (eth) "-O" else "")
  out$C <- sum(vapply(p1$chains, function(x) x$c, 1))
  out$D <- sum(vapply(p1$chains, function(x) x$d, 1))
  out$ox <- sum(vapply(p1$chains, function(x) x$ox, 1))
  out$key <- paste(out$cls, out$C, out$D, out$ox)
  nexp <- if (base %in% names(NCHAIN)) NCHAIN[[base]] else 1
  pm <- parse_part(parts[length(parts)])
  if (!is.null(pm) && length(pm$chains) == nexp && nexp > 1)
    out$species <- list(species_key(pm$chains, base))
  out
}

# ---- fragment rules (positive mode) -----------------------------------------
class_products <- function(lip, MH) {
  b <- lip$base
  r <- switch(b,
    PC = , SM = list(c(PCHO, "phosphocholine m/z 184")),
    LPC = list(c(PCHO, "phosphocholine m/z 184"), c(CHOL, "choline m/z 104"), c(MH - H2O, "[M+H-H2O]+")),
    PE = , LPE = list(c(MH - NL_PE, "NL 141 (phosphoethanolamine)")),
    PS = , LPS = list(c(MH - NL_PS, "NL 185 (phosphoserine)")),
    PI = , LPI = list(c(MH - NL_PI, "NL 260 (inositol phosphate)")),
    PG = , LPG = list(c(MH - NL_PG, "NL 172 (glycerophosphate)")),
    PA = list(c(MH - NL_PA, "NL 98 (H3PO4)")),
    DG = , MG = list(c(MH - H2O, "[M+H-H2O]+")),
    CE = list(c(CHOLESTADIENE, "cholestadiene m/z 369")),
    Cer = list(c(MH - H2O, "[M+H-H2O]+")),
    HexCer = list(c(MH - HEX - H2O, "[M+H-Hex-H2O]+"), c(MH - H2O, "[M+H-H2O]+")),
    NAE = list(c(MH - H2O, "[M+H-H2O]+"), c(62.06004, "ethanolamine m/z 62")),
    CAR = list(c(85.02841, "carnitine fragment m/z 85")),
    list())
  if (!length(r)) return(data.frame(mz = numeric(0), ann = character(0), stringsAsFactors = FALSE))
  data.frame(mz = as.numeric(vapply(r, `[`, "", 1)), ann = vapply(r, `[`, "", 2), stringsAsFactors = FALSE)
}

# chain-specific products for one molecular species (vector of chain tokens)
chain_products <- function(lip, MH, toks) {
  b <- lip$base; res <- list()
  add <- function(mz, ann, chain) res[[length(res) + 1]] <<- data.frame(mz = mz, ann = ann, chain = chain, stringsAsFactors = FALSE)
  for (i in seq_along(toks)) {
    t <- toks[i]; x <- tok_parse(t)
    if (b == "CL") { add(dg_like(x$c, x$d), paste0("DG-like ion ", t), t); next }
    if (b %in% c("Cer", "HexCer", "Hex2Cer")) {
      if (i == 1 && x$ox >= 2) add(lcb_ion(x$c, x$d, x$ox), paste0("LCB ", t, " [-2H2O]"), t)
      next
    }
    if (x$ether) {
      if (b == "PE" && x$d >= 1) add(plas_sn1(x$c, x$d), paste0("plasmalogen sn-1 ion P-", x$c, ":", x$d - 1), t)
      next
    }
    f <- fa_mass(x$c, x$d)
    if (b == "PC") {
      add(MH - f, paste0("[M+H-FA ", t, "]+"), t)
      add(MH - f + H2O, paste0("[M+H-ketene ", t, "]+"), t)
    } else if (b %in% c("PE", "PS", "PI", "PG", "PA")) {
      add(acylium(x$c, x$d) + C3H6O2, paste0("acyl+C3H6O2 ", t), t)
      add(acylium(x$c, x$d), paste0("acylium ", t), t)
    } else if (b %in% c("TG", "DG")) {
      add(MH - f, paste0("[M+H-FA ", t, "]+ (loss of FA", if (TRUE) "+NH3" else "", ")"), t)
    }
  }
  if (!length(res)) return(data.frame(mz = numeric(0), ann = character(0), chain = character(0), stringsAsFactors = FALSE))
  unique(do.call(rbind, res))
}

# ---- mining molecular species from MS2 --------------------------------------
combos <- function(L, C, D, n, ether = FALSE, ether_maxd = 3) {
  out <- list(); if (!nrow(L)) return(out)
  ok_ether <- function(c, d) c >= 12 && c <= 24 && d >= 0 && d <= ether_maxd
  inL <- function(c, d) { i <- which(L$c == c & L$d == d); if (length(i)) i[1] else NA }
  if (n == 1 && ether) {
    for (i in seq_len(nrow(L))) {
      rc <- C - L$c[i]; rd <- D - L$d[i]
      if (ok_ether(rc, rd)) out[[length(out) + 1]] <- list(toks = c(paste0("O-", rc, ":", rd), L$tok[i]), score = L$int[i])
    }
  } else if (n == 2) {
    for (i in seq_len(nrow(L))) for (j in i:nrow(L)) {
      rc <- C - L$c[i] - L$c[j]; rd <- D - L$d[i] - L$d[j]
      if (!ether) { if (rc == 0 && rd == 0) out[[length(out) + 1]] <- list(toks = c(L$tok[i], L$tok[j]), score = min(L$int[c(i, j)])) }
      else if (ok_ether(rc, rd)) out[[length(out) + 1]] <- list(toks = c(paste0("O-", rc, ":", rd), L$tok[i], L$tok[j]), score = min(L$int[c(i, j)]))
    }
  } else if (n == 3) {
    for (i in seq_len(nrow(L))) for (j in i:nrow(L)) {
      k <- inL(C - L$c[i] - L$c[j], D - L$d[i] - L$d[j])
      if (!is.na(k) && k >= j) out[[length(out) + 1]] <- list(toks = c(L$tok[i], L$tok[j], L$tok[k]), score = min(L$int[c(i, j, k)]))
    }
  }
  out
}

evidence_list <- function(int, relthr, absthr) {
  L <- FA_TAB; L$int <- int
  if (!any(int > 0)) return(L[0, ])
  L <- L[L$int > 0 & L$int >= relthr * max(int) & L$int >= absthr, , drop = FALSE]
  L[order(L$c, L$d), ]
}

mine_pos <- function(lip, MH, spec, o) {
  if (!lip$ok || lip$labelled || is.na(MH) || !nrow(spec)) return(list())
  b <- lip$base; C <- lip$C; D <- lip$D; tol <- o$tol
  ev <- function(int) evidence_list(int, o$relthr, o$absthr)
  res <- list()
  if (b == "PC") {
    L <- ev(pmax(pk(spec, MH - FA_TAB$m, tol), pk(spec, MH - FA_TAB$m + H2O, tol)))
    res <- if (lip$ether) combos(L, C, D, 1, TRUE) else combos(L, C, D, 2)
  } else if (b %in% c("PE", "PS", "PI", "PG", "PA")) {
    L <- ev(pmax(pk(spec, acylium(FA_TAB$c, FA_TAB$d), tol), pk(spec, acylium(FA_TAB$c, FA_TAB$d) + C3H6O2, tol)))
    res <- if (lip$ether) combos(L, C, D, 1, TRUE) else combos(L, C, D, 2)
    if (lip$ether && b == "PE" && length(res)) {     # plasmalogen: require / reward the sn-1 ion
      for (i in seq_along(res)) {
        x <- tok_parse(res[[i]]$toks[1])
        s1 <- if (x$d >= 1) pk(spec, plas_sn1(x$c, x$d), tol) else 0
        res[[i]]$score <- res[[i]]$score + s1; res[[i]]$sn1 <- s1
      }
      if (any(vapply(res, function(r) r$sn1 > 0, TRUE))) res <- res[vapply(res, function(r) r$sn1 > 0, TRUE)]
    }
  } else if (b %in% c("TG", "DG")) {
    L <- ev(pk(spec, MH - FA_TAB$m, tol))
    n <- if (b == "TG") 3 else 2
    res <- if (lip$ether) combos(L, C, D, n - 1, TRUE, 2) else combos(L, C, D, n)
  } else if (b == "CL") {
    U <- expand.grid(c = 24:48, d = 0:10); U$int <- pk(spec, dg_like(U$c, U$d), tol)
    U <- U[U$int > 0 & U$int >= o$relthr * max(U$int, 1), , drop = FALSE]
    for (i in seq_len(nrow(U))) {
      j <- which(U$c == C - U$c[i] & U$d == D - U$d[i])
      if (length(j) && (U$c[i] < U$c[j] || (U$c[i] == U$c[j] && U$d[i] <= U$d[j])))
        res[[length(res) + 1]] <- list(toks = c(paste0(U$c[i], ":", U$d[i]), paste0(U$c[j], ":", U$d[j])), score = min(U$int[c(i, j)]))
    }
  } else if (b %in% c("Cer", "HexCer") && lip$ox == 2) {
    B <- data.frame(c = c(16, 17, 18, 18, 18, 19, 20), d = c(1, 1, 0, 1, 2, 1, 1))
    B$int <- pk(spec, lcb_ion(B$c, B$d, 2), tol)
    B <- B[B$int > 0 & B$int >= o$relthr * max(B$int, 1), , drop = FALSE]
    for (i in seq_len(nrow(B))) {
      rc <- C - B$c[i]; rd <- D - B$d[i]
      if (rc >= 12 && rc <= 28 && rd >= 0 && rd <= 2)
        res[[length(res) + 1]] <- list(toks = c(paste0(B$c[i], ":", B$d[i], ";O2"), paste0(rc, ":", rd)), score = B$int[i])
    }
  }
  res
}

mine_neg <- function(lip, spec, o) {
  if (!lip$ok || lip$labelled || !nrow(spec)) return(list())
  if (!(lip$base %in% c("PC", "PE", "PS", "PI", "PG", "PA"))) return(list())
  L <- evidence_list(pk(spec, FA_TAB$m - PROTON, o$tol), o$relthr, o$absthr)
  if (lip$ether) combos(L, lip$C, lip$D, 1, TRUE) else combos(L, lip$C, lip$D, 2)
}

# ---- feature preparation ----------------------------------------------------
prep <- function(df, o, positive = TRUE) {
  num <- function(x) suppressWarnings(as.numeric(x))
  # mean intensity over the samples in use (default: file type 'Sample'; o$exclude overrides)
  sm <- attr(df, "samples"); use <- integer(0)
  if (!is.null(sm) && nrow(sm)) use <- if (is.null(o$exclude)) sm$idx[sm$type == "Sample"] else sm$idx[!(sm$name %in% o$exclude)]
  df[[".abund"]] <- if (length(use)) rowMeans(vapply(df[use], function(x) num(as.character(x)), numeric(nrow(df))), na.rm = TRUE) else NA_real_
  nm <- as.character(df[["Metabolite name"]])
  keep <- !is.na(nm) & nzchar(nm) & !grepl("^(Unknown|null)$|^RIKEN|no MS2", nm, ignore.case = TRUE)
  if (!isTRUE(o$keep_wo_ms2)) keep <- keep & !grepl("^(w/o MS2:|low score:)", nm)
  df <- df[keep, , drop = FALSE]
  ref <- if ("Reference m/z" %in% names(df)) num(df[["Reference m/z"]]) else rep(NA_real_, nrow(df))
  avg <- num(df[["Average Mz"]])
  mz <- if (isTRUE(o$use_ref_mz)) ifelse(is.na(ref) | abs(ref - avg) > 0.05, avg, ref) else avg
  sn <- if ("S/N average" %in% names(df)) num(df[["S/N average"]]) else rep(NA_real_, nrow(df))
  lapply(seq_len(nrow(df)), function(i) {
    lip <- parse_lipid(df[["Metabolite name"]][i])
    ad <- as.character(df[["Adduct type"]][i])
    MH <- if (positive && ad %in% names(ADDUCT)) mz[i] - ADDUCT[[ad]] + PROTON else NA_real_
    list(id = as.character(df[[1]][i]), name = df[["Metabolite name"]][i], adduct = ad,
         rt = num(df[["Average Rt(min)"]][i]), mz = mz[i], MH = MH, sn = sn[i],
         abund = { a <- df[[".abund"]][i]; if (is.null(a) || !is.finite(a)) 0 else a },
         spec = parse_spec(as.character(df[["MS/MS spectrum"]][i])), lip = lip)
  })
}

add_species <- function(tab, toks, src, score = NA) {
  k <- paste(toks, collapse = "_")
  if (k %in% names(tab)) { tab[[k]]$src <- union(tab[[k]]$src, src); tab[[k]]$score <- max(tab[[k]]$score, score, na.rm = TRUE) }
  else tab[[k]] <- list(toks = toks, src = src, score = if (is.na(score)) 0 else score)
  tab
}

# ---- main engine ------------------------------------------------------------
build_method <- function(pos_df, neg_df = NULL, o) {
  P <- prep(pos_df, o, TRUE)
  N <- if (!is.null(neg_df)) prep(neg_df, o, FALSE) else list()
  if (!length(P)) stop("No annotated features found in the positive file.")
  n_low <- 0
  if (isTRUE(o$min_abund > 0)) {     # features too weak in the samples in use are left out (internal standards always kept)
    ok <- vapply(P, function(f) f$lip$istd || f$abund >= o$min_abund, TRUE); n_low <- sum(!ok); P <- P[ok]
    if (!length(P)) stop("No features reach the minimum abundance in the samples in use.")
  }
  nkey <- vapply(N, function(f) as.character(f$lip$key), ""); nrt <- vapply(N, function(f) f$rt, 0)
  neg_used <- rep(FALSE, length(N))

  # 1. molecular species per positive feature
  for (i in seq_along(P)) {
    f <- P[[i]]; lip <- f$lip; tab <- list(); negmatch <- character(0)
    pc_only <- isTRUE(o$pc184) && lip$ok && lip$base == "PC"
    if (lip$ok && !lip$labelled) {
      for (s in lip$species) tab <- add_species(tab, s, "MS-DIAL pos")
      if (isTRUE(o$mine_pos) && !pc_only) for (r in mine_pos(lip, f$MH, f$spec, o)) tab <- add_species(tab, species_key(lapply(r$toks, tok_parse), lip$base), "pos MS2", r$score)
      j <- which(nkey == lip$key & abs(nrt - f$rt) <= o$rt_tol)
      for (jj in j) {
        g <- N[[jj]]; neg_used[jj] <- TRUE; negmatch <- c(negmatch, g$name)
        if (lip$base != "CL") for (s in g$lip$species) tab <- add_species(tab, s, "MS-DIAL neg")
        if (isTRUE(o$mine_neg)) for (r in mine_neg(g$lip, g$spec, o)) tab <- add_species(tab, species_key(lapply(r$toks, tok_parse), lip$base), "neg MS2", r$score)
      }
    }
    if (pc_only && length(tab)) {      # PC identities come from negative mode whenever a negative match exists
      fromneg <- vapply(tab, function(s) any(grepl("neg", s$src)), TRUE)
      if (any(fromneg)) tab <- tab[fromneg]
    }
    if (length(tab) > o$max_species) {
      named <- vapply(tab, function(s) any(grepl("MS-DIAL", s$src)), TRUE)
      ord <- order(-named, -vapply(tab, function(s) length(s$src), 1), -vapply(tab, function(s) s$score, 1))
      tab <- tab[ord][seq_len(o$max_species)]
    }
    # 2. candidate products with observed intensities
    base_int <- { s <- f$spec; s <- s[s[, 1] < f$mz - 1.5, , drop = FALSE]; if (nrow(s)) max(s[, 2]) else 0 }
    cand <- data.frame(mz = numeric(0), ann = character(0), type = character(0), species = character(0), chain = character(0), stringsAsFactors = FALSE)
    if (lip$ok && !is.na(f$MH)) {
      cp <- class_products(lip, f$MH)
      if (nrow(cp)) cand <- rbind(cand, data.frame(cp, type = "class", species = NA_character_, chain = NA_character_, stringsAsFactors = FALSE))
      if (!pc_only) for (k in names(tab)) {
        ch <- chain_products(lip, f$MH, tab[[k]]$toks)
        if (nrow(ch)) cand <- rbind(cand, data.frame(mz = ch$mz, ann = ch$ann, type = "chain", species = k, chain = ch$chain, stringsAsFactors = FALSE))
      }
    }
    cand <- cand[cand$mz > 50 & cand$mz < f$mz - 1.5, , drop = FALSE]
    cand$int <- pk(f$spec, cand$mz, o$tol)
    P[[i]]$pc_only <- pc_only; P[[i]]$species <- tab; P[[i]]$cand <- cand; P[[i]]$base_int <- base_int; P[[i]]$negmatch <- unique(negmatch)
  }

  # 3. isobar / isomer neighbours: same Q1 window and not chromatographically resolved
  mzv <- vapply(P, function(f) f$mz, 0); rtv <- vapply(P, function(f) f$rt, 0)
  nb <- lapply(seq_along(P), function(i) setdiff(which(abs(mzv - mzv[i]) <= o$q1 / 2 & abs(rtv - rtv[i]) <= o$coelute), i))
  iso2 <- lapply(seq_along(P), function(i) which(abs(mzv + 2 * 1.003355 - mzv[i]) <= o$q1 / 2 & abs(rtv - rtv[i]) <= o$coelute))
  conflicts <- function(mz, i) {
    hit <- integer(0)
    for (j in nb[[i]]) {
      g <- P[[j]]
      # theory decides when the neighbour has fragment rules (co-isolated DDA spectra are chimeric);
      # its observed MS2 peaks are only used when it has no rules
      bad <- if (nrow(g$cand)) any(abs(g$cand$mz - mz) <= o$q3 / 2)
             else nrow(g$spec) && any(abs(g$spec[, 1] - mz) <= o$q3 / 2 & g$spec[, 2] >= 0.05 * max(g$base_int, 1))
      if (bad) hit <- c(hit, j)
    }
    hit
  }

  # 4. transition selection
  TR <- list(); TG <- list()
  for (i in seq_along(P)) {
    f <- P[[i]]; lip <- f$lip; cand <- f$cand; sp <- f$species
    plas <- lip$cls == "PE-O"
    spname <- function(k) species_label(lip$base, sp[[k]]$toks, plas)
    sumname <- if (lip$ok) lip$sumname else trimws(gsub("\\[IS\\]", "", f$name))
    cname <- if (length(sp) == 1) spname(names(sp)[1]) else sumname
    if (lip$labelled) cname <- trimws(gsub("\\[IS\\]|\\|.*$", "", f$name))
    if (lip$istd) cname <- paste(cname, "[IS]")
    rows <- list(); notes <- character(0); main_name <- cname
    put <- function(r, role, compound, flag = "") rows[[length(rows) + 1]] <<- data.frame(
      compound = compound, product = r$mz, ann = r$ann, type = r$type, role = role, int = r$int,
      flag = flag, stringsAsFactors = FALSE)
    obs <- cand[cand$int > 0, , drop = FALSE]; obs <- obs[order(-obs$int), , drop = FALSE]
    co <- obs[obs$type == "class", , drop = FALSE]; ch <- obs[obs$type == "chain", , drop = FALSE]
    quant <- NULL
    if (nrow(co) && (!nrow(ch) || co$int[1] >= o$class_pref * ch$int[1])) quant <- co[1, ]
    else if (nrow(ch)) quant <- ch[1, ]
    else if (nrow(co)) quant <- co[1, ]
    qflag <- ""
    if (is.null(quant)) {            # nothing rule-based observed -> fall back
      s <- f$spec; s <- s[s[, 1] < f$mz - 1.5 & s[, 1] > 50, , drop = FALSE]
      cl <- cand[cand$type == "class", , drop = FALSE]
      if (nrow(cl) && (!nrow(s) || cl$mz[1] < min(f$spec[, 1]) - 1)) { quant <- cl[1, ]; qflag <- "predicted - below MS2 scan range or no MS2" }
      else if (nrow(s)) {
        top <- s[which.max(s[, 2]), ]
        if (nrow(cl) && top[2] < 3 * o$absthr) { quant <- cl[1, ]; qflag <- "predicted - not seen in MS2" }
        else { quant <- data.frame(mz = top[1], ann = "most intense MS2 fragment (unassigned)", type = "ms2", species = NA_character_, chain = NA_character_, int = top[2], stringsAsFactors = FALSE)
               qflag <- "from MS2, no class rule" }
      }
    }
    status <- "OK"
    if (is.null(quant)) { status <- "no usable fragment"; }
    else {
      put(quant, "Quantifier", cname, qflag)
      used <- quant$mz; n_main <- 1
      # 4a. co-eluting isobars in the target list
      if (length(nb[[i]])) {
        hit <- conflicts(quant$mz, i)
        if (length(hit)) {
          alt <- obs[abs(obs$mz - quant$mz) > o$q3 / 2, , drop = FALSE]; found <- FALSE
          for (a in seq_len(nrow(alt))) if (!length(conflicts(alt$mz[a], i))) {
            put(alt[a, ], "Qualifier (discriminating)", cname); used <- c(used, alt$mz[a]); n_main <- n_main + 1; found <- TRUE; break
          }
          who <- paste(unique(vapply(P[hit], function(g) g$lip$sumname, "")), collapse = "; ")
          if (found) { status <- "isobar resolved by 2nd transition"; notes <- c(notes, paste("quantifier shared with", who)) }
          else { status <- "RT separation only"
                 notes <- c(notes, sprintf("no unique fragment vs %s (closest %.2f min away)", who, min(abs(rtv[hit] - f$rt)))) }
        } else notes <- c(notes, "co-eluting isobar present, quantifier is unique")
      }
      # 4b. isomeric molecular species inside this feature
      if (length(sp) > 1 && isTRUE(o$split_isomers) && !f$pc_only) {
        status <- if (status == "OK") "isomers split" else status
        alltok <- lapply(sp, function(s) unique(s$toks))
        for (k in names(sp)) {
          others <- unique(unlist(alltok[names(alltok) != k]))
          c1 <- cand[cand$type == "chain" & !is.na(cand$species) & cand$species == k, , drop = FALSE]
          uq <- c1[!(c1$chain %in% others), , drop = FALSE]
          pool <- if (nrow(uq)) uq else c1
          fl <- if (nrow(uq)) "" else "chain shared with another isomer"
          if (!nrow(pool)) { notes <- c(notes, paste("no chain fragment rule for", spname(k))); next }
          pool <- pool[order(-pool$int), , drop = FALSE]
          if (pool$int[1] == 0) { if (!isTRUE(o$predicted)) { notes <- c(notes, paste(spname(k), "chain ion not seen in pos MS2")); next }
                                  fl <- trimws(paste(fl, "predicted - not seen in MS2")) }
          if (any(abs(used - pool$mz[1]) < 1e-4)) {    # already monitored: label the quantifier with the species
            if (abs(used[1] - pool$mz[1]) < 1e-4) { rows[[1]]$compound <- paste(spname(k), if (lip$istd) "[IS]" else ""); main_name <- trimws(rows[[1]]$compound) }
            next }
          put(pool[1, ], "Isomer-specific", paste0(spname(k), if (lip$istd) " [IS]" else ""), fl); used <- c(used, pool$mz[1])
        }
      }
      # 4c. user-requested minimum number of transitions
      if (n_main < o$min_tr) {
        alt <- obs[!vapply(obs$mz, function(m) any(abs(used - m) <= o$q3 / 2), TRUE), , drop = FALSE]
        for (a in seq_len(min(nrow(alt), o$min_tr - n_main))) { put(alt[a, ], "Qualifier", main_name); used <- c(used, alt$mz[a]) }
      }
    }
    if (length(iso2[[i]])) notes <- c(notes, paste("M+2 isotope overlap from", paste(unique(vapply(P[iso2[[i]]], function(g) g$lip$sumname, "")), collapse = "; ")))
    if (length(rows)) {
      r <- do.call(rbind, rows); r$compound <- trimws(r$compound)
      r$compound[r$role == "Qualifier (discriminating)"] <- main_name
      r <- r[seq_len(min(nrow(r), o$max_tr)), , drop = FALSE]
      r$fid <- i; r$class <- lip$cls; r$precursor <- f$mz; r$adduct <- f$adduct; r$rt <- f$rt
      r$abund <- f$abund
      r$identities <- if (length(sp)) paste(vapply(names(sp), spname, ""), collapse = "; ") else ""
      r$istd <- lip$istd; r$rel <- if (f$base_int > 0) round(100 * r$int / f$base_int, 1) else NA
      TR[[length(TR) + 1]] <- r
    }
    TG[[length(TG) + 1]] <- data.frame(
      ID = f$id, `MS-DIAL name` = f$name, Class = lip$cls, Adduct = f$adduct, Precursor = round(f$mz, 4), RT = round(f$rt, 3),
      Abundance = signif(f$abund, 4),
      `Molecular species` = if (length(sp)) paste(vapply(names(sp), function(k) paste0(spname(k), " [", paste(sp[[k]]$src, collapse = ", "), "]"), ""), collapse = "; ") else "",
      `Neg match` = paste(f$negmatch, collapse = "; "),
      `Co-eluting isobars` = paste(vapply(P[nb[[i]]], function(g) g$name, ""), collapse = "; "),
      Transitions = if (length(rows)) min(length(rows), o$max_tr) else 0L, Status = status, Notes = paste(notes, collapse = " | "),
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  tr <- do.call(rbind, TR); tg <- do.call(rbind, TG)
  # unique compound names (same lipid annotated at >1 RT)
  ck <- unique(tr[, c("fid", "compound", "rt")]); dup <- ck$compound %in% ck$compound[duplicated(ck$compound)]
  ck$newname <- ifelse(dup, sprintf("%s RT%.2f", ck$compound, ck$rt), ck$compound)
  ck$newname <- make.unique(ck$newname, sep = "_")
  tr$compound <- ck$newname[match(paste(tr$fid, tr$compound), paste(ck$fid, ck$compound))]
  neg_only <- if (length(N)) data.frame(Name = vapply(N[!neg_used], function(f) f$name, ""), Adduct = vapply(N[!neg_used], function(f) f$adduct, ""),
                                         RT = vapply(N[!neg_used], function(f) f$rt, 0), stringsAsFactors = FALSE) else NULL
  list(transitions = tr, targets = tg, neg_only = neg_only, n_low = n_low, n_neg = length(N), n_neg_used = sum(neg_used))
}

# ---- dMRM timing ------------------------------------------------------------
dmrm_timing <- function(rt, delta, cycle_ms, isd_ms) {
  st <- rt - delta / 2; en <- rt + delta / 2
  grid <- seq(max(0, min(st)), max(en), by = 0.01)
  n <- findInterval(grid, sort(st)) - findInterval(grid - 1e-9, sort(en))
  worst <- vapply(seq_along(rt), function(i) max(n[grid >= st[i] & grid <= en[i]], 1), 0)
  list(grid = grid, n = n, dwell = ifelse(n > 0, cycle_ms / n - isd_ms, NA), worst = worst, min_dwell = cycle_ms / worst - isd_ms)
}

# Drop the least abundant TG transitions wherever more MRMs are concurrent than the cycle can hold.
# Priority = mean sample intensity of the feature (MS-DIAL), scaled by relative fragment intensity
# for isomer-specific / qualifier rows, so minor isomers of a TG go before its main transition.
trim_tg <- function(tr, delta, cap, classes = c("TG", "TG-O")) {
  st <- tr$rt_method - delta / 2; en <- tr$rt_method + delta / 2
  grid <- seq(max(0, min(st)), max(en), by = 0.01)
  cover <- lapply(seq_len(nrow(tr)), function(i) which(grid >= st[i] & grid <= en[i]))
  n <- integer(length(grid)); for (i in seq_along(cover)) n[cover[[i]]] <- n[cover[[i]]] + 1L
  fmax <- ave(tr$int, tr$fid, FUN = max)
  prio <- ifelse(tr$role == "Quantifier", tr$abund, tr$abund * ifelse(fmax > 0, tr$int / fmax, 0))
  cand <- tr$class %in% classes & !tr$istd
  keep <- rep(TRUE, nrow(tr)); blocked <- rep(FALSE, length(grid))
  repeat {
    over <- which(n > cap & !blocked); if (!length(over)) break
    g <- over[which.max(n[over])]
    idx <- which(keep & cand & st <= grid[g] & en >= grid[g])
    if (!length(idx)) { blocked[g] <- TRUE; next }
    d <- idx[which.min(prio[idx])]
    if (tr$role[d] == "Quantifier") d <- which(keep & tr$fid == tr$fid[d] & tr$compound == tr$compound[d])
    for (j in d) { keep[j] <- FALSE; n[cover[[j]]] <- n[cover[[j]]] - 1L }
  }
  list(keep = keep, prio = prio)
}

finalise <- function(res, ce, o) {
  tr <- res$transitions
  tr$rt_method <- round(tr$rt * o$rt_slope + o$rt_offset, 3)
  cap <- max(1, floor(o$cycle_ms / (o$min_dwell + o$isd_ms)))
  tr$priority <- NA_real_; tr$concurrent_before <- NA_integer_
  dropped <- cbind(tr[0, ], kept_n = integer(0), tg_rank = character(0), reason = character(0))
  if (isTRUE(o$trim_tg)) {
    tt <- trim_tg(tr, o$delta_rt, cap); tr$priority <- tt$prio
    tr$concurrent_before <- dmrm_timing(tr$rt_method, o$delta_rt, o$cycle_ms, o$isd_ms)$worst
    istg <- tr$class %in% c("TG", "TG-O") & !tr$istd
    d <- tr[!tt$keep, , drop = FALSE]; k <- tr[tt$keep, , drop = FALSE]
    if (nrow(d)) {     # say why each one went
      fa <- tapply(tr$abund[istg], tr$fid[istg], max); rk <- rank(-fa, ties.method = "min")
      d$kept_n <- vapply(d$fid, function(x) sum(k$fid == x), 1L)
      d$tg_rank <- sprintf("%d of %d", as.integer(rk[as.character(d$fid)]), length(fa))
      d$reason <- ifelse(d$kept_n > 0,
        sprintf("Enough other ions: TG still monitored by %d transition%s; this was its lowest-intensity one", d$kept_n, ifelse(d$kept_n == 1, "", "s")),
        "Low abundance: whole TG removed (no transitions left)")
      dropped <- d
    }
    tr <- k
  }
  cls <- ifelse(tr$class %in% ce$Class, tr$class, ifelse(sub("-O$", "", tr$class) %in% ce$Class, sub("-O$", "", tr$class), "other"))
  i <- match(cls, ce$Class)
  tr$CE <- ifelse(tr$type == "class", ce$CE_headgroup[i], ce$CE_chain[i])
  tm <- dmrm_timing(tr$rt_method, o$delta_rt, o$cycle_ms, o$isd_ms)
  tr$concurrent <- tm$worst; tr$min_dwell <- round(tm$min_dwell, 1)
  list(tr = tr, timing = tm, dropped = dropped, cap = cap)
}

agilent_table <- function(tr, o) data.frame(
  `Compound Group` = tr$class, `Compound Name` = tr$compound, `ISTD?` = ifelse(tr$istd, "True", "False"),
  `Precursor Ion` = round(tr$precursor, o$digits), `MS1 Res` = o$ms1res,
  `Product Ion` = round(tr$product, o$digits), `MS2 Res` = o$ms2res,
  Primary = ifelse(tr$role == "Quantifier" | tr$role == "Isomer-specific", "True", "False"), Trigger = "False",
  `Ret Time (min)` = tr$rt_method, `Delta Ret Time` = o$delta_rt, Fragmentor = o$fragmentor,
  `Collision Energy` = tr$CE, `Cell Accelerator Voltage` = o$cav, Polarity = "Positive",
  check.names = FALSE, stringsAsFactors = FALSE)

# ---- UI ---------------------------------------------------------------------
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body{font-size:13px} .well{padding:12px} h4{margin-top:14px;font-size:14px;font-weight:600}
    .stat{display:inline-block;margin:0 18px 8px 0} .stat b{font-size:20px;display:block}
    .help-block{font-size:11.5px} .lead-in{color:#555;margin:-6px 0 14px}
    .app-footer{margin:24px 0 16px;padding-top:12px;border-top:1px solid #ddd;color:#666;font-size:12px;line-height:1.7}"))),
  titlePanel("MS-DIAL to Agilent dMRM builder"),
  tags$p(class = "lead-in", "Builds a positive-mode dynamic MRM transition list for Agilent triple quadrupoles from MS-DIAL lipidomics alignment exports. ",
         tags$a(href = REPO_URL, target = "_blank", rel = "noopener", "Source code and documentation on GitHub")),
  sidebarLayout(
    sidebarPanel(width = 3,
      fileInput("pos", "Positive MS-DIAL alignment export (required)", accept = c(".txt", ".tsv", ".csv")),
      fileInput("neg", "Negative MS-DIAL alignment export (identities)", accept = c(".txt", ".tsv", ".csv")),
      h4("Samples"),
      selectizeInput("excl", "Samples excluded (positive file)", choices = NULL, multiple = TRUE),
      helpText("Abundance is the mean over the samples not listed here. It ranks TGs for trimming, sets the simulated peak heights and feeds the filter below. Blanks, QCs and standards are excluded by default."),
      numericInput("min_abund", "Leave out features with mean abundance below", 0, 0, NA, 1000),
      h4("Identification"),
      numericInput("tol", "MS2 fragment tolerance (Da)", 0.01, 0.001, 0.5, 0.005),
      numericInput("rt_tol", "Pos/neg RT match tolerance (min)", 0.15, 0.01, 2, 0.05),
      numericInput("relthr", "Min chain-fragment intensity (% of strongest chain fragment)", 10, 0, 100, 5),
      numericInput("absthr", "Min fragment intensity (absolute)", 0, 0, NA, 1000),
      checkboxInput("mine_neg", "Mine negative MS2 for extra acyl species", TRUE),
      checkboxInput("mine_pos", "Mine positive MS2 for extra acyl species", TRUE),
      numericInput("max_species", "Max molecular species per feature", 4, 1, 12, 1),
      h4("Transition selection"),
      numericInput("q1", "Q1 isolation width (Da)", 0.7, 0.1, 3, 0.1),
      numericInput("q3", "Q3 isolation width (Da)", 0.7, 0.1, 3, 0.1),
      numericInput("coelute", "Treat as co-eluting if RT differs by less than (min)", 0.25, 0.01, 5, 0.05),
      checkboxInput("pc184", "PC / PC-O: m/z 184 only (identities from negative)", TRUE),
      checkboxInput("split_isomers", "Add one transition per isomeric molecular species", TRUE),
      selectizeInput("nosplit", "...but not for these classes", choices = NULL, multiple = TRUE),
      checkboxInput("predicted", "Include predicted chain transitions not seen in pos MS2", TRUE),
      numericInput("min_tr", "Minimum transitions per compound", 1, 1, 5, 1),
      numericInput("max_tr", "Maximum transitions per feature", 6, 1, 20, 1),
      checkboxInput("use_ref_mz", "Use theoretical (Reference) m/z for precursor", TRUE),
      checkboxInput("inc_istd", "Include internal standards", TRUE),
      selectizeInput("classes", "Classes to include (blank = all)", choices = NULL, multiple = TRUE),
      h4("dMRM / instrument"),
      fluidRow(column(6, numericInput("rt_slope", "RT scale", 1, 0.1, 10, 0.01)),
               column(6, numericInput("rt_offset", "RT offset (min)", 0, -30, 30, 0.05))),
      numericInput("delta_rt", "Delta RT - full window (min)", 1.0, 0.1, 10, 0.1),
      numericInput("peak_w", "Peak width (s)", 12.5, 1, 120, 0.5),
      numericInput("cycle_ms", "Cycle time (ms)", 1000, 50, 5000, 50),
      numericInput("isd_ms", "Inter-scan delay per MRM (ms)", 3.5, 0, 20, 0.5),
      numericInput("min_dwell", "Minimum dwell (ms)", 3, 0.5, 100, 0.5),
      checkboxInput("trim_tg", "Drop least abundant TG transitions where a window is over capacity", TRUE),
      fluidRow(column(6, numericInput("fragmentor", "Fragmentor (V)", 380, 0, 600, 10)),
               column(6, numericInput("cav", "Cell acc. (V)", 5, 0, 10, 1))),
      fluidRow(column(6, selectInput("ms1res", "MS1 res", c("Unit", "Wide", "Widest", "Unit/Enh (6490)"))),
               column(6, selectInput("ms2res", "MS2 res", c("Unit", "Wide", "Widest", "Unit/Enh (6490)")))),
      numericInput("digits", "m/z decimals in export", 1, 0, 4, 1),
      downloadButton("dl_agilent", "Agilent dMRM list (.csv)"), br(), br(),
      downloadButton("dl_full", "Annotated transitions (.csv)")
    ),
    mainPanel(width = 9,
      uiOutput("summary"),
      tabsetPanel(
        tabPanel("Transitions", br(), DTOutput("tr_tab")),
        tabPanel("Targets and isomers", br(),
                 helpText("One row per positive feature: which molecular species were found, from which evidence, and whether anything co-elutes in the same Q1 window."),
                 DTOutput("tg_tab")),
        tabPanel("dMRM timing", br(), plotOutput("timing_plot", height = 420), verbatimTextOutput("timing_txt")),
        tabPanel("Simulated chromatogram", br(),
                 fluidRow(column(3, checkboxInput("xic_log", "Log intensity axis", TRUE)),
                          column(3, checkboxInput("xic_dropped", "Show dropped TG (grey, dashed)", FALSE)),
                          column(3, numericInput("xic_labels", "Label the tallest peaks", 10, 0, 50, 1)),
                          column(3, selectizeInput("xic_class", "Classes shown (blank = all)", choices = NULL, multiple = TRUE))),
                 helpText("Each trace is one transition, drawn only inside its dMRM window. Drag to zoom the time axis, double-click to reset, click to list what is being monitored at that time.",
                          "Peak height is the MS-DIAL sample abundance of the feature scaled by the fragment's relative MS2 intensity; peak shape uses the peak width setting with the measured tailing. It is a picture of the method, not a prediction of QQQ response."),
                 plotOutput("xic", height = 480, click = "xic_click", dblclick = "xic_dbl",
                            brush = brushOpts("xic_brush", direction = "x", resetOnNew = TRUE)),
                 h4(textOutput("xic_at")), DTOutput("xic_tab")),
        tabPanel("Dropped TG", br(),
                 helpText("Nothing is dropped unless more MRMs are concurrent than the cycle can hold (capacity = cycle time / (minimum dwell + inter-scan delay)). At the busiest moment the TG transition with the lowest priority score is removed, and this repeats until every moment fits."),
                 helpText("Priority score = the TG's mean abundance in the samples in use, multiplied for isomer-specific ions by that fragment's MS2 intensity relative to the TG's strongest fragment. So a TG's weaker isomer ions go first ('enough other ions'), and a TG disappears completely only when its last transition is the lowest score left in an overloaded window ('low abundance')."),
                 downloadButton("dl_dropped", "Dropped TG list (.csv)"), br(), br(),
                 DTOutput("drop_tab")),
        tabPanel("Collision energies", br(),
                 helpText("Starting values only - double-click a cell to edit. 'Headgroup' applies to class-specific products, 'chain' to acyl-specific products."),
                 DTOutput("ce_tab")),
        tabPanel("Negative-only", br(),
                 helpText("Annotated in negative but with no matching positive feature (same class, sum composition and RT). These are not in the method."),
                 DTOutput("neg_tab")),
        tabPanel("Rules", br(), htmlOutput("rules"))
      )
    )
  ),
  tags$footer(class = "app-footer",
    tags$div(tags$b("MSDIAL-dMRM-builder"), " - source code, documentation and updates: ",
             tags$a(href = REPO_URL, target = "_blank", rel = "noopener", REPO_URL)),
    tags$div("Found a problem or have a request? ",
             tags$a(href = paste0(REPO_URL, "/issues"), target = "_blank", rel = "noopener", "Open an issue on GitHub"), "."),
    tags$div("Copyright 2026 Sarah Hancock. Free software released under the ",
             tags$a(href = "https://www.gnu.org/licenses/gpl-3.0.html", target = "_blank", rel = "noopener", "GNU General Public License v3.0"),
             ", with no warranty. Check every transition list against your own instrument and standards before acquiring data."))
)

# ---- server -----------------------------------------------------------------
server <- function(input, output, session) {
  ce <- reactiveVal(default_ce())
  pos_df <- reactive({ req(input$pos); read_msdial(input$pos$datapath) })
  neg_df <- reactive({ if (is.null(input$neg)) NULL else read_msdial(input$neg$datapath) })

  samples <- reactive(attr(pos_df(), "samples"))
  default_excl <- reactive({ s <- samples(); s$name[s$type != "Sample"] })
  observeEvent(samples(), updateSelectizeInput(session, "excl", choices = samples()$name, selected = default_excl()))
  excl_ready <- reactiveVal("")      # the selector only takes effect once it has caught up with a newly loaded file
  observe({ req(input$pos); if (setequal(as.character(input$excl), default_excl())) excl_ready(input$pos$datapath) })
  excl <- reactive(if (identical(excl_ready(), input$pos$datapath)) as.character(input$excl) else default_excl())

  opts <- reactive(list(
    tol = input$tol, rt_tol = input$rt_tol, relthr = input$relthr / 100, absthr = input$absthr,
    mine_neg = input$mine_neg, mine_pos = input$mine_pos, max_species = input$max_species,
    q1 = input$q1, q3 = input$q3, coelute = input$coelute, split_isomers = input$split_isomers,
    predicted = input$predicted, pc184 = input$pc184, exclude = excl(), min_abund = input$min_abund, min_tr = input$min_tr, max_tr = input$max_tr, class_pref = 0.2,
    use_ref_mz = input$use_ref_mz, keep_wo_ms2 = FALSE))

  built <- reactive({
    o <- opts(); validate(need(all(!vapply(o[names(o) != "exclude"], function(x) length(x) == 0 || is.na(x), TRUE)), "Fill in all settings."))
    validate(need(length(setdiff(samples()$name, o$exclude)) > 0 || nrow(samples()) == 0, "Every sample is excluded - leave at least one in."))
    withProgress(message = "Building transitions...", build_method(pos_df(), neg_df(), o))
  })
  observeEvent(built(), {
    cl <- sort(unique(built()$transitions$class))
    updateSelectizeInput(session, "classes", choices = cl, selected = intersect(isolate(input$classes), cl))
    updateSelectizeInput(session, "xic_class", choices = cl, selected = intersect(isolate(input$xic_class), cl))
    updateSelectizeInput(session, "nosplit", choices = cl, selected = intersect(isolate(input$nosplit), cl))
  })
  final <- reactive({
    res <- built(); tr <- res$transitions
    if (!isTRUE(input$inc_istd)) tr <- tr[!tr$istd, , drop = FALSE]
    if (length(input$classes)) tr <- tr[tr$class %in% input$classes, , drop = FALSE]
    if (length(input$nosplit)) tr <- tr[!(tr$role == "Isomer-specific" & tr$class %in% input$nosplit), , drop = FALSE]
    validate(need(nrow(tr) > 0, "No transitions left after filtering."))
    res$transitions <- tr
    o <- list(rt_slope = input$rt_slope, rt_offset = input$rt_offset, delta_rt = input$delta_rt,
              cycle_ms = input$cycle_ms, isd_ms = input$isd_ms, min_dwell = input$min_dwell, trim_tg = input$trim_tg)
    finalise(res, ce(), o)
  })
  exp_opts <- reactive(list(digits = input$digits, ms1res = sub(" .*", "", input$ms1res), ms2res = sub(" .*", "", input$ms2res),
                            delta_rt = input$delta_rt, fragmentor = input$fragmentor, cav = input$cav))

  output$summary <- renderUI({
    if (is.null(input$pos)) return(div(class = "alert alert-info", "Load a positive-mode MS-DIAL alignment export to begin. Add the negative-mode export to resolve acyl chains."))
    f <- final(); tr <- f$tr; tg <- built()$targets
    stat <- function(v, l) div(class = "stat", tags$b(v), l)
    div(stat(length(unique(tr$fid)), "features"), stat(length(unique(tr$compound)), "compounds"), stat(nrow(tr), "transitions"),
        stat(max(tr$concurrent), "max concurrent"), stat(sprintf("%.1f ms", min(tr$min_dwell)), "min dwell"),
        stat(sum(tg$Status == "RT separation only"), "rely on RT only"), stat(nrow(f$dropped), "TG transitions dropped"),
        stat(sprintf("%d / %d", built()$n_neg_used, built()$n_neg), "neg features matched"),
        stat(sprintf("%d / %d", length(setdiff(samples()$name, excl())), nrow(samples())), "samples used"),
        if (built()$n_low > 0) stat(built()$n_low, "features below min abundance"), hr())
  })

  output$tr_tab <- renderDT({
    tr <- final()$tr
    d <- data.frame(Class = tr$class, Compound = tr$compound, Adduct = tr$adduct, Precursor = round(tr$precursor, 4),
                    Product = round(tr$product, 4), Fragment = tr$ann, Role = tr$role, `MS2 int` = round(tr$int),
                    `% base` = tr$rel, Abundance = signif(tr$abund, 4), RT = tr$rt_method, CE = tr$CE, Concurrent = tr$concurrent, `Min dwell (ms)` = tr$min_dwell,
                    Flag = tr$flag, check.names = FALSE)
    datatable(d, rownames = FALSE, options = list(pageLength = 25, scrollX = TRUE)) |>
      formatStyle("Min dwell (ms)", color = styleInterval(input$min_dwell, c("#b00020", "inherit"))) |>
      formatStyle("Flag", color = "#b06000")
  })
  output$tg_tab <- renderDT(datatable(built()$targets, rownames = FALSE, options = list(pageLength = 25, scrollX = TRUE)) |>
    formatStyle("Status", color = styleEqual("RT separation only", "#b00020"), fontWeight = styleEqual("RT separation only", "bold")))
  # ---- simulated chromatogram ----
  xic_range <- reactiveVal(NULL); xic_time <- reactiveVal(NULL)
  observeEvent(input$xic_brush, { b <- input$xic_brush; if (b$xmax - b$xmin > 0.02) xic_range(c(b$xmin, b$xmax)); session$resetBrush("xic_brush") })
  observeEvent(input$xic_dbl, xic_range(NULL))
  observeEvent(input$xic_click, xic_time(input$xic_click$x))
  xic_data <- reactive({
    f <- final(); cols <- c("compound", "class", "precursor", "product", "ann", "role", "int", "abund", "fid", "rt_method", "istd")
    d <- cbind(f$tr[cols], dropped = FALSE)
    if (isTRUE(input$xic_dropped) && nrow(f$dropped)) d <- rbind(d, cbind(f$dropped[cols], dropped = TRUE))
    if (length(input$xic_class)) d <- d[d$class %in% input$xic_class, , drop = FALSE]
    validate(need(nrow(d) > 0, "No transitions to draw."))
    fmax <- ave(d$int, d$fid, FUN = max)
    rel <- ifelse(d$role == "Quantifier", 1, ifelse(fmax > 0 & d$int > 0, d$int / fmax, 0.02))
    ab <- if (any(d$abund > 0)) pmax(d$abund, min(d$abund[d$abund > 0])) else rep(1, nrow(d))
    d$h <- ab * rel
    d
  })
  output$xic <- renderPlot({
    d <- xic_data(); delta <- input$delta_rt
    s <- input$peak_w / 60 / sqrt(2 * pi); sf <- s * 2 / 2.6; sb <- 1.6 * sf     # tailing: back half 1.6x the front
    xl <- xic_range(); if (is.null(xl)) xl <- range(c(d$rt_method - delta / 2, d$rt_method + delta / 2))
    d <- d[d$rt_method + delta / 2 >= xl[1] & d$rt_method - delta / 2 <= xl[2], , drop = FALSE]
    validate(need(nrow(d) > 0, "No transitions in this time range - double-click to reset the zoom."))
    off <- seq(-delta / 2, delta / 2, length.out = 121)
    shape <- exp(-0.5 * (off / ifelse(off < 0, sf, sb))^2)
    X <- outer(off, d$rt_method, "+"); Y <- outer(shape, d$h)
    logy <- isTRUE(input$xic_log); top <- max(d$h); lo <- top * 1e-5
    if (logy) Y <- pmax(Y, lo)
    cls <- sort(unique(d$class)); pal <- setNames(hcl.colors(max(length(cls), 2), "Dark 3")[seq_along(cls)], cls)
    op <- par(mar = c(4.5, 4.8, 1, 9), xpd = FALSE); on.exit(par(op))
    plot(NA, xlim = xl, ylim = if (logy) c(lo, top * 3) else c(0, top * 1.12), log = if (logy) "y" else "", xaxs = "i",
         xlab = "Retention time (min)", ylab = "Simulated intensity", bty = "l", las = 1, yaxt = if (logy) "n" else "s")
    if (logy) { at <- 10^seq(ceiling(log10(lo)), floor(log10(top * 3))); axis(2, at = at, labels = sprintf("1e%d", round(log10(at))), las = 1) }
    draw <- function(i, ...) if (length(i)) lines(as.vector(rbind(X[, i, drop = FALSE], NA)), as.vector(rbind(Y[, i, drop = FALSE], NA)), ...)
    draw(which(d$dropped), col = "grey60", lty = 2)
    for (k in cls) draw(which(d$class == k & !d$dropped), col = adjustcolor(pal[[k]], 0.75), lwd = 1.2)
    if (!is.null(xic_time())) abline(v = xic_time(), lty = 3)
    n <- min(input$xic_labels, sum(!d$dropped), na.rm = TRUE)
    if (isTRUE(n > 0)) { k <- d[!d$dropped, ]; k <- k[order(-k$h)[seq_len(n)], ]
      text(k$rt_method, k$h, k$compound, cex = 0.72, pos = 3, offset = 0.25, col = pal[k$class]) }
    legend("topright", inset = c(-0.13, 0), legend = cls, col = pal, lwd = 2, bty = "n", cex = 0.85, xpd = NA)
  })
  xic_at <- reactive({
    req(xic_time()); d <- xic_data(); t <- xic_time()
    d <- d[!d$dropped & abs(d$rt_method - t) <= input$delta_rt / 2, , drop = FALSE]; d[order(-d$h), , drop = FALSE]
  })
  output$xic_at <- renderText({
    if (is.null(xic_time())) return("")
    n <- nrow(xic_at()); sprintf("%d transitions monitored at %.2f min - dwell %.1f ms each", n, xic_time(), if (n) input$cycle_ms / n - input$isd_ms else NA)
  })
  output$xic_tab <- renderDT({
    d <- xic_at()
    datatable(data.frame(Class = d$class, Compound = d$compound, Precursor = round(d$precursor, 4), Product = round(d$product, 4), Fragment = d$ann,
                         Role = d$role, RT = d$rt_method, `Simulated height` = signif(d$h, 3), check.names = FALSE),
              rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  })

  drop_table <- reactive({
    d <- final()$dropped
    data.frame(Compound = d$compound, Reason = d$reason, `TG abundance rank` = d$tg_rank, Abundance = signif(d$abund, 4),
               `Priority score` = signif(d$priority, 3), `Transitions kept for this TG` = d$kept_n,
               `Concurrent before trimming` = d$concurrent_before, Capacity = rep(final()$cap, nrow(d)), Role = d$role, Fragment = d$ann,
               `% base` = d$rel, Precursor = round(d$precursor, 4), Product = round(d$product, 4), RT = d$rt_method,
               check.names = FALSE, stringsAsFactors = FALSE)[order(d$rt_method), ]
  })
  output$drop_tab <- renderDT(datatable(drop_table(), rownames = FALSE, options = list(pageLength = 25, scrollX = TRUE)))
  output$dl_dropped <- downloadHandler(filename = function() paste0("dMRM_dropped_TG_", Sys.Date(), ".csv"),
                                       content = function(file) write.csv(drop_table(), file, row.names = FALSE))
  output$neg_tab <- renderDT({ req(built()$neg_only); datatable(built()$neg_only, rownames = FALSE, options = list(pageLength = 25)) })
  output$ce_tab <- renderDT(datatable(ce(), rownames = FALSE, editable = list(target = "cell", disable = list(columns = 0)),
                                      options = list(dom = "t", pageLength = 50, ordering = FALSE)), server = FALSE)
  observeEvent(input$ce_tab_cell_edit, {
    e <- input$ce_tab_cell_edit; x <- ce(); v <- suppressWarnings(as.numeric(e$value))
    if (!is.na(v)) { x[e$row, e$col + 1] <- v; ce(x) }
  })

  output$timing_plot <- renderPlot({
    f <- final(); tm <- f$timing
    op <- par(mar = c(4.5, 4.5, 3, 4.5), mfrow = c(1, 1)); on.exit(par(op))
    plot(tm$grid, tm$n, type = "s", lwd = 2, col = "#1f5fa8", xlab = "Retention time (min)", ylab = "Concurrent MRMs", bty = "l", las = 1)
    par(new = TRUE)
    plot(tm$grid, pmin(tm$dwell, 200), type = "s", col = "#c2571a", axes = FALSE, xlab = "", ylab = "", ylim = c(0, min(200, max(tm$dwell, na.rm = TRUE))))
    axis(4, las = 1, col.axis = "#c2571a"); mtext("Dwell time (ms, capped at 200)", 4, line = 3, col = "#c2571a")
    abline(h = input$min_dwell, lty = 2, col = "#b00020")
    legend("top", c("Concurrent MRMs", "Dwell time", "Dwell limit"), col = c("#1f5fa8", "#c2571a", "#b00020"), lty = c(1, 1, 2), lwd = c(2, 1, 1), bty = "n", horiz = TRUE, inset = -0.07, xpd = TRUE)
  })
  output$timing_txt <- renderText({
    f <- final(); tr <- f$tr; w <- tr[which.max(tr$concurrent), ]; low <- sum(tr$min_dwell < input$min_dwell)
    paste0("Capacity: ", f$cap, " concurrent MRMs (cycle time / (minimum dwell + inter-scan delay))",
           "\nTG transitions dropped to fit: ", nrow(f$dropped), " (", length(setdiff(f$dropped$fid, tr$fid)), " TG features removed entirely)",
           "\nTotal transitions: ", nrow(tr), "\nMax concurrent MRMs: ", max(tr$concurrent), " (around ", w$rt_method, " min)",
           "\nMinimum dwell: ", min(tr$min_dwell), " ms at a ", input$cycle_ms, " ms cycle",
           "\nTransitions below the ", input$min_dwell, " ms dwell limit: ", low,
           if (low) "\n-> still over capacity from non-TG lipids: narrow Delta RT, lengthen the cycle time, or exclude classes." else "",
           "\nPoints across a ", input$peak_w, " s peak: ", round(input$peak_w * 1000 / input$cycle_ms, 1),
           "\n\nDwell = cycle time / concurrent MRMs - inter-scan delay. Points per peak = peak width / cycle time.")
  })

  output$rules <- renderUI(HTML("
  <h4>How transitions are chosen</h4>
  <ol>
  <li><b>Molecular species</b> are collected for every positive feature from four sources: the MS-DIAL positive name, the MS-DIAL negative name
  (same class, sum composition and RT), fatty-acid anions mined from the negative MS2 spectrum, and chain fragments mined from the positive MS2 spectrum.
  A species is only accepted when its chains add up to the sum composition and every acyl chain has a fragment.</li>
  <li><b>Quantifier</b> = the most intense class-specific product seen in the positive MS2 spectrum (head group / neutral loss);
  if the class has none (TG, CL) or it is weak, the most intense chain product is used.</li>
  <li><b>One transition is kept</b> unless (a) another target shares the Q1 window and co-elutes - then a second product that the
  other lipid does not produce is added, or, if no such product exists, the feature is marked 'RT separation only'; or (b) more than one isomeric molecular species was found in the feature -
  then each species gets a transition on a chain that the other isomers do not carry.</li>
  </ol>
  <h4>Positive-mode products used</h4>
  <table class='table table-condensed' style='width:auto'>
  <tr><th>Class</th><th>Head group / class</th><th>Chain-specific</th></tr>
  <tr><td>PC, PC-O, LPC, SM</td><td>m/z 184.07</td><td>none by default; [M+H-FA]+, [M+H-ketene]+ for PC if the 184-only box is unticked</td></tr>
  <tr><td>PE, LPE</td><td>NL 141.02</td><td>acylium, acylium+74</td></tr>
  <tr><td>PE-P</td><td>NL 141.02</td><td>plasmalogen sn-1 ion, sn-2 acylium+74</td></tr>
  <tr><td>PS / PI / PG / PA</td><td>NL 185.01 / 260.03 / 172.01 / 97.98</td><td>acylium, acylium+74</td></tr>
  <tr><td>TG, TG-O, DG</td><td>DG: [M+H-H2O]+</td><td>[M+NH4-NH3-FA]+</td></tr>
  <tr><td>CE</td><td>m/z 369.35</td><td>-</td></tr>
  <tr><td>Cer, HexCer</td><td>[M+H-H2O]+, [M+H-Hex-H2O]+</td><td>long-chain base ion (e.g. 264.27 for 18:1;O2)</td></tr>
  <tr><td>CL</td><td>-</td><td>DG-like ions</td></tr>
  <tr><td>anything else</td><td colspan='2'>most intense MS2 fragment below the precursor (flagged)</td></tr>
  </table>
  <p>Ether chains are written O- except PE, where O-x:1 is shown as P-x:0. RTs come from the positive file; use RT scale/offset if the QQQ gradient differs.</p>"))

  output$dl_agilent <- downloadHandler(
    filename = function() paste0("dMRM_positive_", Sys.Date(), ".csv"),
    content = function(file) write.csv(agilent_table(final()$tr, exp_opts()), file, row.names = FALSE))
  output$dl_full <- downloadHandler(
    filename = function() paste0("dMRM_annotated_", Sys.Date(), ".csv"),
    content = function(file) write.csv(final()$tr, file, row.names = FALSE))
}

shinyApp(ui, server)
