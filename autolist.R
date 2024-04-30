library(googlesheets4)
library(gmailr)
library(gargle)
library(stringr)
library(lubridate)
#library(rss)
library(rbiorxiv)
library(easyPubMed)
library(xml2)
library(httr)
library(tibble)
library(purrr)
# Gmail OAuth-----------
gs4_auth("museumofst@gmail.com", token = secret_read_rds("gs4.rds", "GARGLE_KEY"))
gm_auth("museumofst@gmail.com", token = secret_read_rds("gm.rds", "GARGLE_KEY"))

sheet_url <- "https://docs.google.com/spreadsheets/d/1sJDb9B7AtYmfKv4-m8XR7uc3XXw_k4kGSout8cqZ8bY/edit?usp=sharing"

to_check <- read_sheet(sheet_url, sheet = "to_check")
to_check <- to_check[!is.na(to_check$date_published),]

if (!file.exists("last_checked.rds")) saveRDS(as.POSIXct(0), file = "last_checked.rds")
last_checked <- readRDS("last_checked.rds")
last_checked_pubmed <- readRDS("last_checked_pubmed.rds")
# PubMed--------------
geomx <- "https://pubmed.ncbi.nlm.nih.gov/rss/search/1ZSp7ZOQTWb6f7XIhdYoHOybYnbDfOV6dv96LPTgURfTQ2EgWt/?limit=15&utm_campaign=pubmed-2&fc=20220424112459"
st <- "https://pubmed.ncbi.nlm.nih.gov/rss/search/1HK5U4U_QH8LXfanBvIif8FyuJFOCMRe8gokq-75uxSqekzEmG/?limit=15&utm_campaign=pubmed-2&fc=20220424112150"
visium <- "https://pubmed.ncbi.nlm.nih.gov/rss/search/1R__6bbhMkenq1M5NePyMAVwqIH-27yBjh8XiC1LGlM6ICAAn2/?limit=15&utm_campaign=pubmed-2&fc=20220424112410"
proteomics <- "https://pubmed.ncbi.nlm.nih.gov/rss/search/1dSaW42H3lP27KctEmjamkNUtVEyYRuxY8eM9za16qa1Gb7Qhl/?limit=15&utm_campaign=pubmed-2&fc=20240319195805"
metabolomics <- "https://pubmed.ncbi.nlm.nih.gov/rss/search/1BOJ2J2KKlLHU1JzUwOlCd6YPMI6hUyKdaeuvsgoAK603AJUWN/?limit=15&utm_campaign=pubmed-2&fc=20240319195827"
imc <- "https://pubmed.ncbi.nlm.nih.gov/rss/search/1BOJ2J2KKlLHU1JzUwOlCd6YPMI6hUyKdaeuvsgoAK603AJUWN/?limit=15&utm_campaign=pubmed-2&fc=20240319195827"
msi <- "https://pubmed.ncbi.nlm.nih.gov/rss/search/1BOJ2J2KKlLHU1JzUwOlCd6YPMI6hUyKdaeuvsgoAK603AJUWN/?limit=15&utm_campaign=pubmed-2&fc=20240319195827"

urls <- c(geomx, st, visium, proteomics, metabolomics, imc, msi)

.get_article_type <- function(pmid) {
    info_xml <- get_pubmed_ids(pmid) |>
        fetch_pubmed_data(retmax = 1L) |>
        read_xml()
    xml_find_first(info_xml, "//PublicationTypeList") |> xml_text()
}

.get_feed <- function(url) {
    x <- GET(url) |> read_xml(options = "HUGE")
    res_entry <- xml_find_all(x, "//*[name()='item']") |> as_list()
    res_entry <- map(res_entry, function(x) {
        names(x) <- make.unique(names(x))
        x
    })
    tibble(date_published = map(res_entry, "pubDate") |> unlist(),
           title = map(res_entry, "title") |> unlist(),
           pmid = map(res_entry, ~ .x$identifier |> str_remove("^pmid\\:") |> as.integer()) |> unlist(),
           journal = map(res_entry, ~ .x$source |> str_to_title()) |> unlist(),
           URL = map_chr(res_entry, ~ {
               name_use <- names(.x)[max(which(str_detect(names(.x), "^identifier")))]
               str_replace(.x[[name_use]], "^doi\\:", "https://doi.org/")
           }))
}

.get_pubmed_feed <- function(url) {
    cat("Reading PubMed feed\n")
    df <- .get_feed(url)
    df$date_published <- as.POSIXct(df$date_published,
                                    tryFormats = "%a, %d %b %Y %H:%M:%S %z")
    df <- df[df$date_published > last_checked_pubmed,]
    if (nrow(df)) {
        # Remove protocols
        is_protocol <- str_detect(str_to_lower(df$journal), "(methods in molecular biology)|(jove)|(protocol)")
        df <- df[!is_protocol,]
        df$type <- vapply(df$pmid, .get_article_type, FUN.VALUE = character(1))
        df <- df[!str_detect(df$type, "(R|r)eview"),]
        df$type <- NULL
    }
    if (nrow(df)) return(df) else return(NULL)
}
# When the PubMed RSS is down, internal server error
pubmed_res <- try(lapply(urls, .get_pubmed_feed))
if (is(pubmed_res, "error")) pubmed_res <- NULL else saveRDS(Sys.time(), "last_checked_pubmed.rds")

# bio/medRxiv------------------
threads <- gm_threads(num_results = 20)
terms_bio <- c("spatial transcriptomics", "visium", "merfish", "seqfish",
               "GeoMX", "in situ sequenc", "CosMX", "xenium",
               "imaging mass cytometry", "mass spectrometry imaging",
               "spatial proteomics", "spatial metabolomics")

.extract_rxiv_info <- function(m) { # For one message
    entries <- str_split(m, "\\r\\n\\r\\n")[[1]]
    entries <- entries[!str_detect(entries, "Forwarded (M|m)essage")]
    entries <- entries[!str_detect(entries, "Unsubscribe")]
    entries <- entries[str_length(entries) > 4L]
    entries[1] <- str_split(entries[1], "Results\\:\\r\\n")[[1]][2]
    entries <- lapply(entries, function(x) str_split(x, "\\r\\n")[[1]])
    # Remove names
    entries <- lapply(entries, function(x) {
        # Multiple whole names
        reg1 <- "[A-Z][a-z\\.]*( [A-Z][a-zA-Z\\.-]*?)( [A-Z][a-zA-Z\\.-]*?)?( [A-Z][a-zA-Z\\.-]*?)?((, )|( and ))"
        # Single name
        reg2 <- "^[A-Z][a-z\\.]*( [A-Z][a-zA-Z\\.-]*?)? [A-Z][a-zA-Z\\.-]*( [A-Z][a-zA-Z\\.-]*)?$"
        # Multiple names, the last one incomplete
        reg3 <- "^[A-Z][a-z\\.]*( [A-Z][a-zA-Z\\.-]*?)? [A-Z][a-zA-Z\\.-]*( [A-Z][a-zA-Z\\.-]*)?, [A-Z][a-z\\.]*(( [A-Z][a-zA-Z\\.-]*?)? [A-Z][a-zA-Z\\.-]*( [A-Z][a-zA-Z\\.-]*)?)? $"
        name_inds <- which(str_detect(x, reg1) | str_detect(x, reg2) | str_detect(x, reg3))
        # Anything between first line of names and the last 3 lines are names
        if (length(x)-3 > max(name_inds))
            name_inds <- c(name_inds, seq(max(name_inds) + 1, length(x) - 3))
        x[-name_inds]
  })
    # Extract URL
    urls <- vapply(entries, function(x) {
        x[str_detect(x, "^http\\://")] |> str_trim()
    }, FUN.VALUE = character(1))
    # Extract titles
    titles <- vapply(entries, function(x) {
        x <- x[-c((length(x)-2):length(x))]
        paste(x, collapse = "")
    }, FUN.VALUE = character(1))

    # Extract date posted
    date_regex <- "(?<=posted )\\d+ [A-Za-z]+ 20\\d{2}(?=,)"
    dates <- vapply(entries, function(x) {
        x <- x[str_detect(x, date_regex)]
        as.POSIXct(str_extract(x, date_regex), tryFormats = "%e %B %Y")
    }, FUN.VALUE = POSIXct(1))
    data.frame(date_published = dates,
               title = titles,
               pmid = NA,
               URL = urls)
}

# Function to get URL, title, date from one term
.get_rxiv_feed <- function(threads, term, rxiv = "bioRxiv") {
    ids <- gm_id(threads)
    messages <- lapply(ids, function(x) gm_thread(x)$messages[[length(gm_thread(x)$messages)]])
    # Gmail API doesn't pull body of html emails
    # Have to forward to get the body, hacky but works
    subject <- paste0("Fwd: ", rxiv, ": ", term)
    subjects <- vapply(messages, gm_subject, FUN.VALUE = character(1))
    inds <- which(subjects == subject)
    if (length(inds)) {
        messages <- messages[inds]
        subjects <- subjects[inds]
        dates <- vapply(messages, function(x) {
            as.POSIXct(gm_date(x),
                       tryFormats = c("%a, %d %b %Y %T %z",
                                      "%a, %e %b %Y %T %z"))
        }, FUN.VALUE = POSIXct(1))
        df <- data.frame(date = dates, # Here when I got the feed not when published
                         subject = subjects,
                         body = vapply(messages, function(x) gm_body(x)[[1]],
                                       FUN.VALUE = character(1)))
        df <- df[df$date > last_checked,]
        if (nrow(df)) {
            cat("Reading", subject, "\n")
            df2 <- lapply(df$body, function(body) {
                out <- .extract_rxiv_info(body)
                out$journal <- rxiv
                out
            })
            df2 <- do.call(rbind, df2)
        } else return(NULL)
    } else return(NULL)
    df2
}

biorxiv_res <- lapply(terms_bio, .get_rxiv_feed, threads = threads, rxiv = "bioRxiv")
medrxiv_res <- lapply(terms_bio, .get_rxiv_feed, threads = threads, rxiv = "medRxiv")

new_res <- do.call(rbind, c(biorxiv_res, medrxiv_res, pubmed_res))
if (!is.null(new_res)) {
    new_res$existing_sheet <- NA
    rownames(new_res) <- NULL
    new_res <- new_res[!duplicated(new_res$URL),]
}

.get_rxiv_compare <- function(dois, rm_pfx = TRUE) {
    regex_use <- if (rm_pfx) "^(https://doi.org/)?10\\.1101/" else "^https://doi.org/"
    str_remove(dois, regex_use) |>
        str_remove(";$")
}

.simp_str <- function(x) {
    str_remove_all(x, "\\p{P}") |> str_to_lower()
}

if (!is.null(new_res)) {
    # Check for relevance---------
    # I'll note it here since sometimes a new version brings spatial stuff into a
    # previous irrelevant paper, but this is rare.
    irrelevant <- read_sheet(sheet_url, sheet = "irrelevant")
    irrelevant <- irrelevant[!is.na(irrelevant$doi),]
    # This only applies to preprints where further versions can show up in RSS
    irrelevant <- irrelevant[str_detect(irrelevant$doi, "^10\\.1101/\\d"),]
    irrelevant <- .get_rxiv_compare(irrelevant$doi)
    part_match <- str_extract(new_res$URL, "(?<=/)[\\d\\.]+$")
    part_match[!str_detect(new_res$journal, "Rxiv$")] <- NA
    title_simp <- .simp_str(new_res$title)
    new_res$string_match <- part_match
    new_res$URL <- ifelse(is.na(part_match), new_res$URL, paste0("https://doi.org/10.1101/", part_match))
    new_res$existing_sheet[part_match %in% irrelevant] <- "irrelevant"
    new_res$updated <- FALSE
    # Deal with bioRxiv and medRxiv entries from PubMed
    new_res$journal[str_detect(new_res$journal, "Biorxiv")] <- "bioRxiv"
    new_res$journal[str_detect(new_res$journal, "Medrxiv")] <- "medRxiv"

    # Check if it's already present----------
    # What to do? I think I'll note it here since new versions sometimes have new datasets
    # Again only applies to subsequent versions of preprints
    sheets_use <- c("Prequel", "ROI selection", "NGS barcoding", "smFISH", "ISS", "De novo",
                    "Analysis", "Prequel analysis", "Reanalysis", "multiomics")
    for (s in sheets_use) {
        sh <- read_sheet(sheet_url, s, guess_max = 2000)
        sh$date_published <- as.POSIXct(ymd(sh$date_published), tz = "UTC") |> as.numeric()
        ref <- .get_rxiv_compare(sh$URL)
        ref_title <- .simp_str(sh$title)
        inds <- match(part_match, ref)
        inds[is.na(part_match)] <- NA
        indst <- match(title_simp, ref_title)
        for (j in seq_along(inds)) {
            i <- inds[j]
            if (is.na(i)) i <- indst[j]
            if (is.na(i)) next
            s_use <- paste0(s, " (", i+1, ")") # +1 to account for header
            if (is.na(new_res$existing_sheet[i]))
                new_res$existing_sheet[j] <- s_use
            else
                new_res$existing_sheet[j] <- paste(new_res$existing_sheet[j],
                                                    s_use, sep = ", ")
            # Update entries already found
            inds_update <- which((ref == part_match[j]) | (ref_title == title_simp[j]))
            # Also update URL if new entry is in pubmed
            if (str_detect(new_res$journal[j], "Rxiv$")) {
                # Update date, need bioRxiv API
                doi_use <- .get_rxiv_compare(sh$URL[i], FALSE)
                content <- biorxiv_content(server = new_res$journal[j] |> str_to_lower(),
                                           doi = doi_use)
                # Get newest version
                content <- content[[length(content)]]
                new_res$date_published[j] <- as.POSIXct(content$date, tz = "UTC")
                cols_use <- c("date_published", "title", "journal")
            } else
                cols_use <- c("date_published", "title", "pmid", "journal", "URL")
            if (is.na(sh$date_published[i]) || (new_res$date_published[j] > sh$date_published[i])) {
                sh[inds_update, cols_use] <- new_res[j, cols_use]
                new_res$updated[j] <- TRUE
            }
        }
        if (any(new_res$updated & str_detect(new_res$existing_sheet, s))) {
            # Write sheet if some entries have been updated
            sh$date_published <- format(as.POSIXct(sh$date_published, tz = "UTC"), "%Y/%m/%d")
            write_sheet(sh, sheet_url, sheet = s)
        }
    }
    # Remove old entries that are already updated in the database
    new_res <- new_res[!(!new_res$updated & !is.na(new_res$existing_sheet)), ]
    new_res <- new_res[,c("date_published", "title", "pmid", "journal", "URL", "existing_sheet", "string_match")]
    new_res$updated <- NULL
    new_res$date_published <- format(as.POSIXct(new_res$date_published), "%Y/%m/%d")
    to_check <- rbind(to_check, new_res)
    to_check <- to_check[!duplicated(to_check$title),]
}

# Write to sheet------------
write_sheet(to_check, sheet_url, sheet = "to_check")
last_checked <- Sys.time()
saveRDS(last_checked, "last_checked.rds")
