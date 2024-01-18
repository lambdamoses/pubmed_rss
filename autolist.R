library(googlesheets4)
library(gmailr)
library(gargle)
library(stringr)
library(lubridate)
library(rss)

# Gmail OAuth-----------
gs4_auth("museumofst@gmail.com", token = secret_read_rds("gs4.rds", "GARGLE_KEY"))
gm_auth("museumofst@gmail.com", token = secret_read_rds("gm.rds", "GARGLE_KEY"))

sheet_url <- "https://docs.google.com/spreadsheets/d/1sJDb9B7AtYmfKv4-m8XR7uc3XXw_k4kGSout8cqZ8bY/edit?usp=sharing"

to_check <- read_sheet(sheet_url, sheet = "to_check")
to_check <- to_check[!is.na(to_check$date_published),]

if (!file.exists("last_checked.rds")) saveRDS(as.POSIXct(0), file = "last_checked.rds")
last_checked <- readRDS("last_checked.rds")

# PubMed--------------
geomx <- "https://pubmed.ncbi.nlm.nih.gov/rss/search/1ZSp7ZOQTWb6f7XIhdYoHOybYnbDfOV6dv96LPTgURfTQ2EgWt/?limit=15&utm_campaign=pubmed-2&fc=20220424112459"
st <- "https://pubmed.ncbi.nlm.nih.gov/rss/search/1HK5U4U_QH8LXfanBvIif8FyuJFOCMRe8gokq-75uxSqekzEmG/?limit=15&utm_campaign=pubmed-2&fc=20220424112150"
visium <- "https://pubmed.ncbi.nlm.nih.gov/rss/search/1R__6bbhMkenq1M5NePyMAVwqIH-27yBjh8XiC1LGlM6ICAAn2/?limit=15&utm_campaign=pubmed-2&fc=20220424112410"
urls <- c(geomx, st, visium)

.get_pubmed_feed <- function(url) {
    feed <- getFeed(url)$items
    df <- data.frame(date_published = vapply(feed, function(x) as.POSIXct(x$pubDate), FUN.VALUE = POSIXct(1)),
                     title = vapply(feed, function(x) x$title, FUN.VALUE = character(1)),
                     pmid = vapply(feed, function(x) x$identifier |> str_remove("^pmid\\:") |> as.integer(), FUN.VALUE = character(1)),
                     journal = vapply(feed, function(x) x$source |> str_to_title(), FUN.VALUE = character(1)),
                     URL = vapply(feed, function(x) {
                         names(x) <- make.unique(names(x))
                         # should be the last identifier
                         name_use <- names(x)[max(which(str_detect(names(x), "^identifier")))]
                         paste0("https://doi.org/", str_remove(x[[name_use]], "^doi\\:"))
                     }, FUN.VALUE = character(1)))
    df <- df[df$date_published > last_checked,]
    if (nrow(df)) {
        df$date_published <- df$date_published |> as.POSIXct() |> format("%Y-%m-%d")
        return(df)
    } else return(NULL)
}

pubmed_res <- lapply(urls, .get_pubmed_feed)

# bio/medRxiv------------------
threads <- gm_threads(num_results = 20)
terms_bio <- c("spatial transcriptomics", "visium", "merfish", "seqfish", "GeoMX", "in situ sequenc", "CosMX", "xenium")

.extract_rxiv_info <- function(m) { # For one message
    entries <- str_split(m, "\\r\\n\\r\\n")[[1]]
    entries <- entries[!str_detect(entries, "Forwarded M|message")]
    entries <- entries[!str_detect(entries, "Unsubscribe")]
    entries <- entries[str_length(entries) > 4L]
    entries[1] <- str_split(entries[1], "Results\\:\\r\\n")[[1]][2]
    entries <- lapply(entries, function(x) str_split(x, "\\r\\n")[[1]])
    # Remove names
    entries <- lapply(entries, function(x) {
    name_inds <- which(str_detect(x, "[A-Z][a-z\\.]+ [A-Z][a-z\\.]+( ?:[A-Z][a-z\\.]+)?((, )|( and ))"))
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
        dmy(str_extract(x, date_regex)) |> format("%Y-%m-%d")
    }, FUN.VALUE = character(1))
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
    title_simp <- .simp_str(new_res$title)
    new_res$string_match <- part_match
    new_res$existing_sheet[part_match %in% irrelevant] <- "irrelevant"

    # Check if it's already present----------
    # What to do? I think I'll note it here since new versions sometimes have new datasets
    # Again only applies to subsequent versions of preprints
    sheets_use <- c("Prequel", "ROI selection", "NGS barcoding", "smFISH", "ISS", "De novo",
                    "Analysis", "Prequel analysis", "Reanalysis")
    for (s in sheets_use) {
        sh <- read_sheet(sheet_url, s)
        sh$date_published <- as.character(sh$date_published)
        ref <- .get_rxiv_compare(sh$URL)
        ref_title <- .simp_str(sh$title)
        inds <- match(part_match, ref)
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
            if (str_detect(new_res$journal[j], "Rxiv$"))
                cols_use <- c("date_published", "title")
            else
                cols_use <- c("date_published", "title", "pmid", "journal", "URL")
            sh[inds_update, cols_use] <- new_res[j, cols_use]
        }
        if (!all(is.na(inds)) && !all(is.na(indst))) {
            write_sheet(sh, sheet_url, sheet = s)
        }
    }

    to_check <- rbind(to_check, new_res)
    to_check <- to_check[!duplicated(to_check$title),]
}

# Write to sheet------------
write_sheet(to_check, sheet_url, sheet = "to_check")
last_checked <- Sys.time()
saveRDS(last_checked, "last_checked.rds")

# Delete tokens
file.remove(list.files(pattern = "museumofst"))
