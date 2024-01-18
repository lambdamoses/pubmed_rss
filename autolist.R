library(googlesheets4)
library(gmailr)
library(gargle)
library(stringr)
library(lubridate)
library(rbiorxiv)
library(rss)

# Gmail OAuth-----------
gs4_auth("museumofst@gmail.com", token = secret_read_rds("gs4.rds", "GARGLE_KEY"))
gm_auth("museumofst@gmail.com", token = secret_read_rds("gm.rds", "GARGLE_KEY"))

sheet_url <- "https://docs.google.com/spreadsheets/d/1sJDb9B7AtYmfKv4-m8XR7uc3XXw_k4kGSout8cqZ8bY/edit?usp=sharing"

to_check <- read_sheet(sheet_url, sheet = "to_check")
to_check <- to_check[!is.na(to_check$date_added),]

if (!file.exists("last_checked.rds")) saveRDS(as.POSIXct(0), file = "last_checked.rds")
last_checked <- readRDS("last_checked.rds")

# PubMed--------------
geomx <- "https://pubmed.ncbi.nlm.nih.gov/rss/search/1ZSp7ZOQTWb6f7XIhdYoHOybYnbDfOV6dv96LPTgURfTQ2EgWt/?limit=15&utm_campaign=pubmed-2&fc=20220424112459"
st <- "https://pubmed.ncbi.nlm.nih.gov/rss/search/1HK5U4U_QH8LXfanBvIif8FyuJFOCMRe8gokq-75uxSqekzEmG/?limit=15&utm_campaign=pubmed-2&fc=20220424112150"
visium <- "https://pubmed.ncbi.nlm.nih.gov/rss/search/1R__6bbhMkenq1M5NePyMAVwqIH-27yBjh8XiC1LGlM6ICAAn2/?limit=15&utm_campaign=pubmed-2&fc=20220424112410"
urls <- c(geomx, st, visium)

.get_pubmed_feed2 <- function(url, to_check) {
    feed <- getFeed(url)$items
    df <- data.frame(date_added = vapply(feed, function(x) as.POSIXct(x$pubDate), FUN.VALUE = POSIXct(1)),
                     title = vapply(feed, function(x) x$title, FUN.VALUE = character(1)),
                     pmid = vapply(feed, function(x) x$identifier |> str_remove("^pmid\\:"), FUN.VALUE = character(1)),
                     journal = vapply(feed, function(x) x$source |> str_to_title(), FUN.VALUE = character(1)),
                     URL = vapply(feed, function(x) x$link, FUN.VALUE = character(1)))
    df <- df[df$date_added > last_checked,]
    if (nrow(df)) {
        df$date_added <- df$date_added |> as.POSIXct() |> format("%Y-%m-%d")
        df$string_match <- NA
        df$existing_sheet <- NA
        to_check <- rbind(to_check, df)
    }
    to_check
}

for (u in urls) {
  to_check <- .get_pubmed_feed(u, to_check)
}

# bio/medRxiv------------------
threads <- gm_threads(num_results = 20)
terms_bio <- c("spatial transcriptomics", "visium", "merfish", "seqfish", "GeoMX", "ISS", "CosMX", "xenium")

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
  data.frame(title = titles,
             pmid = NA,
             URL = urls)
}

# Function to get URL, title, date from one term
.get_rxiv_feed <- function(threads, term, rxiv = "bioRxiv") {
    ids <- gm_id(threads)
    messages <- lapply(ids, function(x) gm_thread(x)$messages[[2]])
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
        df <- data.frame(date = dates,
                         subject = subjects,
                         body = vapply(messages, function(x) gm_body(x)[[1]],
                                       FUN.VALUE = character(1)))
        df <- df[df$date > last_checked,]
        df$date <- vapply(df$date, function(x) format(as.POSIXct(x), "%Y-%m-%d"),
                          FUN.VALUE = character(1))
        if (nrow(df)) {
            cat("Reading", subject, "\n")
            df2 <- mapply(function(body, date) {
                out <- .extract_rxiv_info(body)
                out$date_added <- date
                out$journal <- rxiv
                out
            }, body = df$body, date = df$date, SIMPLIFY = FALSE)
            df2 <- do.call(rbind, df2)
        } else return(NULL)
    } else return(NULL)
    df2
}

biorxiv_res <- lapply(terms_bio, .get_rxiv_feed, threads = threads, rxiv = "bioRxiv")
medrxiv_res <- lapply(terms_bio, .get_rxiv_feed, threads = threads, rxiv = "medRxiv")

rxiv_res <- do.call(rbind, c(biorxiv_res, medrxiv_res))
if (!is.null(rxiv_res)) {
    rxiv_res$existing_sheet <- NA
    rownames(rxiv_res) <- NULL
    rxiv_res <- rxiv_res[!duplicated(rxiv_res$URL),]
}

.get_rxiv_compare <- function(dois, rm_pfx = TRUE) {
    regex_use <- if (rm_pfx) "^(https://doi.org/)?10\\.1101/" else "^https://doi.org/"
    str_remove(dois, regex_use) |>
        str_remove(";$")
}
if (!is.null(rxiv_res)) {
    # Check for relevance---------
    # I'll note it here since sometimes a new version brings spatial stuff into a
    # previous irrelevant paper, but this is rare.
    irrelevant <- read_sheet(sheet_url, sheet = "irrelevant")
    irrelevant <- irrelevant[!is.na(irrelevant$doi),]
    # This only applies to preprints where further versions can show up in RSS
    irrelevant <- irrelevant[str_detect(irrelevant$doi, "^10\\.1101/\\d"),]
    irrelevant <- .get_rxiv_compare(irrelevant$doi)
    part_match <- str_extract(rxiv_res$URL, "(?<=/)[\\d\\.]+$")
    rxiv_res$string_match <- part_match
    rxiv_res$existing_sheet[part_match %in% irrelevant] <- "irrelevant"

    # Check if it's already present----------
    # What to do? I think I'll note it here since new versions sometimes have new datasets
    # Again only applies to subsequent versions of preprints
    sheets_use <- c("Prequel", "ROI selection", "NGS barcoding", "smFISH", "ISS", "De novo",
                    "Analysis", "Prequel analysis")
    for (s in sheets_use) {
        sh <- read_sheet(sheet_url, s)
        sh$date_published <- as.character(sh$date_published)
        ref <- .get_rxiv_compare(sh$URL)
        inds <- match(part_match, ref)
        for (j in seq_along(inds)) {
            i <- inds[j]
            if (is.na(i)) next
            s_use <- paste0(s, " (", i+1, ")") # +1 to account for header
            if (is.na(rxiv_res$existing_sheet[i]))
                rxiv_res$existing_sheet[j] <- s_use
            else
                rxiv_res$existing_sheet[j] <- paste(rxiv_res$existing_sheet[j],
                                                    s_use, sep = ", ")
            # Update titles
            title_old <- sh$title[i]
            title_new <- rxiv_res$title[j]
            inds_update <- ref == part_match
            if (title_new != title_old) {
                sh$title[inds_update] <- title_new
            }
            # Update date, need bioRxiv API
            doi_use <- .get_rxiv_compare(sh$URL[i], FALSE)
            content <- biorxiv_content(server = rxiv_res$journal[j] |> str_to_lower(),
                                       doi = doi_use)
            # Get newest version
            content <- content[[length(content)]]
            sh$date_published[inds_update] <- content$date
        }
        if (!all(is.na(inds))) {
            write_sheet(sh, sheet_url, sheet = s)
        }
    }

    to_check <- rbind(to_check, rxiv_res)
    to_check <- to_check[!duplicated(to_check$URL),]
}

# Write to sheet------------
write_sheet(to_check, sheet_url, sheet = "to_check")
last_checked <- Sys.time()
saveRDS(last_checked, "last_checked.rds")

# Delete tokens
file.remove(list.files(pattern = "museumofst"))
