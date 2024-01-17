args <- commandArgs(trailingOnly = TRUE)

library(tidyRSS)
library(googlesheets4)
library(gmailr)
library(gargle)
library(stringr)
library(lubridate)

# Gmail OAuth-----------
client <- gargle_oauth_client(id = args[1], secret = args[2],
                              redirect_uris = "http://localhost", type = "installed",
                              name = "autolist")
gs4_auth_configure(client)
gm_auth_configure(client)
sheet_url <- "https://docs.google.com/spreadsheets/d/1sJDb9B7AtYmfKv4-m8XR7uc3XXw_k4kGSout8cqZ8bY/edit?usp=sharing"

to_check <- read_sheet(sheet_url, sheet = "to_check")
to_check <- to_check[!is.na(to_check$date_published),]

# PubMed--------------
geomx <- "https://pubmed.ncbi.nlm.nih.gov/rss/search/1ZSp7ZOQTWb6f7XIhdYoHOybYnbDfOV6dv96LPTgURfTQ2EgWt/?limit=15&utm_campaign=pubmed-2&fc=20220424112459"
st <- "https://pubmed.ncbi.nlm.nih.gov/rss/search/1HK5U4U_QH8LXfanBvIif8FyuJFOCMRe8gokq-75uxSqekzEmG/?limit=15&utm_campaign=pubmed-2&fc=20220424112150"
visium <- "https://pubmed.ncbi.nlm.nih.gov/rss/search/1R__6bbhMkenq1M5NePyMAVwqIH-27yBjh8XiC1LGlM6ICAAn2/?limit=15&utm_campaign=pubmed-2&fc=20220424112410"
urls <- c(geomx, st, visium)

.get_pubmed_feed <- function(url, to_check) {
  df <- tidyfeed(url)
  df <- df[df$item_pub_date > (Sys.time() - dhours(24)),]

  if (nrow(df)) {
    df$item_guid <- str_remove(df$item_guid, "^pubmed\\:")
    df <- df[, c("item_pub_date", "item_title", "item_guid", "item_link")]
    df$item_pub_date <- format(df$item_pub_date, "%Y-%m-%d")
    df$item_link <- str_remove(df$item_link, "/\\?utm_source.+")
    df <- df[!df$item_link %in% to_check$URL,]
    df$existing_sheet <- NA
    names(df) <- names(to_check)
    to_check <- rbind(to_check, df)
  }
  to_check
}

for (u in urls) {
  to_check <- .get_pubmed_feed(u, to_check)
}

# bio/medRxiv------------------
my_threads <- gm_threads(num_results = 20)
terms_bio <- c("spatial transcriptomics", "visium", "merfish", "seqfish", "GeoMX", "ISS", "CosMX", "xenium")

.extract_rxiv_info <- function(m) { # For one message
  entries <- str_split(m, "\\r\\n\\r\\n")[[1]]
  entries <- entries[!str_detect(entries, "Forwarded M|message")]
  entries <- entries[!str_detect(entries, "Unsubscribe")]
  entries <- entries[str_length(entries) > 1L]
  entries[1] <- str_split(entries[1], "Results\\:\\r\\n")[[1]][2]
  entries <- lapply(entries, function(x) str_split(x, "\\r\\n")[[1]])
  # Remove names
  entries <- lapply(entries, function(x) {
    name_inds <- str_detect(x, "[A-Z][a-z]+ [A-Z][a-z]+( ?:[A-Z][a-z]+)?(, )|( and )")
    x[!name_inds]
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
    messages <- lapply(ids, function(x) gm_thread(x)$messages[[1]])
    # Gmail API doesn't pull body of html emails
    # Have to forward to get the body, hacky but works
    subject <- paste0("Fwd: ", rxiv, ": ", term)
    subjects <- vapply(messages, gm_subject, FUN.VALUE = character(1))
    inds <- which(subjects == subject)
    if (length(inds)) {
        messages <- messages[inds]
        subjects <- subjects[inds]
        dates <- vapply(messages, function(x) as.POSIXct(gm_date(x),
                                                         tryFormats = c("%a, %d %b %Y %T %z",
                                                                        "%a, %e %b %Y %T %z")),
                        FUN.VALUE = Date(1))
        df <- data.frame(date = dates,
                         subject = subjects,
                         body = vapply(messages, function(x) gm_body(x)[[1]],
                                       FUN.VALUE = character(1)))
        df <- df[df$date > (Sys.time() - dhours(24)),]
        if (nrow(df)) {
            df2 <- mapply(function(body, date) {
                out <- .extract_rxiv_info(body)
                out$date_published <- date
                out[,c(4, 1:3)]
            }, body = df$body, date = df$date, SIMPLIFY = FALSE)
            df2 <- do.call(rbind, df2)
        } else return(NULL)
    } else return(NULL)
    df2
}

biorxiv_res <- lapply(terms_bio, .get_rxiv_feed, threads = my_threads, rxiv = "bioRxiv")
medrxiv_res <- lapply(terms_bio, .get_rxiv_feed, threads = my_threads, rxiv = "medRxiv")

rxiv_res <- do.call(rbind, c(biorxiv_res, medrxiv_res))
rxiv_res$existing_sheet <- NA
rxiv_res <- rxiv_res[!duplicated(rxiv_res$URL)]

# Check for relevance---------
# I'll note it here since sometimes a new version brings spatial stuff into a
# previous irrelevant paper, but this is rare.
.get_rxiv_compare <- function(dois) {
    str_remove(dois, "^(https://doi.org/)?10\\.1101/")
}
irrelevant <- read_sheet(sheet_url, sheet = "irrelevant")
irrelevant <- irrelevant[!is.na(irrelevant$doi),]
# This only applies to preprints where further versions can show up in RSS
irrelevant <- irrelevant[str_detect(irrelevant$doi, "^10\\.1101/\\d"),]
irrelevant <- .get_rxiv_compare(irrelevant$doi)
part_match <- str_extract(rxiv_res$URL, "(?<=/)[\\d\\.]+$")
rxiv_res$existing_sheet[part_match %in% irrelevant] <- "irrelevant"

# Check if it's already present----------
# What to do? I think I'll note it here since new versions sometimes have new datasets
# Again only applies to subsequent versions of preprints
sheets_use <- c("Prequel", "ROI selection", "NGS barcoding", "smFISH", "ISS", "De novo",
                "Analysis", "Prequel analysis")
for (s in sheets_use) {
    sh <- read_sheet(sheet_url, s)
    ref <- .get_rxiv_compare(unique(sh$URL))
    inds <- part_match %in% ref
    if (any(inds)) {
        rxiv_res$existing_sheet[inds] <- vapply(rxiv_res$existing_sheet[inds], function(x) {
            if (is.na(x)) s else paste(x, s, sep = ", ")
        }, FUN.VALUE = character(1))
    }
}

to_check <- rbind(to_check, rxiv_res)
to_check <- to_check[!duplicated(to_check$URL)]

# Write to sheet------------
write_sheet(to_check, sheet_url, sheet = "to_check")
