# Just found out that the state/province column in ROI selection and NGS barcoding
# sheets are lost, probably due to max_guess = 10 and the first 10 entries are empty
# I think I can recover most of them by aggregating state info from other sheets
library(googlesheets4)
library(tidyverse)
library(gargle)

client <- gargle_oauth_client_from_json("~/museumst/client_secret_703487854956-6spmmh2nh0c82icoq2dr0r8hv9cia7o9.apps.googleusercontent.com.json")
gs4_auth_configure(client)
gs4_auth()

sheet_url <- "https://docs.google.com/spreadsheets/d/1sJDb9B7AtYmfKv4-m8XR7uc3XXw_k4kGSout8cqZ8bY/edit?usp=sharing"

oth_names <- c("Prequel", "smFISH", "ISS", "De novo",
               "Analysis", "Prequel analysis", "Reanalysis")
oth_sheets <- map_dfr(oth_names, function(x) {
    df <- read_sheet(sheet_url, x)
    df[,c("country", "state/province", "city")] |>
        distinct()
}) |> distinct()

# Which country + city match more than 1 state
oth_sheets <- oth_sheets |>
    add_count(country, city) |>
    arrange(desc(n))

# Oh, so there're two different Columbia cities here, but anyway, it's not in sheets of interest


roi <- read_sheet(sheet_url, "ROI selection")
ngs <- read_sheet(sheet_url, "NGS barcoding")

"Columbia" %in% ngs$city

oth_sheets <- semi_join(oth_sheets, roi, by = "city")
oth_sheets <- semi_join(oth_sheets, ngs, by = "city")

oth_multi <- oth_sheets |>
    filter(n > 1) |>
    arrange(country, city, `state/province`)
oth_first <- oth_multi |>
    group_by(country, city) |>
    slice_head(n = 1)
oth_one <- oth_sheets |> filter(n == 1L)
oth_sheets <- bind_rows(oth_one, oth_first)
oth_sheets$n <- NULL
roi$`state/province` <- NULL
roi <- roi |>
    left_join(oth_sheets, by = c("country", "city"))

ngs$`state/province` <- NULL
ngs <- ngs |>
    left_join(oth_sheets, by = c("country", "city"))

write_sheet(roi, sheet_url, "ROI selection")
write_sheet(ngs, sheet_url, "NGS barcoding")

