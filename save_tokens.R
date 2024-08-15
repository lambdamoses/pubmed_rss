library(gargle)
library(gmailr)
library(googlesheets4)
# run once a week since the tokens expire
client <- gargle_oauth_client_from_json("client_secret_703487854956-6spmmh2nh0c82icoq2dr0r8hv9cia7o9.apps.googleusercontent.com.json")

gs4t <- gargle2.0_token(client = client, scope = "https://www.googleapis.com/auth/spreadsheets")
secret_write_rds(gs4t, "gs4.rds", key = "GARGLE_KEY")
gmt <- gargle2.0_token(client = client, scope = "https://www.googleapis.com/auth/gmail.readonly")
secret_write_rds(gmt, "gm.rds", key = "GARGLE_KEY")
