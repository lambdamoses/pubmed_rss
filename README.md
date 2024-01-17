The code in this repo automatically finds new potentially relevant spatial transcriptomics literature from RSS feeds and adds them to [a Google sheet](https://docs.google.com/spreadsheets/d/1sJDb9B7AtYmfKv4-m8XR7uc3XXw_k4kGSout8cqZ8bY/edit?pli=1#gid=305027376) so multiple curators of the Museum of Spatial Transcriptomics can divide the task of screening and adding the new entries to the database.

## To use the code

1. Make sure to have either `devtools` or `pak` installed
2. Clone this repo
3. Run either `devtools::install_deps()` or `pak::local_install_deps()` to install dependencies
4. This code must be run locally, because any attempt to run it on GitHub Actions failed because Gmail will fail to authenticate. The only way I know of to make the authentication work is to use a service account and delegate domain-wide authority to the Gmail account specifically dedicated to the Museum project that can edit the database sheet and receive RSS from bio/medRxiv. Delegate domain-side authority requires Google Workspace which has to be paid for, which is why I'm not using it.
5. Sign up with you own email address for keyword alerts for the following keywords [on bioRxiv](https://www.biorxiv.org/alerts): `c("spatial transcriptomics", "visium", "merfish", "seqfish", "GeoMX", "ISS", "CosMX", "xenium")`; put the keyword in both the "Name" field and in "Full Text or Abstract or Title", keep everything else default. It's important to put the keyword in the "Name" field because the name is used to get the relevant emails with the Gmail API. Do the same for medRxiv.
6. In order for this code to work, you need to use a desktop email client to automatically convert your bio/medRxiv feed into plain text and forward to museumofst\@gmail.com. This is hacky but necessary because the Gmail API doesn't download the html email bodies, at least when used with the `gmailr` package. At present, I'm already doing this so you don't need to, but in case you want to take over, please let me know and here's how: I use Thunderbird: 
    - I created a filter: click the 3 horizontal line icon at the top right (Mac version) to op  en a menu
    - Go to "Tools", then "Message Filters"
    - Click "New..." to create a filter, select "Match all of the following", put From is cshljnls-mailer\@alerts.highwire.org, to get both bioRxiv and medRxiv emails
    - Uncheck "Manually Run" to make it run automatically
    - Under "Perform these actions:", you can move all matching emails to a separate folder. Mo  st importantly, Forward messages to museumofst\@gmail.com. Automatic forwarding on the Gm  ail website in the browser won't work because it won't convert the html into plain text.
    - You can add other actions such as Mark as Read or delete the messages
7. At present I'm running the code at least once a day to add new entries to the ["to_check" sheet](https://docs.google.com/spreadsheets/d/1sJDb9B7AtYmfKv4-m8XR7uc3XXw_k4kGSout8cqZ8bY/edit?pli=1#gid=305027376). If you wish to take over, please message me and we can discuss. Tokens are required and they are not in this repo.
8. The hard part is over. On the command line, in the directory where you cloned the repo, run `Rscript autolist.R`. Shouldn't take long but it won't finish instantly. That's it.

## About the `to_check` sheet
`date_added`: Not the same as `date_published` in the other sheets in the Museum. This is the date when the feed is received; the papers themselves are usually a bit older.

`pmid`: PubMed ID for the PubMed feed, to make it a bit easier to add the entry to other sheets if relevant.

`URL`: The URL from the feed. This is only used to direct the curator to the paper. The URL in other sheets must be the DOI URL for easier search.

`string_match`: Part of the URL in the `URL` column that can be used to match each entry to existing entries in other sheets, including `irrelevant`, in order to populate the next column. This saves a menial step to manually select the part from the URL.

`existing_sheet`: Indicating the existing sheet(s) where this entry has been found, only applicable to followup versions of preprints. The number in the parenthesis is the first row where this entry is found in the sheet of interest. The date and title are automatically updated. These entries are not removed because occasionally later versions have new datasets and entries previous deemed irrelevant becomes relevant, thus requiring manual intervention. 
