# Turn one of the itb_review file from excel to .tsv format 

library(here)
library(readxl)
library(data.table)
library(dplyr)

itb_excel_path <- here("data", "itb_review", "TWJF-5120-02.revd.Annotated.Neoantigen_Candidates-rev3.xlsx")
sheet_name <- "report-median"

itb_excel_data <- read_excel(itb_excel_path, sheet = sheet_name, col_types = "text")

# remove extra columns
df <- itb_excel_data %>% 
  rename(Evaluation = `Evaluation by lowest`) %>% 
  select(-c(row, review, note)) %>% 
  mutate_all(~ ifelse(. == "--", NA, .)) # change "--" into NA

# there is a row with all NAs, remove this row, but keep other rows that may also have some NA
df_cleaned <- df[rowSums(is.na(df)) < ncol(df), ]

# export to .tsv
outdir <- here("data", "itb_review")
write.table(df_cleaned, file.path(outdir, "TWJF-5120-02.revd.Annotated.Neoantigen_Candidates.tsv"), 
            sep = "\t", row.names = FALSE)