capture log close
set more off
log using "create_CPS_CBSA_xwalk.log", replace

tempfile holder1
insheet using ipums_codes.csv, comma case clear
save `holder1'

insheet using cbsa_codes.csv, comma case clear
merge 1:1 name using `holder1'

log close
exit
