# NOTE: For execution, type "nohup bash data_create.bash &"
# sleep 10h

cd [REDACTED]
gunzip [REDACTED]county_city_chars.csv.gz
nohup sas                   citysize.sas
gzip   [REDACTED]county_city_chars.csv

nohup sas                   create2004jmp.sas
nohup sas                   create2008jmp.sas
nohup sas                   export_stata.sas

test -f [REDACTED]longcombinedpuf2004.sas7bdat && gzip -f [REDACTED]longcombinedpuf2004.sas7bdat || echo "longcombinedpuf2004.sas7bdat Not Found"
# sleep 30
test -f [REDACTED]sipp2004NHWmale.csv && gzip -f [REDACTED]sipp2004NHWmale.csv || echo "sipp2004NHWmale.csv Not Found"
# sleep 30
test -f [REDACTED]longcombinedpuf2008.sas7bdat && gzip -f [REDACTED]longcombinedpuf2008.sas7bdat || echo "longcombinedpuf2008.sas7bdat Not Found"
# sleep 30
test -f [REDACTED]sipp2008NHWmale.csv && gzip -f [REDACTED]sipp2008NHWmale.csv || echo "sipp2008NHWmale.csv Not Found"
# sleep 30
statamphup cr_clean_SIPP.do
nohup stata-mp -b do cr_clean_SIPP_2008.do
nohup stata-mp -b do panel_stacker.do
mathup [REDACTED]dataImportCombinedAnnual.m