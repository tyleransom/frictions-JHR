# Data Sources to create data set on county and city characteristics

The final output dataset is `AccraData/finaldata/county_city_chars.dta`

Note: ACCRA source data files are removed in order to 

## ACCRA Data:
* 1990q1 - 2008q3
* Generously provided by Chris Timmins 
* For each urban area (of approximatley 290), contains local price indices for Grocery, Housing, Utilities, Transportation, Health Services and Misc
* Contains weights for the basket (i.e. consumption shares) --- these are computed by ACCRA using the CEX

## Housing prices:
* Follow Winters (2008) and use quality-adjusted gross rents from ACS
* Map METAREA variable in ACS to CBSA by case-by-case matching `metarea_CBSA_mapping.do` (see spreadsheet `SA_CBSA_ACS_compare.xlsx` which does a vlookup on a P/MSA to CBSA crosswalk and assigns the MSA in ACS to the prinicpal CBSA)
* Use state average in rural areas for counties not identified in ACS
* Plug in to cobb-douglas formula (see Baum-Snow and Pavan 2012 for details)

## County population:
* 1900-1990 Decennial population by county
* Retrieved from Census Bureau (http://www.census.gov/population/cencounts/1900-90.txt)
* See https://www.census.gov/population/cencounts/00-90doc.txt for technical notes on county transitions, VA Independent Cities, etc.

## County unemployment:
* 1990 - 2011
* Annualized labor force, employed, unemployed for all counties from 1990-2012
* Retrieved from BLS-LAUS (http://www.bls.gov/lau/)

## MSA-level GDP:
* 2001 - 2010
* Annual GDP (same one published in news articles and derived from National Income Product Accounts) measured at MSA level for all Metropolitan and some Micropolitan Statistical Areas
* based on national prices for the goods and services produced within the metropolitan area
* real GDP by metropolitan area does not capture geographic differences in the prices of goods and services that are produced and sold locally
* Listed in Thousands of Constant-2005 dollars (for per-capita levels)
* Listed in Millions of Constant-2005 dollars (for overall levels)
* Unobserved for rural counties
* Retrieved from BEA (http://bea.gov/regional/downloadzip.cfm)

## Land Use:
* Percent land in urbanized area by county
* Also gives percent of population in those urbanized areas
* Retrieved from Census (http://www.census.gov/geo/reference/ua/urban-rural-2010.html), subsection "Percent urban and rural in 2010 by state and county"
* Another file gives Population change in urbanized area from 2000 to 2010 [urban/rural breakdown by county is not available in 2000]
* Retrieve from same URL as above, but the subsection entitled "Changes in urbanized areas from 2000 to 2010"

## Land Use Regulation:
* Wharton Land Use Data
* Gyourko, Joseph, Albert Saiz, and Anita A. Summers, "A New Measure of the Local Regulatory Environment for Housing Markets: The Wharton Residential Land Use Regulatory Index", Urban Studies, 45(3) 693–729, March 2008. 
* Contains a Land Use regulatory index that summarizes how regulated land use is in the metropolitan area

## Latitude/Longitude of MSA (to get distance moved):
* ZIP-code-level coordinates
* Mapped to counties using MABLE/GeoCorr, using the coordinates of the most-populated ZIP code as the coordinates of the county, and the coordinates of the most-populated county as the coordinates of the MSA
* Retrieved from http://federalgovernmentzipcodes.us/

## Crosswalks:
* Source: MABLE/GeoCorr (U of Missouri --- http://mcdc2.missouri.edu/websas/geocorr2k.html)
* CBSAs are as of 2009 delineation
    * CPS city geography variable converted to FIPS county (`county_cps_msa_xwalk.csv`)
    * CPS city geography variable converted to CBSA code (`cbsa_cps_msa_xwalk.csv`)
    * County to CBSA code, **showing land area** (`countyxwalkarea.csv`)
    * County to CBSA code, **showing 2000 population** (`countyxwalkpop00.csv`)
* also available: # of housing units in 2000 Census [could be useful for housing supply analysis]
