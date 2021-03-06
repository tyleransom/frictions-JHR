version 14.1
clear all
capture log close
set more off
log using "ACS_data_create.log", replace

! gunzip ACS2005raw.dat.gz

clear
quietly infix                   ///
  int     year         1-4      ///
  byte    datanum      5-6      ///
  double  serial       7-14     ///
  float   hhwt         15-24    ///
  byte    region       25-26    ///
  byte    stateicp     27-28    ///
  byte    statefip     29-30    ///
  int     county       31-34    ///
  byte    metro        35-35    ///
  int     metarea      36-38    ///
  int     metaread     39-42    ///
  int     city         43-46    ///
  long    citypop      47-51    ///
  long    puma         52-56    ///
  long    pumares2mig  57-61    ///
  long    pumasupr     62-66    ///
  int     conspuma     67-69    ///
  byte    appal        70-70    ///
  byte    appald       71-72    ///
  byte    homeland     73-73    ///
  int     cntry        74-76    ///
  byte    gq           77-77    ///
  byte    acrehous     78-78    ///
  int     rentgrs      79-82    ///
  byte    kitchen      83-83    ///
  byte    rooms        84-85    ///
  byte    plumbing     86-87    ///
  byte    builtyr2     88-89    ///
  byte    unitsstr     90-91    ///
  byte    bedrooms     92-93    ///
  byte    qacrehou     94-94    ///
  byte    qbedroom     95-95    ///
  byte    qbuilty2     96-96    ///
  byte    qkitchen     97-97    ///
  byte    qrooms       98-98    ///
  byte    qunitsst     99-99    ///
  int     pernum       100-103  ///
  float   perwt        104-113  ///
  using `"ACS2005raw.dat"'

! gzip -f ACS2005raw.dat.gz

replace hhwt        = hhwt        / 100
replace perwt       = perwt       / 100

format serial      %8.0f
format hhwt        %10.2f
format perwt       %10.2f

label var year        `"Census year"'
label var datanum     `"Data set number"'
label var serial      `"Household serial number"'
label var hhwt        `"Household weight"'
label var region      `"Census region and division"'
label var stateicp    `"State (ICPSR code)"'
label var statefip    `"State (FIPS code)"'
label var county      `"County"'
label var metro       `"Metropolitan status"'
label var metarea     `"Metropolitan area [general version]"'
label var metaread    `"Metropolitan area [detailed version]"'
label var city        `"City"'
label var citypop     `"City population"'
label var puma        `"Public Use Microdata Area"'
label var pumares2mig `"Public Use Microdata Area matching MIGPUMA"'
label var pumasupr    `"Super Public Use Microdata Area"'
label var conspuma    `"Consistent Public Use Microdata Area"'
label var appal       `"Appalachian region [general version]"'
label var appald      `"Appalachian region [detailed version]"'
label var homeland    `"American Indian, Alaska Native, or Native Hawaiian homeland area"'
label var cntry       `"Country"'
label var gq          `"Group quarters status"'
label var acrehous    `"House acreage"'
label var rentgrs     `"Monthly gross rent"'
label var kitchen     `"Kitchen or cooking facilities"'
label var rooms       `"Number of rooms"'
label var plumbing    `"Plumbing facilities"'
label var builtyr2    `"Age of structure, decade"'
label var unitsstr    `"Units in structure"'
label var bedrooms    `"Number of bedrooms"'
label var qacrehou    `"Flag for Acrehous"'
label var qbedroom    `"Flag for Bedrooms"'
label var qbuilty2    `"Flag for Builtyr2"'
label var qkitchen    `"Flag for Kitchen"'
label var qrooms      `"Flag for Rooms"'
label var qunitsst    `"Flag for Unitsstr"'
label var pernum      `"Person number in sample unit"'
label var perwt       `"Person weight"'

label define year_lbl 1850 `"1850"'
label define year_lbl 1860 `"1860"', add
label define year_lbl 1870 `"1870"', add
label define year_lbl 1880 `"1880"', add
label define year_lbl 1900 `"1900"', add
label define year_lbl 1910 `"1910"', add
label define year_lbl 1920 `"1920"', add
label define year_lbl 1930 `"1930"', add
label define year_lbl 1940 `"1940"', add
label define year_lbl 1950 `"1950"', add
label define year_lbl 1960 `"1960"', add
label define year_lbl 1970 `"1970"', add
label define year_lbl 1980 `"1980"', add
label define year_lbl 1990 `"1990"', add
label define year_lbl 2000 `"2000"', add
label define year_lbl 2001 `"2001"', add
label define year_lbl 2002 `"2002"', add
label define year_lbl 2003 `"2003"', add
label define year_lbl 2004 `"2004"', add
label define year_lbl 2005 `"2005"', add
label define year_lbl 2006 `"2006"', add
label define year_lbl 2007 `"2007"', add
label define year_lbl 2008 `"2008"', add
label define year_lbl 2009 `"2009"', add
label define year_lbl 2010 `"2010"', add
label define year_lbl 2011 `"2011"', add
label values year year_lbl

label define region_lbl 11 `"New England Division"'
label define region_lbl 12 `"Middle Atlantic Division"', add
label define region_lbl 13 `"Mixed Northeast Divisions (1970 Metro)"', add
label define region_lbl 21 `"East North Central Div."', add
label define region_lbl 22 `"West North Central Div."', add
label define region_lbl 23 `"Mixed Midwest Divisions (1970 Metro)"', add
label define region_lbl 31 `"South Atlantic Division"', add
label define region_lbl 32 `"East South Central Div."', add
label define region_lbl 33 `"West South Central Div."', add
label define region_lbl 34 `"Mixed Southern Divisions (1970 Metro)"', add
label define region_lbl 41 `"Mountain Division"', add
label define region_lbl 42 `"Pacific Division"', add
label define region_lbl 43 `"Mixed Western Divisions (1970 Metro)"', add
label define region_lbl 91 `"Military/Military reservations"', add
label define region_lbl 92 `"PUMA boundaries cross state lines-1% sample"', add
label define region_lbl 97 `"State not identified"', add
label define region_lbl 99 `"Not identified"', add
label values region region_lbl

label define stateicp_lbl 01 `"Connecticut"'
label define stateicp_lbl 02 `"Maine"', add
label define stateicp_lbl 03 `"Massachusetts"', add
label define stateicp_lbl 04 `"New Hampshire"', add
label define stateicp_lbl 05 `"Rhode Island"', add
label define stateicp_lbl 06 `"Vermont"', add
label define stateicp_lbl 11 `"Delaware"', add
label define stateicp_lbl 12 `"New Jersey"', add
label define stateicp_lbl 13 `"New York"', add
label define stateicp_lbl 14 `"Pennsylvania"', add
label define stateicp_lbl 21 `"Illinois"', add
label define stateicp_lbl 22 `"Indiana"', add
label define stateicp_lbl 23 `"Michigan"', add
label define stateicp_lbl 24 `"Ohio"', add
label define stateicp_lbl 25 `"Wisconsin"', add
label define stateicp_lbl 31 `"Iowa"', add
label define stateicp_lbl 32 `"Kansas"', add
label define stateicp_lbl 33 `"Minnesota"', add
label define stateicp_lbl 34 `"Missouri"', add
label define stateicp_lbl 35 `"Nebraska"', add
label define stateicp_lbl 36 `"North Dakota"', add
label define stateicp_lbl 37 `"South Dakota"', add
label define stateicp_lbl 40 `"Virginia"', add
label define stateicp_lbl 41 `"Alabama"', add
label define stateicp_lbl 42 `"Arkansas"', add
label define stateicp_lbl 43 `"Florida"', add
label define stateicp_lbl 44 `"Georgia"', add
label define stateicp_lbl 45 `"Louisiana"', add
label define stateicp_lbl 46 `"Mississippi"', add
label define stateicp_lbl 47 `"North Carolina"', add
label define stateicp_lbl 48 `"South Carolina"', add
label define stateicp_lbl 49 `"Texas"', add
label define stateicp_lbl 51 `"Kentucky"', add
label define stateicp_lbl 52 `"Maryland"', add
label define stateicp_lbl 53 `"Oklahoma"', add
label define stateicp_lbl 54 `"Tennessee"', add
label define stateicp_lbl 56 `"West Virginia"', add
label define stateicp_lbl 61 `"Arizona"', add
label define stateicp_lbl 62 `"Colorado"', add
label define stateicp_lbl 63 `"Idaho"', add
label define stateicp_lbl 64 `"Montana"', add
label define stateicp_lbl 65 `"Nevada"', add
label define stateicp_lbl 66 `"New Mexico"', add
label define stateicp_lbl 67 `"Utah"', add
label define stateicp_lbl 68 `"Wyoming"', add
label define stateicp_lbl 71 `"California"', add
label define stateicp_lbl 72 `"Oregon"', add
label define stateicp_lbl 73 `"Washington"', add
label define stateicp_lbl 81 `"Alaska"', add
label define stateicp_lbl 82 `"Hawaii"', add
label define stateicp_lbl 83 `"Puerto Rico"', add
label define stateicp_lbl 96 `"State groupings (1980 Urban/rural sample)"', add
label define stateicp_lbl 97 `"Military/Mil. Reservations"', add
label define stateicp_lbl 98 `"District of Columbia"', add
label define stateicp_lbl 99 `"State not identified"', add
label values stateicp stateicp_lbl

label define statefip_lbl 01 `"Alabama"'
label define statefip_lbl 02 `"Alaska"', add
label define statefip_lbl 04 `"Arizona"', add
label define statefip_lbl 05 `"Arkansas"', add
label define statefip_lbl 06 `"California"', add
label define statefip_lbl 08 `"Colorado"', add
label define statefip_lbl 09 `"Connecticut"', add
label define statefip_lbl 10 `"Delaware"', add
label define statefip_lbl 11 `"District of Columbia"', add
label define statefip_lbl 12 `"Florida"', add
label define statefip_lbl 13 `"Georgia"', add
label define statefip_lbl 15 `"Hawaii"', add
label define statefip_lbl 16 `"Idaho"', add
label define statefip_lbl 17 `"Illinois"', add
label define statefip_lbl 18 `"Indiana"', add
label define statefip_lbl 19 `"Iowa"', add
label define statefip_lbl 20 `"Kansas"', add
label define statefip_lbl 21 `"Kentucky"', add
label define statefip_lbl 22 `"Louisiana"', add
label define statefip_lbl 23 `"Maine"', add
label define statefip_lbl 24 `"Maryland"', add
label define statefip_lbl 25 `"Massachusetts"', add
label define statefip_lbl 26 `"Michigan"', add
label define statefip_lbl 27 `"Minnesota"', add
label define statefip_lbl 28 `"Mississippi"', add
label define statefip_lbl 29 `"Missouri"', add
label define statefip_lbl 30 `"Montana"', add
label define statefip_lbl 31 `"Nebraska"', add
label define statefip_lbl 32 `"Nevada"', add
label define statefip_lbl 33 `"New Hampshire"', add
label define statefip_lbl 34 `"New Jersey"', add
label define statefip_lbl 35 `"New Mexico"', add
label define statefip_lbl 36 `"New York"', add
label define statefip_lbl 37 `"North Carolina"', add
label define statefip_lbl 38 `"North Dakota"', add
label define statefip_lbl 39 `"Ohio"', add
label define statefip_lbl 40 `"Oklahoma"', add
label define statefip_lbl 41 `"Oregon"', add
label define statefip_lbl 42 `"Pennsylvania"', add
label define statefip_lbl 44 `"Rhode Island"', add
label define statefip_lbl 45 `"South Carolina"', add
label define statefip_lbl 46 `"South Dakota"', add
label define statefip_lbl 47 `"Tennessee"', add
label define statefip_lbl 48 `"Texas"', add
label define statefip_lbl 49 `"Utah"', add
label define statefip_lbl 50 `"Vermont"', add
label define statefip_lbl 51 `"Virginia"', add
label define statefip_lbl 53 `"Washington"', add
label define statefip_lbl 54 `"West Virginia"', add
label define statefip_lbl 55 `"Wisconsin"', add
label define statefip_lbl 56 `"Wyoming"', add
label define statefip_lbl 61 `"Maine-New Hampshire-Vermont"', add
label define statefip_lbl 62 `"Massachusetts-Rhode Island"', add
label define statefip_lbl 63 `"Minnesota-Iowa-Missouri-Kansas-Nebraska-S.Dakota-N.Dakota"', add
label define statefip_lbl 64 `"Maryland-Delaware"', add
label define statefip_lbl 65 `"Montana-Idaho-Wyoming"', add
label define statefip_lbl 66 `"Utah-Nevada"', add
label define statefip_lbl 67 `"Arizona-New Mexico"', add
label define statefip_lbl 68 `"Alaska-Hawaii"', add
label define statefip_lbl 72 `"Puerto Rico"', add
label define statefip_lbl 97 `"Military/Mil. Reservation"', add
label define statefip_lbl 99 `"State not identified"', add
label values statefip statefip_lbl

label define county_lbl 0010 `"0010"'
label define county_lbl 0030 `"0030"', add
label define county_lbl 0050 `"0050"', add
label define county_lbl 0070 `"0070"', add
label define county_lbl 0090 `"0090"', add
label define county_lbl 0110 `"0110"', add
label define county_lbl 0130 `"0130"', add
label define county_lbl 0150 `"0150"', add
label define county_lbl 0170 `"0170"', add
label define county_lbl 0190 `"0190"', add
label define county_lbl 0210 `"0210"', add
label define county_lbl 0230 `"0230"', add
label define county_lbl 0250 `"0250"', add
label define county_lbl 0270 `"0270"', add
label define county_lbl 0290 `"0290"', add
label define county_lbl 0310 `"0310"', add
label define county_lbl 0330 `"0330"', add
label define county_lbl 0350 `"0350"', add
label define county_lbl 0360 `"0360"', add
label define county_lbl 0370 `"0370"', add
label define county_lbl 0390 `"0390"', add
label define county_lbl 0410 `"0410"', add
label define county_lbl 0430 `"0430"', add
label define county_lbl 0450 `"0450"', add
label define county_lbl 0455 `"0455"', add
label define county_lbl 0470 `"0470"', add
label define county_lbl 0490 `"0490"', add
label define county_lbl 0510 `"0510"', add
label define county_lbl 0530 `"0530"', add
label define county_lbl 0550 `"0550"', add
label define county_lbl 0570 `"0570"', add
label define county_lbl 0590 `"0590"', add
label define county_lbl 0605 `"0605"', add
label define county_lbl 0610 `"0610"', add
label define county_lbl 0630 `"0630"', add
label define county_lbl 0650 `"0650"', add
label define county_lbl 0670 `"0670"', add
label define county_lbl 0690 `"0690"', add
label define county_lbl 0710 `"0710"', add
label define county_lbl 0730 `"0730"', add
label define county_lbl 0750 `"0750"', add
label define county_lbl 0770 `"0770"', add
label define county_lbl 0790 `"0790"', add
label define county_lbl 0810 `"0810"', add
label define county_lbl 0830 `"0830"', add
label define county_lbl 0850 `"0850"', add
label define county_lbl 0870 `"0870"', add
label define county_lbl 0890 `"0890"', add
label define county_lbl 0910 `"0910"', add
label define county_lbl 0930 `"0930"', add
label define county_lbl 0950 `"0950"', add
label define county_lbl 0970 `"0970"', add
label define county_lbl 0990 `"0990"', add
label define county_lbl 1010 `"1010"', add
label define county_lbl 1030 `"1030"', add
label define county_lbl 1050 `"1050"', add
label define county_lbl 1070 `"1070"', add
label define county_lbl 1090 `"1090"', add
label define county_lbl 1110 `"1110"', add
label define county_lbl 1130 `"1130"', add
label define county_lbl 1150 `"1150"', add
label define county_lbl 1170 `"1170"', add
label define county_lbl 1190 `"1190"', add
label define county_lbl 1210 `"1210"', add
label define county_lbl 1230 `"1230"', add
label define county_lbl 1250 `"1250"', add
label define county_lbl 1270 `"1270"', add
label define county_lbl 1290 `"1290"', add
label define county_lbl 1310 `"1310"', add
label define county_lbl 1330 `"1330"', add
label define county_lbl 1350 `"1350"', add
label define county_lbl 1370 `"1370"', add
label define county_lbl 1390 `"1390"', add
label define county_lbl 1410 `"1410"', add
label define county_lbl 1430 `"1430"', add
label define county_lbl 1450 `"1450"', add
label define county_lbl 1470 `"1470"', add
label define county_lbl 1490 `"1490"', add
label define county_lbl 1510 `"1510"', add
label define county_lbl 1530 `"1530"', add
label define county_lbl 1550 `"1550"', add
label define county_lbl 1570 `"1570"', add
label define county_lbl 1590 `"1590"', add
label define county_lbl 1610 `"1610"', add
label define county_lbl 1630 `"1630"', add
label define county_lbl 1650 `"1650"', add
label define county_lbl 1670 `"1670"', add
label define county_lbl 1690 `"1690"', add
label define county_lbl 1710 `"1710"', add
label define county_lbl 1730 `"1730"', add
label define county_lbl 1750 `"1750"', add
label define county_lbl 1770 `"1770"', add
label define county_lbl 1790 `"1790"', add
label define county_lbl 1810 `"1810"', add
label define county_lbl 1830 `"1830"', add
label define county_lbl 1850 `"1850"', add
label define county_lbl 1870 `"1870"', add
label define county_lbl 1875 `"1875"', add
label define county_lbl 1890 `"1890"', add
label define county_lbl 1910 `"1910"', add
label define county_lbl 1930 `"1930"', add
label define county_lbl 1950 `"1950"', add
label define county_lbl 1970 `"1970"', add
label define county_lbl 1990 `"1990"', add
label define county_lbl 2010 `"2010"', add
label define county_lbl 2030 `"2030"', add
label define county_lbl 2050 `"2050"', add
label define county_lbl 2070 `"2070"', add
label define county_lbl 2090 `"2090"', add
label define county_lbl 2110 `"2110"', add
label define county_lbl 2130 `"2130"', add
label define county_lbl 2150 `"2150"', add
label define county_lbl 2170 `"2170"', add
label define county_lbl 2190 `"2190"', add
label define county_lbl 2210 `"2210"', add
label define county_lbl 2230 `"2230"', add
label define county_lbl 2250 `"2250"', add
label define county_lbl 2270 `"2270"', add
label define county_lbl 2290 `"2290"', add
label define county_lbl 2310 `"2310"', add
label define county_lbl 2330 `"2330"', add
label define county_lbl 2350 `"2350"', add
label define county_lbl 2370 `"2370"', add
label define county_lbl 2390 `"2390"', add
label define county_lbl 2410 `"2410"', add
label define county_lbl 2430 `"2430"', add
label define county_lbl 2450 `"2450"', add
label define county_lbl 2470 `"2470"', add
label define county_lbl 2490 `"2490"', add
label define county_lbl 2510 `"2510"', add
label define county_lbl 2530 `"2530"', add
label define county_lbl 2550 `"2550"', add
label define county_lbl 2570 `"2570"', add
label define county_lbl 2590 `"2590"', add
label define county_lbl 2610 `"2610"', add
label define county_lbl 2630 `"2630"', add
label define county_lbl 2650 `"2650"', add
label define county_lbl 2670 `"2670"', add
label define county_lbl 2690 `"2690"', add
label define county_lbl 2710 `"2710"', add
label define county_lbl 2730 `"2730"', add
label define county_lbl 2750 `"2750"', add
label define county_lbl 2770 `"2770"', add
label define county_lbl 2790 `"2790"', add
label define county_lbl 2810 `"2810"', add
label define county_lbl 2830 `"2830"', add
label define county_lbl 2850 `"2850"', add
label define county_lbl 2870 `"2870"', add
label define county_lbl 2890 `"2890"', add
label define county_lbl 2910 `"2910"', add
label define county_lbl 2930 `"2930"', add
label define county_lbl 2950 `"2950"', add
label define county_lbl 2970 `"2970"', add
label define county_lbl 2990 `"2990"', add
label define county_lbl 3010 `"3010"', add
label define county_lbl 3030 `"3030"', add
label define county_lbl 3050 `"3050"', add
label define county_lbl 3070 `"3070"', add
label define county_lbl 3090 `"3090"', add
label define county_lbl 3110 `"3110"', add
label define county_lbl 3130 `"3130"', add
label define county_lbl 3150 `"3150"', add
label define county_lbl 3170 `"3170"', add
label define county_lbl 3190 `"3190"', add
label define county_lbl 3210 `"3210"', add
label define county_lbl 3230 `"3230"', add
label define county_lbl 3250 `"3250"', add
label define county_lbl 3270 `"3270"', add
label define county_lbl 3290 `"3290"', add
label define county_lbl 3310 `"3310"', add
label define county_lbl 3330 `"3330"', add
label define county_lbl 3350 `"3350"', add
label define county_lbl 3370 `"3370"', add
label define county_lbl 3390 `"3390"', add
label define county_lbl 3410 `"3410"', add
label define county_lbl 3430 `"3430"', add
label define county_lbl 3450 `"3450"', add
label define county_lbl 3470 `"3470"', add
label define county_lbl 3490 `"3490"', add
label define county_lbl 3510 `"3510"', add
label define county_lbl 3530 `"3530"', add
label define county_lbl 3550 `"3550"', add
label define county_lbl 3570 `"3570"', add
label define county_lbl 3590 `"3590"', add
label define county_lbl 3610 `"3610"', add
label define county_lbl 3630 `"3630"', add
label define county_lbl 3650 `"3650"', add
label define county_lbl 3670 `"3670"', add
label define county_lbl 3690 `"3690"', add
label define county_lbl 3710 `"3710"', add
label define county_lbl 3730 `"3730"', add
label define county_lbl 3750 `"3750"', add
label define county_lbl 3770 `"3770"', add
label define county_lbl 3790 `"3790"', add
label define county_lbl 3810 `"3810"', add
label define county_lbl 3830 `"3830"', add
label define county_lbl 3850 `"3850"', add
label define county_lbl 3870 `"3870"', add
label define county_lbl 3890 `"3890"', add
label define county_lbl 3910 `"3910"', add
label define county_lbl 3930 `"3930"', add
label define county_lbl 3950 `"3950"', add
label define county_lbl 3970 `"3970"', add
label define county_lbl 3990 `"3990"', add
label define county_lbl 4010 `"4010"', add
label define county_lbl 4030 `"4030"', add
label define county_lbl 4050 `"4050"', add
label define county_lbl 4070 `"4070"', add
label define county_lbl 4090 `"4090"', add
label define county_lbl 4110 `"4110"', add
label define county_lbl 4130 `"4130"', add
label define county_lbl 4150 `"4150"', add
label define county_lbl 4170 `"4170"', add
label define county_lbl 4190 `"4190"', add
label define county_lbl 4210 `"4210"', add
label define county_lbl 4230 `"4230"', add
label define county_lbl 4250 `"4250"', add
label define county_lbl 4270 `"4270"', add
label define county_lbl 4290 `"4290"', add
label define county_lbl 4310 `"4310"', add
label define county_lbl 4330 `"4330"', add
label define county_lbl 4350 `"4350"', add
label define county_lbl 4370 `"4370"', add
label define county_lbl 4390 `"4390"', add
label define county_lbl 4410 `"4410"', add
label define county_lbl 4430 `"4430"', add
label define county_lbl 4450 `"4450"', add
label define county_lbl 4470 `"4470"', add
label define county_lbl 4490 `"4490"', add
label define county_lbl 4510 `"4510"', add
label define county_lbl 4530 `"4530"', add
label define county_lbl 4550 `"4550"', add
label define county_lbl 4570 `"4570"', add
label define county_lbl 4590 `"4590"', add
label define county_lbl 4610 `"4610"', add
label define county_lbl 4630 `"4630"', add
label define county_lbl 4650 `"4650"', add
label define county_lbl 4670 `"4670"', add
label define county_lbl 4690 `"4690"', add
label define county_lbl 4710 `"4710"', add
label define county_lbl 4730 `"4730"', add
label define county_lbl 4750 `"4750"', add
label define county_lbl 4770 `"4770"', add
label define county_lbl 4790 `"4790"', add
label define county_lbl 4810 `"4810"', add
label define county_lbl 4830 `"4830"', add
label define county_lbl 4850 `"4850"', add
label define county_lbl 4870 `"4870"', add
label define county_lbl 4890 `"4890"', add
label define county_lbl 4910 `"4910"', add
label define county_lbl 4930 `"4930"', add
label define county_lbl 4950 `"4950"', add
label define county_lbl 4970 `"4970"', add
label define county_lbl 4990 `"4990"', add
label define county_lbl 5010 `"5010"', add
label define county_lbl 5030 `"5030"', add
label define county_lbl 5050 `"5050"', add
label define county_lbl 5070 `"5070"', add
label define county_lbl 5100 `"5100"', add
label define county_lbl 5200 `"5200"', add
label define county_lbl 5300 `"5300"', add
label define county_lbl 5400 `"5400"', add
label define county_lbl 5500 `"5500"', add
label define county_lbl 5600 `"5600"', add
label define county_lbl 5700 `"5700"', add
label define county_lbl 5800 `"5800"', add
label define county_lbl 5900 `"5900"', add
label define county_lbl 6100 `"6100"', add
label define county_lbl 6300 `"6300"', add
label define county_lbl 6400 `"6400"', add
label define county_lbl 6500 `"6500"', add
label define county_lbl 6600 `"6600"', add
label define county_lbl 6700 `"6700"', add
label define county_lbl 6800 `"6800"', add
label define county_lbl 6900 `"6900"', add
label define county_lbl 7000 `"7000"', add
label define county_lbl 7100 `"7100"', add
label define county_lbl 7200 `"7200"', add
label define county_lbl 7300 `"7300"', add
label define county_lbl 7400 `"7400"', add
label define county_lbl 7500 `"7500"', add
label define county_lbl 7600 `"7600"', add
label define county_lbl 7700 `"7700"', add
label define county_lbl 7800 `"7800"', add
label define county_lbl 7850 `"7850"', add
label define county_lbl 7900 `"7900"', add
label define county_lbl 8000 `"8000"', add
label define county_lbl 8100 `"8100"', add
label define county_lbl 8200 `"8200"', add
label define county_lbl 8300 `"8300"', add
label define county_lbl 8400 `"8400"', add
label values county county_lbl

label define metro_lbl 0 `"Not identifiable"'
label define metro_lbl 1 `"Not in metro area"', add
label define metro_lbl 2 `"In metro area, central city"', add
label define metro_lbl 3 `"In metro, area, outside central city"', add
label define metro_lbl 4 `"Central city status unknown"', add
label values metro metro_lbl

label define metarea_lbl 000 `"Not identifiable or not in an MSA"'
label define metarea_lbl 004 `"Abilene, TX"', add
label define metarea_lbl 006 `"Aguadilla, PR"', add
label define metarea_lbl 008 `"Akron, OH"', add
label define metarea_lbl 012 `"Albany, GA"', add
label define metarea_lbl 016 `"Albany-Schenectady-Troy, NY"', add
label define metarea_lbl 020 `"Albuquerque, NM"', add
label define metarea_lbl 022 `"Alexandria, LA"', add
label define metarea_lbl 024 `"Allentown-Bethlehem-Easton, PA/NJ"', add
label define metarea_lbl 028 `"Altoona, PA"', add
label define metarea_lbl 032 `"Amarillo, TX"', add
label define metarea_lbl 038 `"Anchorage, AK"', add
label define metarea_lbl 040 `"Anderson, IN"', add
label define metarea_lbl 044 `"Ann Arbor, MI"', add
label define metarea_lbl 045 `"Anniston, AL"', add
label define metarea_lbl 046 `"Appleton-Oshkosh-Neenah, WI"', add
label define metarea_lbl 047 `"Arecibo, PR"', add
label define metarea_lbl 048 `"Asheville, NC"', add
label define metarea_lbl 050 `"Athens, GA"', add
label define metarea_lbl 052 `"Atlanta, GA"', add
label define metarea_lbl 056 `"Atlantic City, NJ"', add
label define metarea_lbl 058 `"Auburn-Opekika, AL"', add
label define metarea_lbl 060 `"Augusta-Aiken, GA-SC"', add
label define metarea_lbl 064 `"Austin, TX"', add
label define metarea_lbl 068 `"Bakersfield, CA"', add
label define metarea_lbl 072 `"Baltimore, MD"', add
label define metarea_lbl 073 `"Bangor, ME"', add
label define metarea_lbl 074 `"Barnstable-Yarmouth, MA"', add
label define metarea_lbl 076 `"Baton Rouge, LA"', add
label define metarea_lbl 078 `"Battle Creek, MI"', add
label define metarea_lbl 084 `"Beaumont-Port Arthur-Orange,TX"', add
label define metarea_lbl 086 `"Bellingham, WA"', add
label define metarea_lbl 087 `"Benton Harbor, MI"', add
label define metarea_lbl 088 `"Billings, MT"', add
label define metarea_lbl 092 `"Biloxi-Gulfport, MS"', add
label define metarea_lbl 096 `"Binghamton, NY"', add
label define metarea_lbl 100 `"Birmingham, AL"', add
label define metarea_lbl 102 `"Bloomington, IN"', add
label define metarea_lbl 104 `"Bloomington-Normal, IL"', add
label define metarea_lbl 108 `"Boise City, ID"', add
label define metarea_lbl 112 `"Boston, MA-NH"', add
label define metarea_lbl 114 `"Bradenton, FL"', add
label define metarea_lbl 115 `"Bremerton, WA"', add
label define metarea_lbl 116 `"Bridgeport, CT"', add
label define metarea_lbl 120 `"Brockton, MA"', add
label define metarea_lbl 124 `"Brownsville-Harlingen-San Benito, TX"', add
label define metarea_lbl 126 `"Bryan-College Station, TX"', add
label define metarea_lbl 128 `"Buffalo-Niagara Falls, NY"', add
label define metarea_lbl 130 `"Burlington, NC"', add
label define metarea_lbl 131 `"Burlington, VT"', add
label define metarea_lbl 132 `"Canton, OH"', add
label define metarea_lbl 133 `"Caguas, PR"', add
label define metarea_lbl 135 `"Casper, WY"', add
label define metarea_lbl 136 `"Cedar Rapids, IA"', add
label define metarea_lbl 140 `"Champaign-Urbana-Rantoul, IL"', add
label define metarea_lbl 144 `"Charleston-N.Charleston,SC"', add
label define metarea_lbl 148 `"Charleston, WV"', add
label define metarea_lbl 152 `"Charlotte-Gastonia-Rock Hill, NC-SC"', add
label define metarea_lbl 154 `"Charlottesville, VA"', add
label define metarea_lbl 156 `"Chattanooga, TN/GA"', add
label define metarea_lbl 158 `"Cheyenne, WY"', add
label define metarea_lbl 160 `"Chicago, IL"', add
label define metarea_lbl 162 `"Chico, CA"', add
label define metarea_lbl 164 `"Cincinnati-Hamilton, OH/KY/IN"', add
label define metarea_lbl 166 `"Clarksville- Hopkinsville, TN/KY"', add
label define metarea_lbl 168 `"Cleveland, OH"', add
label define metarea_lbl 172 `"Colorado Springs, CO"', add
label define metarea_lbl 174 `"Columbia, MO"', add
label define metarea_lbl 176 `"Columbia, SC"', add
label define metarea_lbl 180 `"Columbus, GA/AL"', add
label define metarea_lbl 184 `"Columbus, OH"', add
label define metarea_lbl 188 `"Corpus Christi, TX"', add
label define metarea_lbl 190 `"Cumberland, MD/WV"', add
label define metarea_lbl 192 `"Dallas-Fort Worth, TX"', add
label define metarea_lbl 193 `"Danbury, CT"', add
label define metarea_lbl 195 `"Danville, VA"', add
label define metarea_lbl 196 `"Davenport, IA-Rock Island -Moline, IL"', add
label define metarea_lbl 200 `"Dayton-Springfield, OH"', add
label define metarea_lbl 202 `"Daytona Beach, FL"', add
label define metarea_lbl 203 `"Decatur, AL"', add
label define metarea_lbl 204 `"Decatur, IL"', add
label define metarea_lbl 208 `"Denver-Boulder, CO"', add
label define metarea_lbl 212 `"Des Moines, IA"', add
label define metarea_lbl 216 `"Detroit, MI"', add
label define metarea_lbl 218 `"Dothan, AL"', add
label define metarea_lbl 219 `"Dover, DE"', add
label define metarea_lbl 220 `"Dubuque, IA"', add
label define metarea_lbl 224 `"Duluth-Superior, MN/WI"', add
label define metarea_lbl 228 `"Dutchess Co., NY"', add
label define metarea_lbl 229 `"Eau Claire, WI"', add
label define metarea_lbl 231 `"El Paso, TX"', add
label define metarea_lbl 232 `"Elkhart-Goshen, IN"', add
label define metarea_lbl 233 `"Elmira, NY"', add
label define metarea_lbl 234 `"Enid, OK"', add
label define metarea_lbl 236 `"Erie, PA"', add
label define metarea_lbl 240 `"Eugene-Springfield, OR"', add
label define metarea_lbl 244 `"Evansville, IN/KY"', add
label define metarea_lbl 252 `"Fargo-Morehead, ND/MN"', add
label define metarea_lbl 256 `"Fayetteville, NC"', add
label define metarea_lbl 258 `"Fayetteville-Springdale, AR"', add
label define metarea_lbl 260 `"Fitchburg-Leominster, MA"', add
label define metarea_lbl 262 `"Flagstaff, AZ-UT"', add
label define metarea_lbl 264 `"Flint, MI"', add
label define metarea_lbl 265 `"Florence, AL"', add
label define metarea_lbl 266 `"Florence, SC"', add
label define metarea_lbl 267 `"Fort Collins-Loveland, CO"', add
label define metarea_lbl 268 `"Fort Lauderdale-Hollywood-Pompano Beach, FL"', add
label define metarea_lbl 270 `"Fort Myers-Cape Coral, FL"', add
label define metarea_lbl 271 `"Fort Pierce, FL"', add
label define metarea_lbl 272 `"Fort Smith, AR/OK"', add
label define metarea_lbl 275 `"Fort Walton Beach, FL"', add
label define metarea_lbl 276 `"Fort Wayne, IN"', add
label define metarea_lbl 284 `"Fresno, CA"', add
label define metarea_lbl 288 `"Gadsden, AL"', add
label define metarea_lbl 290 `"Gainesville, FL"', add
label define metarea_lbl 292 `"Galveston-Texas City, TX"', add
label define metarea_lbl 297 `"Glens Falls, NY"', add
label define metarea_lbl 298 `"Goldsboro, NC"', add
label define metarea_lbl 299 `"Grand Forks, ND"', add
label define metarea_lbl 300 `"Grand Rapids, MI"', add
label define metarea_lbl 301 `"Grand Junction, CO"', add
label define metarea_lbl 304 `"Great Falls, MT"', add
label define metarea_lbl 306 `"Greeley, CO"', add
label define metarea_lbl 308 `"Green Bay, WI"', add
label define metarea_lbl 312 `"Greensboro-Winston Salem-High Point, NC"', add
label define metarea_lbl 315 `"Greenville, NC"', add
label define metarea_lbl 316 `"Greenville-Spartanburg-Anderson SC"', add
label define metarea_lbl 318 `"Hagerstown, MD"', add
label define metarea_lbl 320 `"Hamilton-Middleton, OH"', add
label define metarea_lbl 324 `"Harrisburg-Lebanon--Carlisle, PA"', add
label define metarea_lbl 328 `"Hartford-Bristol-Middleton- New Britain, CT"', add
label define metarea_lbl 329 `"Hickory-Morgantown, NC"', add
label define metarea_lbl 330 `"Hattiesburg, MS"', add
label define metarea_lbl 332 `"Honolulu, HI"', add
label define metarea_lbl 335 `"Houma-Thibodoux, LA"', add
label define metarea_lbl 336 `"Houston-Brazoria, TX"', add
label define metarea_lbl 340 `"Huntington-Ashland, WV/KY/OH"', add
label define metarea_lbl 344 `"Huntsville, AL"', add
label define metarea_lbl 348 `"Indianapolis, IN"', add
label define metarea_lbl 350 `"Iowa City, IA"', add
label define metarea_lbl 352 `"Jackson, MI"', add
label define metarea_lbl 356 `"Jackson, MS"', add
label define metarea_lbl 358 `"Jackson, TN"', add
label define metarea_lbl 359 `"Jacksonville, FL"', add
label define metarea_lbl 360 `"Jacksonville, NC"', add
label define metarea_lbl 361 `"Jamestown-Dunkirk, NY"', add
label define metarea_lbl 362 `"Janesville-Beloit, WI"', add
label define metarea_lbl 366 `"Johnson City-Kingsport--Bristol, TN/VA"', add
label define metarea_lbl 368 `"Johnstown, PA"', add
label define metarea_lbl 371 `"Joplin, MO"', add
label define metarea_lbl 372 `"Kalamazoo-Portage, MI"', add
label define metarea_lbl 374 `"Kankakee, IL"', add
label define metarea_lbl 376 `"Kansas City, MO-KS"', add
label define metarea_lbl 380 `"Kenosha, WI"', add
label define metarea_lbl 381 `"Kileen-Temple, TX"', add
label define metarea_lbl 384 `"Knoxville, TN"', add
label define metarea_lbl 385 `"Kokomo, IN"', add
label define metarea_lbl 387 `"LaCrosse, WI"', add
label define metarea_lbl 388 `"Lafayette, LA"', add
label define metarea_lbl 392 `"Lafayette-W. Lafayette, IN"', add
label define metarea_lbl 396 `"Lake Charles, LA"', add
label define metarea_lbl 398 `"Lakeland-Winterhaven, FL"', add
label define metarea_lbl 400 `"Lancaster, PA"', add
label define metarea_lbl 404 `"Lansing-E. Lansing, MI"', add
label define metarea_lbl 408 `"Laredo, TX"', add
label define metarea_lbl 410 `"Las Cruces, NM"', add
label define metarea_lbl 412 `"Las Vegas, NV"', add
label define metarea_lbl 415 `"Lawrence, KS"', add
label define metarea_lbl 420 `"Lawton, OK"', add
label define metarea_lbl 424 `"Lewiston-Auburn, ME"', add
label define metarea_lbl 428 `"Lexington-Fayette, KY"', add
label define metarea_lbl 432 `"Lima, OH"', add
label define metarea_lbl 436 `"Lincoln, NE"', add
label define metarea_lbl 440 `"Little Rock--North Little Rock, AR"', add
label define metarea_lbl 441 `"Long Branch-Asbury Park,NJ"', add
label define metarea_lbl 442 `"Longview-Marshall, TX"', add
label define metarea_lbl 444 `"Lorain-Elyria, OH"', add
label define metarea_lbl 448 `"Los Angeles-Long Beach, CA"', add
label define metarea_lbl 452 `"Louisville, KY/IN"', add
label define metarea_lbl 460 `"Lubbock, TX"', add
label define metarea_lbl 464 `"Lynchburg, VA"', add
label define metarea_lbl 468 `"Macon-Warner Robins, GA"', add
label define metarea_lbl 472 `"Madison, WI"', add
label define metarea_lbl 476 `"Manchester, NH"', add
label define metarea_lbl 480 `"Mansfield, OH"', add
label define metarea_lbl 484 `"Mayaguez, PR"', add
label define metarea_lbl 488 `"McAllen-Edinburg-Pharr-Mission, TX"', add
label define metarea_lbl 489 `"Medford, OR"', add
label define metarea_lbl 490 `"Melbourne-Titusville-Cocoa-Palm Bay, FL"', add
label define metarea_lbl 492 `"Memphis, TN/AR/MS"', add
label define metarea_lbl 494 `"Merced, CA"', add
label define metarea_lbl 500 `"Miami-Hialeah, FL"', add
label define metarea_lbl 504 `"Midland, TX"', add
label define metarea_lbl 508 `"Milwaukee, WI"', add
label define metarea_lbl 512 `"Minneapolis-St. Paul, MN"', add
label define metarea_lbl 514 `"Missoula, MT"', add
label define metarea_lbl 516 `"Mobile, AL"', add
label define metarea_lbl 517 `"Modesto, CA"', add
label define metarea_lbl 519 `"Monmouth-Ocean, NJ"', add
label define metarea_lbl 520 `"Monroe, LA"', add
label define metarea_lbl 524 `"Montgomery, AL"', add
label define metarea_lbl 528 `"Muncie, IN"', add
label define metarea_lbl 532 `"Muskegon-Norton Shores-Muskegon Heights, MI"', add
label define metarea_lbl 533 `"Myrtle Beach, SC"', add
label define metarea_lbl 534 `"Naples, FL"', add
label define metarea_lbl 535 `"Nashua, NH"', add
label define metarea_lbl 536 `"Nashville, TN"', add
label define metarea_lbl 540 `"New Bedford, MA"', add
label define metarea_lbl 546 `"New Brunswick-Perth Amboy-Sayreville, NJ"', add
label define metarea_lbl 548 `"New Haven-Meriden, CT"', add
label define metarea_lbl 552 `"New London-Norwich, CT/RI"', add
label define metarea_lbl 556 `"New Orleans, LA"', add
label define metarea_lbl 560 `"New York-Northeastern NJ"', add
label define metarea_lbl 564 `"Newark, OH"', add
label define metarea_lbl 566 `"Newburgh-Middletown, NY"', add
label define metarea_lbl 572 `"Norfolk-VA Beach--Newport News, VA"', add
label define metarea_lbl 576 `"Norwalk, CT"', add
label define metarea_lbl 579 `"Ocala, FL"', add
label define metarea_lbl 580 `"Odessa, TX"', add
label define metarea_lbl 588 `"Oklahoma City, OK"', add
label define metarea_lbl 591 `"Olympia, WA"', add
label define metarea_lbl 592 `"Omaha, NE/IA"', add
label define metarea_lbl 595 `"Orange, NY"', add
label define metarea_lbl 596 `"Orlando, FL"', add
label define metarea_lbl 599 `"Owensboro, KY"', add
label define metarea_lbl 601 `"Panama City, FL"', add
label define metarea_lbl 602 `"Parkersburg-Marietta,WV/OH"', add
label define metarea_lbl 603 `"Pascagoula-Moss Point, MS"', add
label define metarea_lbl 608 `"Pensacola, FL"', add
label define metarea_lbl 612 `"Peoria, IL"', add
label define metarea_lbl 616 `"Philadelphia, PA/NJ"', add
label define metarea_lbl 620 `"Phoenix, AZ"', add
label define metarea_lbl 628 `"Pittsburgh, PA"', add
label define metarea_lbl 632 `"Pittsfield, MA"', add
label define metarea_lbl 636 `"Ponse, PR"', add
label define metarea_lbl 640 `"Portland, ME"', add
label define metarea_lbl 644 `"Portland, OR-WA"', add
label define metarea_lbl 645 `"Portsmouth-Dover--Rochester, NH/ME"', add
label define metarea_lbl 646 `"Poughkeepsie, NY"', add
label define metarea_lbl 648 `"Providence-Fall River-Pawtucket, MA/RI"', add
label define metarea_lbl 652 `"Provo-Orem, UT"', add
label define metarea_lbl 656 `"Pueblo, CO"', add
label define metarea_lbl 658 `"Punta Gorda, FL"', add
label define metarea_lbl 660 `"Racine, WI"', add
label define metarea_lbl 664 `"Raleigh-Durham, NC"', add
label define metarea_lbl 666 `"Rapid City, SD"', add
label define metarea_lbl 668 `"Reading, PA"', add
label define metarea_lbl 669 `"Redding, CA"', add
label define metarea_lbl 672 `"Reno, NV"', add
label define metarea_lbl 674 `"Richland-Kennewick-Pasco, WA"', add
label define metarea_lbl 676 `"Richmond-Petersburg, VA"', add
label define metarea_lbl 678 `"Riverside-San Bernardino,CA"', add
label define metarea_lbl 680 `"Roanoke, VA"', add
label define metarea_lbl 682 `"Rochester, MN"', add
label define metarea_lbl 684 `"Rochester, NY"', add
label define metarea_lbl 688 `"Rockford, IL"', add
label define metarea_lbl 689 `"Rocky Mount, NC"', add
label define metarea_lbl 692 `"Sacramento, CA"', add
label define metarea_lbl 696 `"Saginaw-Bay City-Midland, MI"', add
label define metarea_lbl 698 `"St. Cloud, MN"', add
label define metarea_lbl 700 `"St. Joseph, MO"', add
label define metarea_lbl 704 `"St. Louis, MO-IL"', add
label define metarea_lbl 708 `"Salem, OR"', add
label define metarea_lbl 712 `"Salinas-Sea Side-Monterey, CA"', add
label define metarea_lbl 714 `"Salisbury-Concord, NC"', add
label define metarea_lbl 716 `"Salt Lake City-Ogden, UT"', add
label define metarea_lbl 720 `"San Angelo, TX"', add
label define metarea_lbl 724 `"San Antonio, TX"', add
label define metarea_lbl 732 `"San Diego, CA"', add
label define metarea_lbl 736 `"San Francisco-Oakland-Vallejo, CA"', add
label define metarea_lbl 740 `"San Jose, CA"', add
label define metarea_lbl 744 `"San Juan-Bayamon, PR"', add
label define metarea_lbl 746 `"San Luis Obispo-Atascad-P Robles, CA"', add
label define metarea_lbl 747 `"Santa Barbara-Santa Maria-Lompoc, CA"', add
label define metarea_lbl 748 `"Santa Cruz, CA"', add
label define metarea_lbl 749 `"Santa Fe, NM"', add
label define metarea_lbl 750 `"Santa Rosa-Petaluma, CA"', add
label define metarea_lbl 751 `"Sarasota, FL"', add
label define metarea_lbl 752 `"Savannah, GA"', add
label define metarea_lbl 756 `"Scranton-Wilkes-Barre, PA"', add
label define metarea_lbl 760 `"Seattle-Everett, WA"', add
label define metarea_lbl 761 `"Sharon, PA"', add
label define metarea_lbl 762 `"Sheboygan, WI"', add
label define metarea_lbl 764 `"Sherman-Davidson, TX"', add
label define metarea_lbl 768 `"Shreveport, LA"', add
label define metarea_lbl 772 `"Sioux City, IA/NE"', add
label define metarea_lbl 776 `"Sioux Falls, SD"', add
label define metarea_lbl 780 `"South Bend-Mishawaka, IN"', add
label define metarea_lbl 784 `"Spokane, WA"', add
label define metarea_lbl 788 `"Springfield, IL"', add
label define metarea_lbl 792 `"Springfield, MO"', add
label define metarea_lbl 800 `"Springfield-Holyoke-Chicopee, MA"', add
label define metarea_lbl 804 `"Stamford, CT"', add
label define metarea_lbl 805 `"State College, PA"', add
label define metarea_lbl 808 `"Steubenville-Weirton,OH/WV"', add
label define metarea_lbl 812 `"Stockton, CA"', add
label define metarea_lbl 814 `"Sumter, SC"', add
label define metarea_lbl 816 `"Syracuse, NY"', add
label define metarea_lbl 820 `"Tacoma, WA"', add
label define metarea_lbl 824 `"Tallahassee, FL"', add
label define metarea_lbl 828 `"Tampa-St. Petersburg-Clearwater, FL"', add
label define metarea_lbl 832 `"Terre Haute, IN"', add
label define metarea_lbl 836 `"Texarkana, TX/AR"', add
label define metarea_lbl 840 `"Toledo, OH/MI"', add
label define metarea_lbl 844 `"Topeka, KS"', add
label define metarea_lbl 848 `"Trenton, NJ"', add
label define metarea_lbl 852 `"Tucson, AZ"', add
label define metarea_lbl 856 `"Tulsa, OK"', add
label define metarea_lbl 860 `"Tuscaloosa, AL"', add
label define metarea_lbl 864 `"Tyler, TX"', add
label define metarea_lbl 868 `"Utica-Rome, NY"', add
label define metarea_lbl 873 `"Ventura-Oxnard-Simi Valley, CA"', add
label define metarea_lbl 875 `"Victoria, TX"', add
label define metarea_lbl 876 `"Vineland-Milville-Bridgetown, NJ"', add
label define metarea_lbl 878 `"Visalia-Tulare-Porterville, CA"', add
label define metarea_lbl 880 `"Waco, TX"', add
label define metarea_lbl 884 `"Washington, DC/MD/VA"', add
label define metarea_lbl 888 `"Waterbury, CT"', add
label define metarea_lbl 892 `"Waterloo-Cedar Falls, IA"', add
label define metarea_lbl 894 `"Wausau, WI"', add
label define metarea_lbl 896 `"West Palm Beach-Boca Raton-Delray Beach, FL"', add
label define metarea_lbl 900 `"Wheeling, WV/OH"', add
label define metarea_lbl 904 `"Wichita, KS"', add
label define metarea_lbl 908 `"Wichita Falls, TX"', add
label define metarea_lbl 914 `"Williamsport, PA"', add
label define metarea_lbl 916 `"Wilmington, DE/NJ/MD"', add
label define metarea_lbl 920 `"Wilmington, NC"', add
label define metarea_lbl 924 `"Worcester, MA"', add
label define metarea_lbl 926 `"Yakima, WA"', add
label define metarea_lbl 927 `"Yolo, CA"', add
label define metarea_lbl 928 `"York, PA"', add
label define metarea_lbl 932 `"Youngstown-Warren, OH-PA"', add
label define metarea_lbl 934 `"Yuba City, CA"', add
label define metarea_lbl 936 `"Yuma, AZ"', add
label values metarea metarea_lbl

label define metaread_lbl 0000 `"Not identifiable or not in an MSA"'
label define metaread_lbl 0040 `"Abilene, TX"', add
label define metaread_lbl 0060 `"Aguadilla, PR"', add
label define metaread_lbl 0080 `"Akron, OH"', add
label define metaread_lbl 0120 `"Albany, GA"', add
label define metaread_lbl 0160 `"Albany-Schenectady-Troy, NY"', add
label define metaread_lbl 0200 `"Albuquerque, NM"', add
label define metaread_lbl 0220 `"Alexandria, LA"', add
label define metaread_lbl 0240 `"Allentown-Bethlehem-Easton, PA/NJ"', add
label define metaread_lbl 0280 `"Altoona, PA"', add
label define metaread_lbl 0320 `"Amarillo, TX"', add
label define metaread_lbl 0380 `"Anchorage, AK"', add
label define metaread_lbl 0400 `"Anderson, IN"', add
label define metaread_lbl 0440 `"Ann Arbor, MI"', add
label define metaread_lbl 0450 `"Anniston, AL"', add
label define metaread_lbl 0460 `"Appleton-Oshkosh-Neenah, WI"', add
label define metaread_lbl 0470 `"Arecibo, PR"', add
label define metaread_lbl 0480 `"Asheville, NC"', add
label define metaread_lbl 0500 `"Athens, GA"', add
label define metaread_lbl 0520 `"Atlanta, GA"', add
label define metaread_lbl 0560 `"Atlantic City, NJ"', add
label define metaread_lbl 0580 `"Auburn-Opelika, AL"', add
label define metaread_lbl 0600 `"Augusta-Aiken, GA-SC"', add
label define metaread_lbl 0640 `"Austin, TX"', add
label define metaread_lbl 0680 `"Bakersfield, CA"', add
label define metaread_lbl 0720 `"Baltimore, MD"', add
label define metaread_lbl 0730 `"Bangor, ME"', add
label define metaread_lbl 0740 `"Barnstable-Yarmouth, MA"', add
label define metaread_lbl 0760 `"Baton Rouge, LA"', add
label define metaread_lbl 0780 `"Battle Creek, MI"', add
label define metaread_lbl 0840 `"Beaumont-Port Arthur-Orange,TX"', add
label define metaread_lbl 0860 `"Bellingham, WA"', add
label define metaread_lbl 0870 `"Benton Harbor, MI"', add
label define metaread_lbl 0880 `"Billings, MT"', add
label define metaread_lbl 0920 `"Biloxi-Gulfport, MS"', add
label define metaread_lbl 0960 `"Binghamton, NY"', add
label define metaread_lbl 1000 `"Birmingham, AL"', add
label define metaread_lbl 1010 `"Bismarck,ND"', add
label define metaread_lbl 1020 `"Bloomington, IN"', add
label define metaread_lbl 1040 `"Bloomington-Normal, IL"', add
label define metaread_lbl 1080 `"Boise City, ID"', add
label define metaread_lbl 1120 `"Boston, MA"', add
label define metaread_lbl 1121 `"Lawrence-Haverhill, MA/NH"', add
label define metaread_lbl 1122 `"Lowell, MA/NH"', add
label define metaread_lbl 1123 `"Salem-Gloucester, MA"', add
label define metaread_lbl 1140 `"Bradenton, FL"', add
label define metaread_lbl 1150 `"Bremerton, WA"', add
label define metaread_lbl 1160 `"Bridgeport, CT"', add
label define metaread_lbl 1200 `"Brockton, MA"', add
label define metaread_lbl 1240 `"Brownsville-Harlingen-San Benito, TX"', add
label define metaread_lbl 1260 `"Bryan-College Station, TX"', add
label define metaread_lbl 1280 `"Buffalo-Niagara Falls, NY"', add
label define metaread_lbl 1281 `"Niagara Falls, NY"', add
label define metaread_lbl 1300 `"Burlington, NC"', add
label define metaread_lbl 1310 `"Burlington, VT"', add
label define metaread_lbl 1320 `"Canton, OH"', add
label define metaread_lbl 1330 `"Caguas, PR"', add
label define metaread_lbl 1350 `"Casper, WY"', add
label define metaread_lbl 1360 `"Cedar Rapids, IA"', add
label define metaread_lbl 1400 `"Champaign-Urbana-Rantoul, IL"', add
label define metaread_lbl 1440 `"Charleston-N.Charleston,SC"', add
label define metaread_lbl 1480 `"Charleston, WV"', add
label define metaread_lbl 1520 `"Charlotte-Gastonia-Rock Hill, SC"', add
label define metaread_lbl 1521 `"Rock Hill, SC"', add
label define metaread_lbl 1540 `"Charlottesville, VA"', add
label define metaread_lbl 1560 `"Chattanooga, TN/GA"', add
label define metaread_lbl 1580 `"Cheyenne, WY"', add
label define metaread_lbl 1600 `"Chicago-Gary-Lake, IL"', add
label define metaread_lbl 1601 `"Aurora-Elgin, IL"', add
label define metaread_lbl 1602 `"Gary-Hammond-East Chicago, IN"', add
label define metaread_lbl 1603 `"Joliet IL"', add
label define metaread_lbl 1604 `"Lake County, IL"', add
label define metaread_lbl 1620 `"Chico, CA"', add
label define metaread_lbl 1640 `"Cincinnati OH/KY/IN"', add
label define metaread_lbl 1660 `"Clarksville-Hopkinsville, TN/KY"', add
label define metaread_lbl 1680 `"Cleveland, OH"', add
label define metaread_lbl 1720 `"Colorado Springs, CO"', add
label define metaread_lbl 1740 `"Columbia, MO"', add
label define metaread_lbl 1760 `"Columbia, SC"', add
label define metaread_lbl 1800 `"Columbus, GA/AL"', add
label define metaread_lbl 1840 `"Columbus, OH"', add
label define metaread_lbl 1880 `"Corpus Christi, TX"', add
label define metaread_lbl 1900 `"Cumberland, MD/WV"', add
label define metaread_lbl 1920 `"Dallas-Fort Worth, TX"', add
label define metaread_lbl 1921 `"Fort Worth-Arlington, TX"', add
label define metaread_lbl 1930 `"Danbury, CT"', add
label define metaread_lbl 1950 `"Danville, VA"', add
label define metaread_lbl 1960 `"Davenport, IA Rock Island-Moline, IL"', add
label define metaread_lbl 2000 `"Dayton-Springfield, OH"', add
label define metaread_lbl 2001 `"Springfield, OH"', add
label define metaread_lbl 2020 `"Daytona Beach, FL"', add
label define metaread_lbl 2030 `"Decatur, AL"', add
label define metaread_lbl 2040 `"Decatur, IL"', add
label define metaread_lbl 2080 `"Denver-Boulder-Longmont, CO"', add
label define metaread_lbl 2081 `"Boulder-Longmont, CO"', add
label define metaread_lbl 2120 `"Des Moines, IA"', add
label define metaread_lbl 2121 `"Polk, IA"', add
label define metaread_lbl 2160 `"Detroit, MI"', add
label define metaread_lbl 2180 `"Dothan, AL"', add
label define metaread_lbl 2190 `"Dover, DE"', add
label define metaread_lbl 2200 `"Dubuque, IA"', add
label define metaread_lbl 2240 `"Duluth-Superior, MN/WI"', add
label define metaread_lbl 2281 `"Dutchess Co., NY"', add
label define metaread_lbl 2290 `"Eau Claire, WI"', add
label define metaread_lbl 2310 `"El Paso, TX"', add
label define metaread_lbl 2320 `"Elkhart-Goshen, IN"', add
label define metaread_lbl 2330 `"Elmira, NY"', add
label define metaread_lbl 2340 `"Enid, OK"', add
label define metaread_lbl 2360 `"Erie, PA"', add
label define metaread_lbl 2400 `"Eugene-Springfield, OR"', add
label define metaread_lbl 2440 `"Evansville, IN/KY"', add
label define metaread_lbl 2520 `"Fargo-Morehead, ND/MN"', add
label define metaread_lbl 2560 `"Fayetteville, NC"', add
label define metaread_lbl 2580 `"Fayetteville-Springdale, AR"', add
label define metaread_lbl 2600 `"Fitchburg-Leominster, MA"', add
label define metaread_lbl 2620 `"Flagstaff, AZ-UT"', add
label define metaread_lbl 2640 `"Flint, MI"', add
label define metaread_lbl 2650 `"Florence, AL"', add
label define metaread_lbl 2660 `"Florence, SC"', add
label define metaread_lbl 2670 `"Fort Collins-Loveland, CO"', add
label define metaread_lbl 2680 `"Fort Lauderdale-Hollywood-Pompano Beach, FL"', add
label define metaread_lbl 2700 `"Fort Myers-Cape Coral, FL"', add
label define metaread_lbl 2710 `"Fort Pierce, FL"', add
label define metaread_lbl 2720 `"Fort Smith, AR/OK"', add
label define metaread_lbl 2750 `"Fort Walton Beach, FL"', add
label define metaread_lbl 2760 `"Fort Wayne, IN"', add
label define metaread_lbl 2840 `"Fresno, CA"', add
label define metaread_lbl 2880 `"Gadsden, AL"', add
label define metaread_lbl 2900 `"Gainesville, FL"', add
label define metaread_lbl 2920 `"Galveston-Texas City, TX"', add
label define metaread_lbl 2970 `"Glens Falls, NY"', add
label define metaread_lbl 2980 `"Goldsboro, NC"', add
label define metaread_lbl 2990 `"Grand Forks, ND/MN"', add
label define metaread_lbl 3000 `"Grand Rapids, MI"', add
label define metaread_lbl 3010 `"Grand Junction, CO"', add
label define metaread_lbl 3040 `"Great Falls, MT"', add
label define metaread_lbl 3060 `"Greeley, CO"', add
label define metaread_lbl 3080 `"Green Bay, WI"', add
label define metaread_lbl 3120 `"Greensboro-Winston Salem-High Point, NC"', add
label define metaread_lbl 3121 `"Winston-Salem, NC"', add
label define metaread_lbl 3150 `"Greenville, NC"', add
label define metaread_lbl 3160 `"Greenville-Spartanburg-Anderson SC"', add
label define metaread_lbl 3161 `"Anderson, SC"', add
label define metaread_lbl 3180 `"Hagerstown, MD"', add
label define metaread_lbl 3200 `"Hamilton-Middleton, OH"', add
label define metaread_lbl 3240 `"Harrisburg-Lebanon-Carlisle, PA"', add
label define metaread_lbl 3280 `"Hartford-Bristol-Middleton-New Britain, CT"', add
label define metaread_lbl 3281 `"Bristol, CT"', add
label define metaread_lbl 3282 `"Middletown, CT"', add
label define metaread_lbl 3283 `"New Britain, CT"', add
label define metaread_lbl 3290 `"Hickory-Morgantown, NC"', add
label define metaread_lbl 3300 `"Hattiesburg, MS"', add
label define metaread_lbl 3320 `"Honolulu, HI"', add
label define metaread_lbl 3350 `"Houma-Thibodoux, LA"', add
label define metaread_lbl 3360 `"Houston-Brazoria, TX"', add
label define metaread_lbl 3361 `"Brazoria, TX"', add
label define metaread_lbl 3400 `"Huntington-Ashland, WV/KY/OH"', add
label define metaread_lbl 3440 `"Huntsville, AL"', add
label define metaread_lbl 3480 `"Indianapolis, IN"', add
label define metaread_lbl 3500 `"Iowa City, IA"', add
label define metaread_lbl 3520 `"Jackson, MI"', add
label define metaread_lbl 3560 `"Jackson, MS"', add
label define metaread_lbl 3580 `"Jackson, TN"', add
label define metaread_lbl 3590 `"Jacksonville, FL"', add
label define metaread_lbl 3600 `"Jacksonville, NC"', add
label define metaread_lbl 3610 `"Jamestown-Dunkirk, NY"', add
label define metaread_lbl 3620 `"Janesville-Beloit, WI"', add
label define metaread_lbl 3660 `"Johnson City-Kingsport-Bristol, TN/VA"', add
label define metaread_lbl 3680 `"Johnstown, PA"', add
label define metaread_lbl 3710 `"Joplin, MO"', add
label define metaread_lbl 3720 `"Kalamazoo-Portage, MI"', add
label define metaread_lbl 3740 `"Kankakee, IL"', add
label define metaread_lbl 3760 `"Kansas City, MO-KS"', add
label define metaread_lbl 3800 `"Kenosha, WI"', add
label define metaread_lbl 3810 `"Kileen-Temple, TX"', add
label define metaread_lbl 3840 `"Knoxville, TN"', add
label define metaread_lbl 3850 `"Kokomo, IN"', add
label define metaread_lbl 3870 `"LaCrosse, WI"', add
label define metaread_lbl 3880 `"Lafayette, LA"', add
label define metaread_lbl 3920 `"Lafayette-W. Lafayette, IN"', add
label define metaread_lbl 3960 `"Lake Charles, LA"', add
label define metaread_lbl 3980 `"Lakeland-Winterhaven, FL"', add
label define metaread_lbl 4000 `"Lancaster, PA"', add
label define metaread_lbl 4040 `"Lansing-E. Lansing, MI"', add
label define metaread_lbl 4080 `"Laredo, TX"', add
label define metaread_lbl 4100 `"Las Cruces, NM"', add
label define metaread_lbl 4120 `"Las Vegas, NV"', add
label define metaread_lbl 4150 `"Lawrence, KS"', add
label define metaread_lbl 4200 `"Lawton, OK"', add
label define metaread_lbl 4240 `"Lewiston-Auburn, ME"', add
label define metaread_lbl 4280 `"Lexington-Fayette, KY"', add
label define metaread_lbl 4320 `"Lima, OH"', add
label define metaread_lbl 4360 `"Lincoln, NE"', add
label define metaread_lbl 4400 `"Little Rock-North Little Rock, AR"', add
label define metaread_lbl 4410 `"Long Branch-Asbury Park,NJ"', add
label define metaread_lbl 4420 `"Longview-Marshall, TX"', add
label define metaread_lbl 4440 `"Lorain-Elyria, OH"', add
label define metaread_lbl 4480 `"Los Angeles-Long Beach, CA"', add
label define metaread_lbl 4481 `"Anaheim-Santa Ana-Garden Grove, CA"', add
label define metaread_lbl 4482 `"Orange County, CA"', add
label define metaread_lbl 4520 `"Louisville, KY/IN"', add
label define metaread_lbl 4600 `"Lubbock, TX"', add
label define metaread_lbl 4640 `"Lynchburg, VA"', add
label define metaread_lbl 4680 `"Macon-Warner Robins, GA"', add
label define metaread_lbl 4720 `"Madison, WI"', add
label define metaread_lbl 4760 `"Manchester, NH"', add
label define metaread_lbl 4800 `"Mansfield, OH"', add
label define metaread_lbl 4840 `"Mayaguez, PR"', add
label define metaread_lbl 4880 `"McAllen-Edinburg-Pharr-Mission, TX"', add
label define metaread_lbl 4890 `"Medford, OR"', add
label define metaread_lbl 4900 `"Melbourne-Titusville-Cocoa-Palm Bay, FL"', add
label define metaread_lbl 4920 `"Memphis, TN/AR/MS"', add
label define metaread_lbl 4940 `"Merced, CA"', add
label define metaread_lbl 5000 `"Miami-Hialeah, FL"', add
label define metaread_lbl 5040 `"Midland, TX"', add
label define metaread_lbl 5080 `"Milwaukee, WI"', add
label define metaread_lbl 5120 `"Minneapolis-St. Paul, MN"', add
label define metaread_lbl 5140 `"Missoula, MT"', add
label define metaread_lbl 5160 `"Mobile, AL"', add
label define metaread_lbl 5170 `"Modesto, CA"', add
label define metaread_lbl 5190 `"Monmouth-Ocean, NJ"', add
label define metaread_lbl 5200 `"Monroe, LA"', add
label define metaread_lbl 5240 `"Montgomery, AL"', add
label define metaread_lbl 5280 `"Muncie, IN"', add
label define metaread_lbl 5320 `"Muskegon-Norton Shores-Muskegon Heights, MI"', add
label define metaread_lbl 5330 `"Myrtle Beach, SC"', add
label define metaread_lbl 5340 `"Naples, FL"', add
label define metaread_lbl 5350 `"Nashua, NH"', add
label define metaread_lbl 5360 `"Nashville, TN"', add
label define metaread_lbl 5400 `"New Bedford, MA"', add
label define metaread_lbl 5460 `"New Brunswick-Perth Amboy-Sayreville, NJ"', add
label define metaread_lbl 5480 `"New Haven-Meriden, CT"', add
label define metaread_lbl 5481 `"Meriden"', add
label define metaread_lbl 5482 `"New Haven, CT"', add
label define metaread_lbl 5520 `"New London-Norwich, CT/RI"', add
label define metaread_lbl 5560 `"New Orleans, LA"', add
label define metaread_lbl 5600 `"New York-Northeastern NJ"', add
label define metaread_lbl 5601 `"Nassau Co, NY"', add
label define metaread_lbl 5602 `"Bergen-Passaic, NJ"', add
label define metaread_lbl 5603 `"Jersey City, NJ"', add
label define metaread_lbl 5604 `"Middlesex-Somerset-Hunterdon, NJ"', add
label define metaread_lbl 5605 `"Newark, NJ"', add
label define metaread_lbl 5640 `"Newark, OH"', add
label define metaread_lbl 5660 `"Newburgh-Middletown, NY"', add
label define metaread_lbl 5720 `"Norfolk-VA Beach-Newport News, VA"', add
label define metaread_lbl 5721 `"Newport News-Hampton"', add
label define metaread_lbl 5722 `"Norfolk- VA Beach-Portsmouth"', add
label define metaread_lbl 5760 `"Norwalk, CT"', add
label define metaread_lbl 5790 `"Ocala, FL"', add
label define metaread_lbl 5800 `"Odessa, TX"', add
label define metaread_lbl 5880 `"Oklahoma City, OK"', add
label define metaread_lbl 5910 `"Olympia, WA"', add
label define metaread_lbl 5920 `"Omaha, NE/IA"', add
label define metaread_lbl 5950 `"Orange, NY"', add
label define metaread_lbl 5960 `"Orlando, FL"', add
label define metaread_lbl 5990 `"Owensboro, KY"', add
label define metaread_lbl 6010 `"Panama City, FL"', add
label define metaread_lbl 6020 `"Parkersburg-Marietta,WV/OH"', add
label define metaread_lbl 6030 `"Pascagoula-Moss Point, MS"', add
label define metaread_lbl 6080 `"Pensacola, FL"', add
label define metaread_lbl 6120 `"Peoria, IL"', add
label define metaread_lbl 6160 `"Philadelphia, PA/NJ"', add
label define metaread_lbl 6200 `"Phoenix, AZ"', add
label define metaread_lbl 6240 `"Pine Bluff, AR"', add
label define metaread_lbl 6280 `"Pittsburgh-Beaver Valley, PA"', add
label define metaread_lbl 6281 `"Beaver County, PA"', add
label define metaread_lbl 6320 `"Pittsfield, MA"', add
label define metaread_lbl 6360 `"Ponce, PR"', add
label define metaread_lbl 6400 `"Portland, ME"', add
label define metaread_lbl 6440 `"Portland-Vancouver, OR"', add
label define metaread_lbl 6441 `"Vancouver, WA"', add
label define metaread_lbl 6450 `"Portsmouth-Dover-Rochester, NH/ME"', add
label define metaread_lbl 6460 `"Poughkeepsie, NY"', add
label define metaread_lbl 6480 `"Providence-Fall River-Pawtucket, MA/RI"', add
label define metaread_lbl 6481 `"Fall River, MA/RI"', add
label define metaread_lbl 6482 `"Pawtuckett-Woonsocket-Attleboro, RI/MA"', add
label define metaread_lbl 6520 `"Provo-Orem, UT"', add
label define metaread_lbl 6560 `"Pueblo, CO"', add
label define metaread_lbl 6580 `"Punta Gorda, FL"', add
label define metaread_lbl 6600 `"Racine, WI"', add
label define metaread_lbl 6640 `"Raleigh-Durham, NC"', add
label define metaread_lbl 6641 `"Durham, NC"', add
label define metaread_lbl 6660 `"Rapid City, SD"', add
label define metaread_lbl 6680 `"Reading, PA"', add
label define metaread_lbl 6690 `"Redding, CA"', add
label define metaread_lbl 6720 `"Reno, NV"', add
label define metaread_lbl 6740 `"Richland-Kennewick-Pasco, WA"', add
label define metaread_lbl 6760 `"Richmond-Petersburg, VA"', add
label define metaread_lbl 6761 `"Petersburg-Colonial Heights, VA"', add
label define metaread_lbl 6780 `"Riverside-San Bernardino,CA"', add
label define metaread_lbl 6781 `"San Bernardino, CA"', add
label define metaread_lbl 6800 `"Roanoke, VA"', add
label define metaread_lbl 6820 `"Rochester, MN"', add
label define metaread_lbl 6840 `"Rochester, NY"', add
label define metaread_lbl 6880 `"Rockford, IL"', add
label define metaread_lbl 6895 `"Rocky Mount, NC"', add
label define metaread_lbl 6920 `"Sacramento, CA"', add
label define metaread_lbl 6960 `"Saginaw-Bay City-Midland, MI"', add
label define metaread_lbl 6961 `"Bay City, MI"', add
label define metaread_lbl 6980 `"St. Cloud, MN"', add
label define metaread_lbl 7000 `"St. Joseph, MO"', add
label define metaread_lbl 7040 `"St. Louis, MO-IL"', add
label define metaread_lbl 7080 `"Salem, OR"', add
label define metaread_lbl 7120 `"Salinas-Sea Side-Monterey, CA"', add
label define metaread_lbl 7140 `"Salisbury-Concord, NC"', add
label define metaread_lbl 7160 `"Salt Lake City-Ogden, UT"', add
label define metaread_lbl 7161 `"Ogden"', add
label define metaread_lbl 7200 `"San Angelo, TX"', add
label define metaread_lbl 7240 `"San Antonio, TX"', add
label define metaread_lbl 7320 `"San Diego, CA"', add
label define metaread_lbl 7360 `"San Francisco-Oakland-Vallejo, CA"', add
label define metaread_lbl 7361 `"Oakland, CA"', add
label define metaread_lbl 7362 `"Vallejo-Fairfield-Napa, CA"', add
label define metaread_lbl 7400 `"San Jose, CA"', add
label define metaread_lbl 7440 `"San Juan-Bayamon, PR"', add
label define metaread_lbl 7460 `"San Luis Obispo-Atascad-P Robles, CA"', add
label define metaread_lbl 7470 `"Santa Barbara-Santa Maria-Lompoc, CA"', add
label define metaread_lbl 7480 `"Santa Cruz, CA"', add
label define metaread_lbl 7490 `"Santa Fe, NM"', add
label define metaread_lbl 7500 `"Santa Rosa-Petaluma, CA"', add
label define metaread_lbl 7510 `"Sarasota, FL"', add
label define metaread_lbl 7520 `"Savannah, GA"', add
label define metaread_lbl 7560 `"Scranton-Wilkes-Barre, PA"', add
label define metaread_lbl 7561 `"Wilkes-Barre-Hazelton, PA"', add
label define metaread_lbl 7600 `"Seattle-Everett, WA"', add
label define metaread_lbl 7610 `"Sharon, PA"', add
label define metaread_lbl 7620 `"Sheboygan, WI"', add
label define metaread_lbl 7640 `"Sherman-Denison, TX"', add
label define metaread_lbl 7680 `"Shreveport, LA"', add
label define metaread_lbl 7720 `"Sioux City, IA/NE"', add
label define metaread_lbl 7760 `"Sioux Falls, SD"', add
label define metaread_lbl 7800 `"South Bend-Mishawaka, IN"', add
label define metaread_lbl 7840 `"Spokane, WA"', add
label define metaread_lbl 7880 `"Springfield, IL"', add
label define metaread_lbl 7920 `"Springfield, MO"', add
label define metaread_lbl 8000 `"Springfield-Holyoke-Chicopee, MA"', add
label define metaread_lbl 8040 `"Stamford, CT"', add
label define metaread_lbl 8050 `"State College, PA"', add
label define metaread_lbl 8080 `"Steubenville-Weirton,OH/WV"', add
label define metaread_lbl 8120 `"Stockton, CA"', add
label define metaread_lbl 8140 `"Sumter, SC"', add
label define metaread_lbl 8160 `"Syracuse, NY"', add
label define metaread_lbl 8200 `"Tacoma, WA"', add
label define metaread_lbl 8240 `"Tallahassee, FL"', add
label define metaread_lbl 8280 `"Tampa-St. Petersburg-Clearwater, FL"', add
label define metaread_lbl 8320 `"Terre Haute, IN"', add
label define metaread_lbl 8360 `"Texarkana, TX/AR"', add
label define metaread_lbl 8400 `"Toledo, OH/MI"', add
label define metaread_lbl 8440 `"Topeka, KS"', add
label define metaread_lbl 8480 `"Trenton, NJ"', add
label define metaread_lbl 8520 `"Tucson, AZ"', add
label define metaread_lbl 8560 `"Tulsa, OK"', add
label define metaread_lbl 8600 `"Tuscaloosa, AL"', add
label define metaread_lbl 8640 `"Tyler, TX"', add
label define metaread_lbl 8680 `"Utica-Rome, NY"', add
label define metaread_lbl 8730 `"Ventura-Oxnard-Simi Valley, CA"', add
label define metaread_lbl 8750 `"Victoria, TX"', add
label define metaread_lbl 8760 `"Vineland-Milville-Bridgetown, NJ"', add
label define metaread_lbl 8780 `"Visalia-Tulare-Porterville, CA"', add
label define metaread_lbl 8800 `"Waco, TX"', add
label define metaread_lbl 8840 `"Washington, DC/MD/VA"', add
label define metaread_lbl 8880 `"Waterbury, CT"', add
label define metaread_lbl 8920 `"Waterloo-Cedar Falls, IA"', add
label define metaread_lbl 8940 `"Wausau, WI"', add
label define metaread_lbl 8960 `"West Palm Beach-Boca Raton-Delray Beach, FL"', add
label define metaread_lbl 9000 `"Wheeling, WV/OH"', add
label define metaread_lbl 9040 `"Wichita, KS"', add
label define metaread_lbl 9080 `"Wichita Falls, TX"', add
label define metaread_lbl 9140 `"Williamsport, PA"', add
label define metaread_lbl 9160 `"Wilmington, DE/NJ/MD"', add
label define metaread_lbl 9200 `"Wilmington, NC"', add
label define metaread_lbl 9240 `"Worcester, MA"', add
label define metaread_lbl 9260 `"Yakima, WA"', add
label define metaread_lbl 9270 `"Yolo, CA"', add
label define metaread_lbl 9280 `"York, PA"', add
label define metaread_lbl 9320 `"Youngstown-Warren, OH-PA"', add
label define metaread_lbl 9340 `"Yuba City, CA"', add
label define metaread_lbl 9360 `"Yuma, AZ"', add
label values metaread metaread_lbl

label define city_lbl 0000 `"Not in identifiable city (or size group)"'
label define city_lbl 0001 `"Aberdeen, SD"', add
label define city_lbl 0002 `"Aberdeen, WA"', add
label define city_lbl 0003 `"Abilene, TX"', add
label define city_lbl 0004 `"Ada, OK"', add
label define city_lbl 0005 `"Adams, MA"', add
label define city_lbl 0006 `"Adrian, MI"', add
label define city_lbl 0010 `"Akron, OH"', add
label define city_lbl 0030 `"Alameda, CA"', add
label define city_lbl 0050 `"Albany, NY"', add
label define city_lbl 0051 `"Albany, GA"', add
label define city_lbl 0052 `"Albert Lea, MN"', add
label define city_lbl 0070 `"Albuquerque, NM"', add
label define city_lbl 0090 `"Alexandria, VA"', add
label define city_lbl 0091 `"Alexandria, LA"', add
label define city_lbl 0100 `"Alhambra, CA"', add
label define city_lbl 0101 `"Aliquippa, PA"', add
label define city_lbl 0110 `"Allegheny, PA"', add
label define city_lbl 0120 `"Aliquippa, PA"', add
label define city_lbl 0130 `"Allentown, PA"', add
label define city_lbl 0131 `"Alliance, OH"', add
label define city_lbl 0132 `"Alpena, MI"', add
label define city_lbl 0140 `"Alton, IL"', add
label define city_lbl 0150 `"Altoona, PA"', add
label define city_lbl 0160 `"Amarillo, TX"', add
label define city_lbl 0161 `"Ambridge, PA"', add
label define city_lbl 0162 `"Ames, IA"', add
label define city_lbl 0163 `"Amesbury, MA"', add
label define city_lbl 0170 `"Amsterdam, NY"', add
label define city_lbl 0171 `"Anaconda, MT"', add
label define city_lbl 0190 `"Anaheim, CA"', add
label define city_lbl 0210 `"Anchorage, AK"', add
label define city_lbl 0230 `"Anderson, IN"', add
label define city_lbl 0231 `"Anderson, SC"', add
label define city_lbl 0250 `"Andover, MA"', add
label define city_lbl 0270 `"Ann Arbor, MI"', add
label define city_lbl 0271 `"Annapolis, MD"', add
label define city_lbl 0272 `"Anniston, AL"', add
label define city_lbl 0273 `"Ansonia, CT"', add
label define city_lbl 0280 `"Appleton, WI"', add
label define city_lbl 0281 `"Ardmore, OK"', add
label define city_lbl 0282 `"Argenta, AR"', add
label define city_lbl 0283 `"Arkansas, KS"', add
label define city_lbl 0290 `"Arlington, TX"', add
label define city_lbl 0310 `"Arlington, VA"', add
label define city_lbl 0311 `"Arlington, MA"', add
label define city_lbl 0312 `"Arnold, PA"', add
label define city_lbl 0313 `"Asbury Park, NJ"', add
label define city_lbl 0330 `"Asheville, NC"', add
label define city_lbl 0331 `"Ashland, OH"', add
label define city_lbl 0340 `"Ashland, KY"', add
label define city_lbl 0341 `"Ashland, WI"', add
label define city_lbl 0342 `"Ashtabula, OH"', add
label define city_lbl 0343 `"Astoria, OR"', add
label define city_lbl 0344 `"Atchison, KS"', add
label define city_lbl 0345 `"Athens, GA"', add
label define city_lbl 0346 `"Athol, MA"', add
label define city_lbl 0350 `"Atlanta, GA"', add
label define city_lbl 0370 `"Atlantic City, NJ"', add
label define city_lbl 0371 `"Attleboro, MA"', add
label define city_lbl 0390 `"Auburn, NY"', add
label define city_lbl 0391 `"Auburn, ME"', add
label define city_lbl 0410 `"Augusta, GA"', add
label define city_lbl 0430 `"Augusta, ME"', add
label define city_lbl 0450 `"Aurora, CO"', add
label define city_lbl 0470 `"Aurora, IL"', add
label define city_lbl 0490 `"Austin, TX"', add
label define city_lbl 0491 `"Austin, MN"', add
label define city_lbl 0510 `"Bakersfield, CA"', add
label define city_lbl 0530 `"Baltimore, MD"', add
label define city_lbl 0550 `"Bangor, ME"', add
label define city_lbl 0551 `"Barberton, OH"', add
label define city_lbl 0552 `"Barre, VT"', add
label define city_lbl 0553 `"Bartlesville, OK"', add
label define city_lbl 0554 `"Batavia, NY"', add
label define city_lbl 0570 `"Bath, ME"', add
label define city_lbl 0590 `"Baton Rouge, LA"', add
label define city_lbl 0610 `"Battle Creek, MI"', add
label define city_lbl 0630 `"Bay City, MI"', add
label define city_lbl 0640 `"Bayamon, PR"', add
label define city_lbl 0650 `"Bayonne, NJ"', add
label define city_lbl 0651 `"Beacon, NY"', add
label define city_lbl 0652 `"Beatrice, NE"', add
label define city_lbl 0660 `"Belleville, IL"', add
label define city_lbl 0670 `"Beaumont, TX"', add
label define city_lbl 0671 `"Beaver Falls, PA"', add
label define city_lbl 0672 `"Bedford, IN"', add
label define city_lbl 0673 `"Bellaire, OH"', add
label define city_lbl 0680 `"Bellevue, WA"', add
label define city_lbl 0690 `"Bellingham, WA"', add
label define city_lbl 0700 `"Belleville, NJ"', add
label define city_lbl 0701 `"Bellevue, PA"', add
label define city_lbl 0702 `"Belmont, OH"', add
label define city_lbl 0703 `"Belmont, MA"', add
label define city_lbl 0704 `"Beloit, WI"', add
label define city_lbl 0705 `"Bennington, VT"', add
label define city_lbl 0706 `"Benton Harbor, MI"', add
label define city_lbl 0710 `"Berkeley, CA"', add
label define city_lbl 0711 `"Berlin, NH"', add
label define city_lbl 0712 `"Berwick, PA"', add
label define city_lbl 0720 `"Berwyn, IL"', add
label define city_lbl 0721 `"Bessemer, AL"', add
label define city_lbl 0730 `"Bethlehem, PA"', add
label define city_lbl 0740 `"Biddeford, ME"', add
label define city_lbl 0741 `"Big Spring, TX"', add
label define city_lbl 0742 `"Billings, MT"', add
label define city_lbl 0743 `"Biloxi, MS"', add
label define city_lbl 0750 `"Binghamton, NY"', add
label define city_lbl 0760 `"Beverly, MA"', add
label define city_lbl 0761 `"Beverly Hills, CA"', add
label define city_lbl 0770 `"Birmingham, AL"', add
label define city_lbl 0771 `"Birmingham, CT"', add
label define city_lbl 0772 `"Bismarck, ND"', add
label define city_lbl 0780 `"Bloomfield, NJ"', add
label define city_lbl 0790 `"Bloomington, IL"', add
label define city_lbl 0791 `"Bloomington, IN"', add
label define city_lbl 0792 `"Blue Island, IL"', add
label define city_lbl 0793 `"Bluefield, WV"', add
label define city_lbl 0794 `"Blytheville, AR"', add
label define city_lbl 0795 `"Bogalusa, LA"', add
label define city_lbl 0800 `"Boise, ID"', add
label define city_lbl 0801 `"Boone, IA"', add
label define city_lbl 0810 `"Boston, MA"', add
label define city_lbl 0811 `"Boulder, CO"', add
label define city_lbl 0812 `"Bowling Green, KY"', add
label define city_lbl 0813 `"Braddock, PA"', add
label define city_lbl 0814 `"Braden, WA"', add
label define city_lbl 0815 `"Bradford, PA"', add
label define city_lbl 0816 `"Brainerd, MN"', add
label define city_lbl 0817 `"Braintree, MA"', add
label define city_lbl 0818 `"Brawley, CA"', add
label define city_lbl 0819 `"Bremerton, WA"', add
label define city_lbl 0830 `"Bridgeport, CT"', add
label define city_lbl 0831 `"Bridgeton, NJ"', add
label define city_lbl 0832 `"Bristol, CT"', add
label define city_lbl 0833 `"Bristol, PA"', add
label define city_lbl 0834 `"Bristol, VA"', add
label define city_lbl 0835 `"Bristol, TN"', add
label define city_lbl 0837 `"Bristol, RI"', add
label define city_lbl 0850 `"Brockton, MA"', add
label define city_lbl 0851 `"Brookfield, IL"', add
label define city_lbl 0870 `"Brookline, MA"', add
label define city_lbl 0880 `"Brownsville, TX"', add
label define city_lbl 0881 `"Brownwood, TX"', add
label define city_lbl 0882 `"Brunswick, GA"', add
label define city_lbl 0883 `"Bucyrus, OH"', add
label define city_lbl 0890 `"Buffalo, NY"', add
label define city_lbl 0900 `"Burlington, IA"', add
label define city_lbl 0905 `"Burlington, VT"', add
label define city_lbl 0906 `"Burlington, NJ"', add
label define city_lbl 0907 `"Bushkill, PA"', add
label define city_lbl 0910 `"Butte, MT"', add
label define city_lbl 0911 `"Butler, PA"', add
label define city_lbl 0920 `"Burbank, CA"', add
label define city_lbl 0921 `"Burlingame, CA"', add
label define city_lbl 0926 `"Cairo, IL"', add
label define city_lbl 0927 `"Calumet City, IL"', add
label define city_lbl 0930 `"Cambridge, MA"', add
label define city_lbl 0931 `"Cambridge, OH"', add
label define city_lbl 0950 `"Camden, NJ"', add
label define city_lbl 0951 `"Campbell, OH"', add
label define city_lbl 0952 `"Canonsburg, PA"', add
label define city_lbl 0970 `"Camden, NY"', add
label define city_lbl 0990 `"Canton, OH"', add
label define city_lbl 0991 `"Canton, IL"', add
label define city_lbl 0992 `"Cape Girardeau, MO"', add
label define city_lbl 0993 `"Carbondale, PA"', add
label define city_lbl 0994 `"Carlisle, PA"', add
label define city_lbl 0995 `"Carnegie, PA"', add
label define city_lbl 0996 `"Carrick, PA"', add
label define city_lbl 0997 `"Carteret, NJ"', add
label define city_lbl 0998 `"Carthage, MO"', add
label define city_lbl 0999 `"Casper, WY"', add
label define city_lbl 1000 `"Cape Coral, FL"', add
label define city_lbl 1010 `"Cedar Rapids, IA"', add
label define city_lbl 1020 `"Central Falls, RI"', add
label define city_lbl 1021 `"Centralia, IL"', add
label define city_lbl 1023 `"Chambersburg, PA"', add
label define city_lbl 1024 `"Champaign, IL"', add
label define city_lbl 1025 `"Chanute, KS"', add
label define city_lbl 1026 `"Charleroi, PA"', add
label define city_lbl 1030 `"Charlestown, MA"', add
label define city_lbl 1050 `"Charleston, SC"', add
label define city_lbl 1060 `"Carolina, PR"', add
label define city_lbl 1070 `"Charleston, WV"', add
label define city_lbl 1090 `"Charlotte, NC"', add
label define city_lbl 1091 `"Charlottesville, VA"', add
label define city_lbl 1110 `"Chattanooga, TN"', add
label define city_lbl 1130 `"Chelsea, MA"', add
label define city_lbl 1150 `"Chesapeake, VA"', add
label define city_lbl 1170 `"Chester, PA"', add
label define city_lbl 1171 `"Cheyenne, WY"', add
label define city_lbl 1190 `"Chicago, IL"', add
label define city_lbl 1191 `"Chicago Heights, IL"', add
label define city_lbl 1192 `"Chickasha, OK"', add
label define city_lbl 1210 `"Chicopee, MA"', add
label define city_lbl 1230 `"Chillicothe, OH"', add
label define city_lbl 1250 `"Chula Vista, CA"', add
label define city_lbl 1270 `"Cicero, IL"', add
label define city_lbl 1290 `"Cincinnati, OH"', add
label define city_lbl 1291 `"Clairton, PA"', add
label define city_lbl 1292 `"Claremont, NH"', add
label define city_lbl 1310 `"Clarksburg, WV"', add
label define city_lbl 1311 `"Clarksdale, MS"', add
label define city_lbl 1312 `"Cleburne, TX"', add
label define city_lbl 1330 `"Cleveland, OH"', add
label define city_lbl 1340 `"Cleveland Heights, OH"', add
label define city_lbl 1341 `"Cliffside Park, NJ"', add
label define city_lbl 1350 `"Clifton, NJ"', add
label define city_lbl 1351 `"Clinton, IN"', add
label define city_lbl 1370 `"Clinton, IA"', add
label define city_lbl 1371 `"Clinton, MA"', add
label define city_lbl 1372 `"Coatesville, PA"', add
label define city_lbl 1373 `"Coffeyville, KS"', add
label define city_lbl 1374 `"Cohoes, NY"', add
label define city_lbl 1375 `"Collingswood, NJ"', add
label define city_lbl 1390 `"Colorado Springs, CO"', add
label define city_lbl 1400 `"Cohoes, NY"', add
label define city_lbl 1410 `"Columbia, SC"', add
label define city_lbl 1411 `"Columbia, PA"', add
label define city_lbl 1412 `"Columbia, MO"', add
label define city_lbl 1420 `"Columbia City, IN"', add
label define city_lbl 1430 `"Columbus, GA"', add
label define city_lbl 1450 `"Columbus, OH"', add
label define city_lbl 1451 `"Columbus, MS"', add
label define city_lbl 1452 `"Compton, CA"', add
label define city_lbl 1470 `"Concord, CA"', add
label define city_lbl 1490 `"Concord, NH"', add
label define city_lbl 1491 `"Concord, NC"', add
label define city_lbl 1492 `"Connellsville, PA"', add
label define city_lbl 1493 `"Connersville, IN"', add
label define city_lbl 1494 `"Conshohocken, PA"', add
label define city_lbl 1495 `"Coraopolis, PA"', add
label define city_lbl 1496 `"Corning, NY"', add
label define city_lbl 1500 `"Corona, CA"', add
label define city_lbl 1510 `"Council Bluffs, IA"', add
label define city_lbl 1520 `"Corpus Christi, TX"', add
label define city_lbl 1521 `"Corsicana, TX"', add
label define city_lbl 1522 `"Cortland, NY"', add
label define city_lbl 1523 `"Coshocton, OH"', add
label define city_lbl 1530 `"Covington, KY"', add
label define city_lbl 1540 `"Costa Mesa, CA"', add
label define city_lbl 1550 `"Cranston, RI"', add
label define city_lbl 1551 `"Crawfordsville, IN"', add
label define city_lbl 1552 `"Cripple Creek, CO"', add
label define city_lbl 1553 `"Cudahy, WI"', add
label define city_lbl 1570 `"Cumberland, MD"', add
label define city_lbl 1571 `"Cumberland, RI"', add
label define city_lbl 1572 `"Cuyahoga Falls, OH"', add
label define city_lbl 1590 `"Dallas, TX"', add
label define city_lbl 1591 `"Danbury, CT"', add
label define city_lbl 1610 `"Danvers, MA"', add
label define city_lbl 1630 `"Danville, IL"', add
label define city_lbl 1631 `"Danville, VA"', add
label define city_lbl 1650 `"Davenport, IA"', add
label define city_lbl 1670 `"Dayton, OH"', add
label define city_lbl 1671 `"Daytona Beach, FL"', add
label define city_lbl 1680 `"Dearborn, MI"', add
label define city_lbl 1690 `"Decatur, IL"', add
label define city_lbl 1691 `"Decatur, AL"', add
label define city_lbl 1692 `"Decatur, GA"', add
label define city_lbl 1693 `"Dedham, MA"', add
label define city_lbl 1694 `"Del Rio, TX"', add
label define city_lbl 1695 `"Denison, TX"', add
label define city_lbl 1710 `"Denver, CO"', add
label define city_lbl 1711 `"Derby, CT"', add
label define city_lbl 1713 `"Derry, PA"', add
label define city_lbl 1730 `"Des Moines, IA"', add
label define city_lbl 1750 `"Detroit, MI"', add
label define city_lbl 1751 `"Dickson City, PA"', add
label define city_lbl 1752 `"Dodge, KS"', add
label define city_lbl 1753 `"Donora, PA"', add
label define city_lbl 1754 `"Dormont, PA"', add
label define city_lbl 1755 `"Dothan, AL"', add
label define city_lbl 1770 `"Dorchester, MA"', add
label define city_lbl 1790 `"Dover, NH"', add
label define city_lbl 1791 `"Dover, NJ"', add
label define city_lbl 1792 `"Du Bois, PA"', add
label define city_lbl 1800 `"Downey, CA"', add
label define city_lbl 1810 `"Dubuque, IA"', add
label define city_lbl 1830 `"Duluth, MN"', add
label define city_lbl 1831 `"Dunkirk, NY"', add
label define city_lbl 1832 `"Dunmore, PA"', add
label define city_lbl 1833 `"Duquesne, PA"', add
label define city_lbl 1850 `"Durham, NC"', add
label define city_lbl 1860 `"1860"', add
label define city_lbl 1870 `"East Chicago, IN"', add
label define city_lbl 1890 `"East Cleveland, OH"', add
label define city_lbl 1891 `"East Hartford, CT"', add
label define city_lbl 1892 `"East Liverpool, OH"', add
label define city_lbl 1893 `"East Moline, IL"', add
label define city_lbl 1910 `"East Los Angeles, CA"', add
label define city_lbl 1930 `"East Orange, NJ"', add
label define city_lbl 1931 `"East Providence, RI"', add
label define city_lbl 1940 `"East Saginaw, MI"', add
label define city_lbl 1950 `"East St. Louis, IL"', add
label define city_lbl 1951 `"East Youngstown, OH"', add
label define city_lbl 1952 `"Easthampton, MA"', add
label define city_lbl 1970 `"Easton, PA"', add
label define city_lbl 1971 `"Eau Claire, WI"', add
label define city_lbl 1972 `"Ecorse, MI"', add
label define city_lbl 1973 `"El Dorado, KS"', add
label define city_lbl 1974 `"El Dorado, AR"', add
label define city_lbl 1990 `"El Monte, CA"', add
label define city_lbl 2010 `"El Paso, TX"', add
label define city_lbl 2030 `"Elgin, IL"', add
label define city_lbl 2040 `"Elyria, OH"', add
label define city_lbl 2050 `"Elizabeth, NJ"', add
label define city_lbl 2051 `"Elizabeth City, NC"', add
label define city_lbl 2060 `"Elkhart, IN"', add
label define city_lbl 2061 `"Ellwood City, PA"', add
label define city_lbl 2062 `"Elmhurst, IL"', add
label define city_lbl 2070 `"Elmira, NY"', add
label define city_lbl 2071 `"Elmwood Park, IL"', add
label define city_lbl 2072 `"Elwood, IN"', add
label define city_lbl 2073 `"Emporia, KS"', add
label define city_lbl 2074 `"Endicott, NY"', add
label define city_lbl 2075 `"Enfield, CT"', add
label define city_lbl 2076 `"Englewood, NJ"', add
label define city_lbl 2080 `"Enid, OK"', add
label define city_lbl 2090 `"Erie, PA"', add
label define city_lbl 2091 `"Escanaba, MI"', add
label define city_lbl 2092 `"Euclid, OH"', add
label define city_lbl 2110 `"Escondido, CA"', add
label define city_lbl 2130 `"Eugene, OR"', add
label define city_lbl 2131 `"Eureka, CA"', add
label define city_lbl 2150 `"Evanston, IL"', add
label define city_lbl 2170 `"Evansville, IN"', add
label define city_lbl 2190 `"Everett, MA"', add
label define city_lbl 2210 `"Everett, WA"', add
label define city_lbl 2211 `"Fairfield, AL"', add
label define city_lbl 2212 `"Fairfield, CT"', add
label define city_lbl 2213 `"Fairhaven, MA"', add
label define city_lbl 2214 `"Fairmont, WV"', add
label define city_lbl 2220 `"Fargo, ND"', add
label define city_lbl 2221 `"Faribault, MN"', add
label define city_lbl 2222 `"Farrell, PA"', add
label define city_lbl 2230 `"Fall River, MA"', add
label define city_lbl 2240 `"Fayetteville, NC"', add
label define city_lbl 2241 `"Ferndale, MI"', add
label define city_lbl 2242 `"Findlay, OH"', add
label define city_lbl 2250 `"Fitchburg, MA"', add
label define city_lbl 2260 `"Fontana, CA"', add
label define city_lbl 2270 `"Flint, MI"', add
label define city_lbl 2271 `"Floral Park, NY"', add
label define city_lbl 2273 `"Florence, AL"', add
label define city_lbl 2274 `"Florence, SC"', add
label define city_lbl 2275 `"Flushing, NY"', add
label define city_lbl 2280 `"Fond du Lac, WI"', add
label define city_lbl 2281 `"Forest Park, IL"', add
label define city_lbl 2290 `"Fort Lauderdale, FL"', add
label define city_lbl 2300 `"Fort Collins, CO"', add
label define city_lbl 2301 `"Fort Dodge, IA"', add
label define city_lbl 2302 `"Fort Madison, IA"', add
label define city_lbl 2303 `"Fort Scott, KS"', add
label define city_lbl 2310 `"Fort Smith, AR"', add
label define city_lbl 2311 `"Fort Thomas, KY"', add
label define city_lbl 2330 `"Fort Wayne, IN"', add
label define city_lbl 2350 `"Fort Worth, TX"', add
label define city_lbl 2351 `"Fostoria, OH"', add
label define city_lbl 2352 `"Framingham, MA"', add
label define city_lbl 2353 `"Frankfort, IN"', add
label define city_lbl 2354 `"Frankfort, KY"', add
label define city_lbl 2355 `"Franklin, PA"', add
label define city_lbl 2356 `"Frederick, MD"', add
label define city_lbl 2357 `"Freeport, NY"', add
label define city_lbl 2358 `"Freeport, IL"', add
label define city_lbl 2359 `"Fremont, OH"', add
label define city_lbl 2360 `"Fremont, NE"', add
label define city_lbl 2370 `"Fresno, CA"', add
label define city_lbl 2390 `"Fullerton, CA"', add
label define city_lbl 2391 `"Fulton, NY"', add
label define city_lbl 2392 `"Gadsden, AL"', add
label define city_lbl 2393 `"Galena, KS"', add
label define city_lbl 2400 `"Galesburg, IL"', add
label define city_lbl 2410 `"Galveston, TX"', add
label define city_lbl 2411 `"Gardner, MA"', add
label define city_lbl 2430 `"Garden Grove, CA"', add
label define city_lbl 2440 `"Garfield, NJ"', add
label define city_lbl 2441 `"Garfield Heights, OH"', add
label define city_lbl 2450 `"Garland, TX"', add
label define city_lbl 2470 `"Gary, IN"', add
label define city_lbl 2471 `"Gastonia, NC"', add
label define city_lbl 2472 `"Geneva, NY"', add
label define city_lbl 2473 `"Glen Cove, NY"', add
label define city_lbl 2490 `"Glendale, CA"', add
label define city_lbl 2491 `"Glens Falls, NY"', add
label define city_lbl 2510 `"Gloucester, MA"', add
label define city_lbl 2511 `"Gloucester, NJ"', add
label define city_lbl 2512 `"Gloversville, NY"', add
label define city_lbl 2513 `"Goldsboro, NC"', add
label define city_lbl 2514 `"Goshen, IN"', add
label define city_lbl 2515 `"Grand Forks, ND"', add
label define city_lbl 2516 `"Grand Island, NE"', add
label define city_lbl 2517 `"Grand Junction, CO"', add
label define city_lbl 2520 `"Granite City, IL"', add
label define city_lbl 2530 `"Grand Rapids, MI"', add
label define city_lbl 2531 `"Grandville, MI"', add
label define city_lbl 2540 `"Great Falls, MT"', add
label define city_lbl 2541 `"Greeley, CO"', add
label define city_lbl 2550 `"Green Bay, WI"', add
label define city_lbl 2551 `"Greenfield, MA"', add
label define city_lbl 2570 `"Greensboro, NC"', add
label define city_lbl 2571 `"Greensburg, PA"', add
label define city_lbl 2572 `"Greenville, MS"', add
label define city_lbl 2573 `"Greenville, SC"', add
label define city_lbl 2574 `"Greenville, TX"', add
label define city_lbl 2575 `"Greenwich, CT"', add
label define city_lbl 2576 `"Greenwood, MS"', add
label define city_lbl 2577 `"Greenwood, SC"', add
label define city_lbl 2578 `"Griffin, GA"', add
label define city_lbl 2579 `"Grosse Pointe Park, MI"', add
label define city_lbl 2580 `"Guynabo, PR"', add
label define city_lbl 2581 `"Groton, CT"', add
label define city_lbl 2582 `"Gulfport, MS"', add
label define city_lbl 2583 `"Guthrie, OK"', add
label define city_lbl 2584 `"Hackensack, NJ"', add
label define city_lbl 2590 `"Hagerstown, MD"', add
label define city_lbl 2591 `"Hamden, CT"', add
label define city_lbl 2610 `"Hamilton, OH"', add
label define city_lbl 2630 `"Hammond, IN"', add
label define city_lbl 2650 `"Hampton, VA"', add
label define city_lbl 2670 `"Hamtramck Village, MI"', add
label define city_lbl 2680 `"Hannibal, MO"', add
label define city_lbl 2681 `"Hanover, PA"', add
label define city_lbl 2682 `"Harlingen, TX"', add
label define city_lbl 2690 `"Harrisburg, PA"', add
label define city_lbl 2691 `"Harrisburg, IL"', add
label define city_lbl 2692 `"Harrison, NJ"', add
label define city_lbl 2710 `"Hartford, CT"', add
label define city_lbl 2711 `"Harvey, IL"', add
label define city_lbl 2712 `"Hastings, NE"', add
label define city_lbl 2713 `"Hattiesburg, MS"', add
label define city_lbl 2730 `"Haverhill, MA"', add
label define city_lbl 2731 `"Hawthorne, NJ"', add
label define city_lbl 2750 `"Hazleton, PA"', add
label define city_lbl 2751 `"Helena, MT"', add
label define city_lbl 2752 `"Hempstead, NY"', add
label define city_lbl 2753 `"Henderson, KY"', add
label define city_lbl 2754 `"Herkimer, NY"', add
label define city_lbl 2755 `"Herrin, IL"', add
label define city_lbl 2756 `"Hibbing, MN"', add
label define city_lbl 2770 `"Hialeah, FL"', add
label define city_lbl 2780 `"High Point, NC"', add
label define city_lbl 2781 `"Highland Park, IL"', add
label define city_lbl 2790 `"Highland Park, MI"', add
label define city_lbl 2791 `"Hilo, HI"', add
label define city_lbl 2810 `"Hoboken, NJ"', add
label define city_lbl 2811 `"Holland, MI"', add
label define city_lbl 2830 `"Hollywood, FL"', add
label define city_lbl 2850 `"Holyoke, MA"', add
label define city_lbl 2851 `"Homestead, PA"', add
label define city_lbl 2870 `"Honolulu, HI"', add
label define city_lbl 2871 `"Hopewell, VA"', add
label define city_lbl 2872 `"Hopkinsville, KY"', add
label define city_lbl 2873 `"Hoquiam, WA"', add
label define city_lbl 2874 `"Hornell, NY"', add
label define city_lbl 2875 `"Hot Springs, AR"', add
label define city_lbl 2890 `"Houston, TX"', add
label define city_lbl 2891 `"Hudson, NY"', add
label define city_lbl 2892 `"Huntington, IN"', add
label define city_lbl 2910 `"Huntington, WV"', add
label define city_lbl 2930 `"Huntington Beach, CA"', add
label define city_lbl 2950 `"Huntsville, AL"', add
label define city_lbl 2951 `"Huron, SD"', add
label define city_lbl 2960 `"Hutchinson, KS"', add
label define city_lbl 2961 `"Hyde Park, MA"', add
label define city_lbl 2962 `"Ilion, NY"', add
label define city_lbl 2963 `"Independence, KS"', add
label define city_lbl 2970 `"Independence, MO"', add
label define city_lbl 2990 `"Indianapolis, IN"', add
label define city_lbl 3010 `"Inglewood, CA"', add
label define city_lbl 3011 `"Iowa City, IA"', add
label define city_lbl 3012 `"Iron Mountain, MI"', add
label define city_lbl 3013 `"Ironton, OH"', add
label define city_lbl 3014 `"Ironwood, MI"', add
label define city_lbl 3020 `"Irvine, CA"', add
label define city_lbl 3030 `"Irving, TX"', add
label define city_lbl 3050 `"Irvington, NJ"', add
label define city_lbl 3051 `"Ishpeming, MI"', add
label define city_lbl 3052 `"Ithaca, NY"', add
label define city_lbl 3070 `"Jackson, MI"', add
label define city_lbl 3071 `"Jackson, MN"', add
label define city_lbl 3090 `"Jackson, MS"', add
label define city_lbl 3091 `"Jackson, TN"', add
label define city_lbl 3110 `"Jacksonville, FL"', add
label define city_lbl 3111 `"Jacksonville, IL"', add
label define city_lbl 3130 `"Jamestown , NY"', add
label define city_lbl 3131 `"Janesville, WI"', add
label define city_lbl 3132 `"Jeannette, PA"', add
label define city_lbl 3133 `"Jefferson City, MO"', add
label define city_lbl 3134 `"Jeffersonville, IN"', add
label define city_lbl 3150 `"Jersey City, NJ"', add
label define city_lbl 3151 `"Johnson City, NY"', add
label define city_lbl 3160 `"Johnson City, TN"', add
label define city_lbl 3161 `"Johnstown, NY"', add
label define city_lbl 3170 `"Johnstown, PA"', add
label define city_lbl 3190 `"Joliet, IL"', add
label define city_lbl 3191 `"Jonesboro, AR"', add
label define city_lbl 3210 `"Joplin, MO"', add
label define city_lbl 3230 `"Kalamazoo, MI"', add
label define city_lbl 3231 `"Kankakee, IL"', add
label define city_lbl 3250 `"Kansas City, KS"', add
label define city_lbl 3260 `"Kansas City, MO"', add
label define city_lbl 3270 `"Kearny, NJ"', add
label define city_lbl 3271 `"Keene, NH"', add
label define city_lbl 3272 `"Kenmore, NY"', add
label define city_lbl 3273 `"Kenmore, OH"', add
label define city_lbl 3290 `"Kenosha, WI"', add
label define city_lbl 3291 `"Keokuk, IA"', add
label define city_lbl 3292 `"Kewanee, IL"', add
label define city_lbl 3293 `"Key West, FL"', add
label define city_lbl 3294 `"Kingsport, TN"', add
label define city_lbl 3310 `"Kingston, NY"', add
label define city_lbl 3311 `"Kingston, PA"', add
label define city_lbl 3312 `"Kinston, NC"', add
label define city_lbl 3313 `"Klamath Falls, OR"', add
label define city_lbl 3330 `"Knoxville, TN"', add
label define city_lbl 3350 `"Kokomo, IN"', add
label define city_lbl 3370 `"La Crosse, WI"', add
label define city_lbl 3380 `"Lafayette, IL"', add
label define city_lbl 3390 `"Lafayette, LA"', add
label define city_lbl 3391 `"La Grange, IL"', add
label define city_lbl 3392 `"La Grange, GA"', add
label define city_lbl 3393 `"La Porte, IN"', add
label define city_lbl 3394 `"La Salle, IL"', add
label define city_lbl 3395 `"Lackawanna, NY"', add
label define city_lbl 3396 `"Laconia, NH"', add
label define city_lbl 3410 `"Lakewood, CO"', add
label define city_lbl 3430 `"Lakewood, OH"', add
label define city_lbl 3440 `"Lancaster, CA"', add
label define city_lbl 3450 `"Lancaster, PA"', add
label define city_lbl 3451 `"Lancaster, OH"', add
label define city_lbl 3470 `"Lansing, MI"', add
label define city_lbl 3471 `"Lansingburgh, NY"', add
label define city_lbl 3480 `"Laredo, TX"', add
label define city_lbl 3481 `"Latrobe, PA"', add
label define city_lbl 3482 `"Laurel, MS"', add
label define city_lbl 3490 `"Las Vegas, NV"', add
label define city_lbl 3510 `"Lawrence, MA"', add
label define city_lbl 3511 `"Lawrence, KS"', add
label define city_lbl 3512 `"Lawton, OK"', add
label define city_lbl 3513 `"Leadville, CO"', add
label define city_lbl 3520 `"Leavenworth, KS"', add
label define city_lbl 3521 `"Lebanon, PA"', add
label define city_lbl 3522 `"Leominster, MA"', add
label define city_lbl 3530 `"Lehigh, PA"', add
label define city_lbl 3540 `"Lebanon, PA"', add
label define city_lbl 3550 `"Lewiston, ME"', add
label define city_lbl 3551 `"Lewistown, PA"', add
label define city_lbl 3570 `"Lexington, KY"', add
label define city_lbl 3590 `"Lexington-Fayette, KY"', add
label define city_lbl 3610 `"Lima, OH"', add
label define city_lbl 3630 `"Lincoln, NE"', add
label define city_lbl 3631 `"Lincoln, IL"', add
label define city_lbl 3632 `"Lincoln Park, MI"', add
label define city_lbl 3633 `"Lincoln, RI"', add
label define city_lbl 3634 `"Linden, NJ"', add
label define city_lbl 3635 `"Little Falls, NY"', add
label define city_lbl 3638 `"Lodi, NJ"', add
label define city_lbl 3639 `"Logansport, IN"', add
label define city_lbl 3650 `"Little Rock, AR"', add
label define city_lbl 3670 `"Livonia, MI"', add
label define city_lbl 3680 `"Lockport, NY"', add
label define city_lbl 3690 `"Long Beach, CA"', add
label define city_lbl 3691 `"Long Branch, NJ"', add
label define city_lbl 3692 `"Long Island City, NY"', add
label define city_lbl 3693 `"Longview, WA"', add
label define city_lbl 3710 `"Lorain, OH"', add
label define city_lbl 3730 `"Los Angeles, CA"', add
label define city_lbl 3750 `"Louisville, KY"', add
label define city_lbl 3770 `"Lowell, MA"', add
label define city_lbl 3771 `"Lubbock, TX"', add
label define city_lbl 3772 `"Lynbrook, NY"', add
label define city_lbl 3790 `"Lynchburg, VA"', add
label define city_lbl 3810 `"Lynn, MA"', add
label define city_lbl 3830 `"Macon, GA"', add
label define city_lbl 3850 `"Madison, IN"', add
label define city_lbl 3870 `"Madison, WI"', add
label define city_lbl 3871 `"Mahanoy City, PA"', add
label define city_lbl 3890 `"Malden, MA"', add
label define city_lbl 3891 `"Mamaroneck, NY"', add
label define city_lbl 3910 `"Manchester, NH"', add
label define city_lbl 3911 `"Manchester, CT"', add
label define city_lbl 3912 `"Manhattan, KS"', add
label define city_lbl 3913 `"Manistee, MI"', add
label define city_lbl 3914 `"Manitowoc, WI"', add
label define city_lbl 3915 `"Mankato, MN"', add
label define city_lbl 3930 `"Mansfield, OH"', add
label define city_lbl 3931 `"Maplewood, MO"', add
label define city_lbl 3932 `"Marietta, OH"', add
label define city_lbl 3933 `"Marinette, WI"', add
label define city_lbl 3934 `"Marion, IN"', add
label define city_lbl 3940 `"Maywood, IL"', add
label define city_lbl 3950 `"Marion, OH"', add
label define city_lbl 3951 `"Marlborough, MA"', add
label define city_lbl 3952 `"Marquette, MI"', add
label define city_lbl 3953 `"Marshall, TX"', add
label define city_lbl 3954 `"Marshalltown, IA"', add
label define city_lbl 3955 `"Martins Ferry, OH"', add
label define city_lbl 3956 `"Martinsburg, WV"', add
label define city_lbl 3957 `"Mason City, IA"', add
label define city_lbl 3958 `"Massena, NY"', add
label define city_lbl 3959 `"Massillon, OH"', add
label define city_lbl 3960 `"McAllen, TX"', add
label define city_lbl 3961 `"Mattoon, IL"', add
label define city_lbl 3962 `"Mcalester, OK"', add
label define city_lbl 3963 `"Mccomb, MS"', add
label define city_lbl 3964 `"Mckees Rocks, PA"', add
label define city_lbl 3970 `"McKeesport, PA"', add
label define city_lbl 3971 `"Meadville, PA"', add
label define city_lbl 3990 `"Medford, MA"', add
label define city_lbl 3991 `"Medford, OR"', add
label define city_lbl 3992 `"Melrose, MA"', add
label define city_lbl 3993 `"Melrose Park, IL"', add
label define city_lbl 4010 `"Memphis, TN"', add
label define city_lbl 4011 `"Menominee, MI"', add
label define city_lbl 4030 `"Meriden, CT"', add
label define city_lbl 4040 `"Meridian, MS"', add
label define city_lbl 4041 `"Methuen, MA"', add
label define city_lbl 4050 `"Mesa, AZ"', add
label define city_lbl 4070 `"Mesquite, TX"', add
label define city_lbl 4090 `"Metairie, LA"', add
label define city_lbl 4110 `"Miami, FL"', add
label define city_lbl 4120 `"Michigan City, IN"', add
label define city_lbl 4121 `"Middlesborough, KY"', add
label define city_lbl 4122 `"Middletown, CT"', add
label define city_lbl 4123 `"Middletown, NY"', add
label define city_lbl 4124 `"Middletown, OH"', add
label define city_lbl 4125 `"Milford, CT"', add
label define city_lbl 4126 `"Milford, MA"', add
label define city_lbl 4127 `"Millville, NJ"', add
label define city_lbl 4128 `"Milton, MA"', add
label define city_lbl 4130 `"Milwaukee, WI"', add
label define city_lbl 4150 `"Minneapolis, MN"', add
label define city_lbl 4151 `"Minot, ND"', add
label define city_lbl 4160 `"Mishawaka, IN"', add
label define city_lbl 4161 `"Missoula, MT"', add
label define city_lbl 4162 `"Mitchell, SD"', add
label define city_lbl 4163 `"Moberly, MO"', add
label define city_lbl 4170 `"Mobile, AL"', add
label define city_lbl 4190 `"Modesto, CA"', add
label define city_lbl 4210 `"Moline, IL"', add
label define city_lbl 4211 `"Monessen, PA"', add
label define city_lbl 4212 `"Monroe, MI"', add
label define city_lbl 4213 `"Monroe, LA"', add
label define city_lbl 4214 `"Monrovia, CA"', add
label define city_lbl 4230 `"Montclair, NJ"', add
label define city_lbl 4250 `"Montgomery, AL"', add
label define city_lbl 4251 `"Morgantown, WV"', add
label define city_lbl 4252 `"Morristown, NJ"', add
label define city_lbl 4253 `"Moundsville, WV"', add
label define city_lbl 4254 `"Mount Arlington, NJ"', add
label define city_lbl 4255 `"Mount Carmel, PA"', add
label define city_lbl 4256 `"Mount Clemens, MI"', add
label define city_lbl 4270 `"Moreno Valley, CA"', add
label define city_lbl 4290 `"Mount Vernon, NY"', add
label define city_lbl 4291 `"Mount Vernon, IL"', add
label define city_lbl 4310 `"Muncie, IN"', add
label define city_lbl 4311 `"Munhall, PA"', add
label define city_lbl 4312 `"Murphysboro, IL"', add
label define city_lbl 4313 `"Muscatine, IA"', add
label define city_lbl 4330 `"Muskegon, MI"', add
label define city_lbl 4331 `"Muskegon Heights, MI"', add
label define city_lbl 4350 `"Muskogee, OK"', add
label define city_lbl 4351 `"Nanticoke, PA"', add
label define city_lbl 4370 `"Nantucket, MA"', add
label define city_lbl 4390 `"Nashua, NH"', add
label define city_lbl 4410 `"Nashville-Davidson, TN"', add
label define city_lbl 4411 `"Nashville, TN"', add
label define city_lbl 4413 `"Natchez, MS"', add
label define city_lbl 4414 `"Natick, MA"', add
label define city_lbl 4415 `"Naugatuck, CT"', add
label define city_lbl 4416 `"Needham, MA"', add
label define city_lbl 4430 `"New Albany, IN"', add
label define city_lbl 4450 `"New Bedford, MA"', add
label define city_lbl 4451 `"New Bern, NC"', add
label define city_lbl 4452 `"New Brighton, NY"', add
label define city_lbl 4470 `"New Britain, CT"', add
label define city_lbl 4490 `"New Brunswick, NJ"', add
label define city_lbl 4510 `"New Castle, PA"', add
label define city_lbl 4511 `"New Castle, IN"', add
label define city_lbl 4530 `"New Haven, CT"', add
label define city_lbl 4550 `"New London, CT"', add
label define city_lbl 4570 `"New Orleans, LA"', add
label define city_lbl 4571 `"New Philadelphia, OH"', add
label define city_lbl 4590 `"New Rochelle, NY"', add
label define city_lbl 4610 `"New York, NY"', add
label define city_lbl 4611 `"Brooklyn (only in census years before 1900)"', add
label define city_lbl 4630 `"Newark, NJ"', add
label define city_lbl 4650 `"Newark, OH"', add
label define city_lbl 4670 `"Newburgh, NY"', add
label define city_lbl 4690 `"Newburyport, MA"', add
label define city_lbl 4710 `"Newport, KY"', add
label define city_lbl 4730 `"Newport, RI"', add
label define city_lbl 4750 `"Newport News, VA"', add
label define city_lbl 4770 `"Newton, MA"', add
label define city_lbl 4771 `"Newton, IA"', add
label define city_lbl 4772 `"Newton, KS"', add
label define city_lbl 4790 `"Niagara Falls, NY"', add
label define city_lbl 4791 `"Niles, MI"', add
label define city_lbl 4792 `"Niles, OH"', add
label define city_lbl 4810 `"Norfolk, VA"', add
label define city_lbl 4811 `"Norfolk, NE"', add
label define city_lbl 4820 `"North Las Vegas, NV"', add
label define city_lbl 4830 `"Norristown Borough, PA"', add
label define city_lbl 4831 `"North Adams, MA"', add
label define city_lbl 4832 `"North Attleborough, MA"', add
label define city_lbl 4833 `"North Bennington, VT"', add
label define city_lbl 4834 `"North Braddock, PA"', add
label define city_lbl 4835 `"North Branford, CT"', add
label define city_lbl 4836 `"North Haven, CT"', add
label define city_lbl 4837 `"North Little Rock, AR"', add
label define city_lbl 4838 `"North Platte, NE"', add
label define city_lbl 4839 `"North Providence, RI"', add
label define city_lbl 4840 `"Northampton, MA"', add
label define city_lbl 4841 `"North Tonawanda, NY"', add
label define city_lbl 4842 `"North Yakima, WA"', add
label define city_lbl 4843 `"Northbridge, MA"', add
label define city_lbl 4850 `"North Providence, RI"', add
label define city_lbl 4860 `"Norwalk, CA"', add
label define city_lbl 4870 `"Norwalk, CT"', add
label define city_lbl 4890 `"Norwich, CT"', add
label define city_lbl 4900 `"Norwood, OH"', add
label define city_lbl 4901 `"Norwood, MA"', add
label define city_lbl 4902 `"Nutley, NJ"', add
label define city_lbl 4910 `"Oak Park Village"', add
label define city_lbl 4930 `"Oakland, CA"', add
label define city_lbl 4950 `"Oceanside, CA"', add
label define city_lbl 4970 `"Ogden, UT"', add
label define city_lbl 4971 `"Ogdensburg, NY"', add
label define city_lbl 4972 `"Oil City, PA"', add
label define city_lbl 4990 `"Oklahoma City, OK"', add
label define city_lbl 4991 `"Okmulgee, OK"', add
label define city_lbl 4992 `"Old Bennington, VT"', add
label define city_lbl 4993 `"Old Forge, PA"', add
label define city_lbl 4994 `"Olean, NY"', add
label define city_lbl 4995 `"Olympia, WA"', add
label define city_lbl 4996 `"Olyphant, PA"', add
label define city_lbl 5010 `"Omaha, NE"', add
label define city_lbl 5011 `"Oneida, NY"', add
label define city_lbl 5012 `"Oneonta, NY"', add
label define city_lbl 5030 `"Ontario, CA"', add
label define city_lbl 5040 `"Orange, CA"', add
label define city_lbl 5050 `"Orange, NJ"', add
label define city_lbl 5051 `"Orange, CT"', add
label define city_lbl 5070 `"Orlando, FL"', add
label define city_lbl 5090 `"Oshkosh, WI"', add
label define city_lbl 5091 `"Oskaloosa, IA"', add
label define city_lbl 5092 `"Ossining, NY"', add
label define city_lbl 5110 `"Oswego, NY"', add
label define city_lbl 5111 `"Ottawa, IL"', add
label define city_lbl 5112 `"Ottumwa, IA"', add
label define city_lbl 5113 `"Owensboro, KY"', add
label define city_lbl 5114 `"Owosso, MI"', add
label define city_lbl 5116 `"Painesville, OH"', add
label define city_lbl 5117 `"Palestine, TX"', add
label define city_lbl 5118 `"Palo Alto, CA"', add
label define city_lbl 5119 `"Pampa, TX"', add
label define city_lbl 5121 `"Paris, TX"', add
label define city_lbl 5122 `"Park Ridge, IL"', add
label define city_lbl 5123 `"Parkersburg, WV"', add
label define city_lbl 5124 `"Parma, OH"', add
label define city_lbl 5125 `"Parsons, KS"', add
label define city_lbl 5130 `"Oxnard, CA"', add
label define city_lbl 5140 `"Palmdale, CA"', add
label define city_lbl 5150 `"Pasadena, CA"', add
label define city_lbl 5170 `"Pasadena, TX"', add
label define city_lbl 5180 `"Paducah, KY"', add
label define city_lbl 5190 `"Passaic, NJ"', add
label define city_lbl 5210 `"Paterson, NJ"', add
label define city_lbl 5230 `"Pawtucket, RI"', add
label define city_lbl 5231 `"Peabody, MA"', add
label define city_lbl 5232 `"Peekskill, NY"', add
label define city_lbl 5233 `"Pekin, IL"', add
label define city_lbl 5250 `"Pensacola, FL"', add
label define city_lbl 5270 `"Peoria, IL"', add
label define city_lbl 5271 `"Peoria Heights, IL"', add
label define city_lbl 5290 `"Perth Amboy, NJ"', add
label define city_lbl 5291 `"Peru, IN"', add
label define city_lbl 5310 `"Petersburg, VA"', add
label define city_lbl 5311 `"Phenix City, AL"', add
label define city_lbl 5330 `"Philadelphia, PA"', add
label define city_lbl 5331 `"Kensington"', add
label define city_lbl 5332 `"Mayamensing"', add
label define city_lbl 5333 `"Northern Liberties"', add
label define city_lbl 5334 `"Southwark"', add
label define city_lbl 5335 `"Spring Garden"', add
label define city_lbl 5341 `"Phillipsburg, NJ"', add
label define city_lbl 5350 `"Phoenix, AZ"', add
label define city_lbl 5351 `"Phoenixville, PA"', add
label define city_lbl 5352 `"Pine Bluff, AR"', add
label define city_lbl 5353 `"Piqua, OH"', add
label define city_lbl 5354 `"Pittsburg, KS"', add
label define city_lbl 5370 `"Pittsburgh, PA"', add
label define city_lbl 5390 `"Pittsfield, MA"', add
label define city_lbl 5391 `"Pittston, PA"', add
label define city_lbl 5410 `"Plainfield, NJ"', add
label define city_lbl 5411 `"Plattsburg, NY"', add
label define city_lbl 5412 `"Pleasantville, NJ"', add
label define city_lbl 5413 `"Plymouth, PA"', add
label define city_lbl 5414 `"Plymouth, MA"', add
label define city_lbl 5415 `"Pocatello, ID"', add
label define city_lbl 5430 `"Plano, TX"', add
label define city_lbl 5450 `"Pomona, CA"', add
label define city_lbl 5451 `"Ponca City, OK"', add
label define city_lbl 5460 `"Ponce, PR"', add
label define city_lbl 5470 `"Pontiac, MI"', add
label define city_lbl 5471 `"Port Angeles, WA"', add
label define city_lbl 5480 `"Port Arthur, TX"', add
label define city_lbl 5481 `"Port Chester, NY"', add
label define city_lbl 5490 `"Port Huron, MI"', add
label define city_lbl 5491 `"Port Jervis, NY"', add
label define city_lbl 5510 `"Portland, ME"', add
label define city_lbl 5511 `"Portland, IL"', add
label define city_lbl 5530 `"Portland, OR"', add
label define city_lbl 5550 `"Portsmouth, NH"', add
label define city_lbl 5570 `"Portsmouth, OH"', add
label define city_lbl 5590 `"Portsmouth, VA"', add
label define city_lbl 5591 `"Pottstown, PA"', add
label define city_lbl 5610 `"Pottsville, PA"', add
label define city_lbl 5630 `"Poughkeepsie, NY"', add
label define city_lbl 5650 `"Providence, RI"', add
label define city_lbl 5660 `"Provo, UT"', add
label define city_lbl 5670 `"Pueblo, CO"', add
label define city_lbl 5671 `"Punxsutawney, PA"', add
label define city_lbl 5690 `"Quincy, IL"', add
label define city_lbl 5710 `"Quincy, MA"', add
label define city_lbl 5730 `"Racine, WI"', add
label define city_lbl 5731 `"Rahway, NJ"', add
label define city_lbl 5750 `"Raleigh, NC"', add
label define city_lbl 5751 `"Ranger, TX"', add
label define city_lbl 5752 `"Rapid City, SD"', add
label define city_lbl 5770 `"Rancho Cucamonga, CA"', add
label define city_lbl 5790 `"Reading, PA"', add
label define city_lbl 5791 `"Red Bank, NJ"', add
label define city_lbl 5792 `"Redlands, CA"', add
label define city_lbl 5810 `"Reno, NV"', add
label define city_lbl 5811 `"Rensselaer, NY"', add
label define city_lbl 5830 `"Revere, MA"', add
label define city_lbl 5850 `"Richmond, IN"', add
label define city_lbl 5870 `"Richmond, VA"', add
label define city_lbl 5871 `"Richmond, CA"', add
label define city_lbl 5872 `"Ridgefield Park, NJ"', add
label define city_lbl 5873 `"Ridgewood, NJ"', add
label define city_lbl 5874 `"River Rouge, MI"', add
label define city_lbl 5890 `"Riverside, CA"', add
label define city_lbl 5910 `"Roanoke, VA"', add
label define city_lbl 5930 `"Rochester, NY"', add
label define city_lbl 5931 `"Rochester, NH"', add
label define city_lbl 5932 `"Rochester, MN"', add
label define city_lbl 5933 `"Rock Hill, SC"', add
label define city_lbl 5950 `"Rock Island, IL"', add
label define city_lbl 5970 `"Rockford, IL"', add
label define city_lbl 5971 `"Rockland, ME"', add
label define city_lbl 5972 `"Rockton, IL"', add
label define city_lbl 5973 `"Rockville Centre, NY"', add
label define city_lbl 5974 `"Rocky Mount, NC"', add
label define city_lbl 5990 `"Rome, NY"', add
label define city_lbl 5991 `"Rome, GA"', add
label define city_lbl 5992 `"Roosevelt, NJ"', add
label define city_lbl 5993 `"Roselle, NJ"', add
label define city_lbl 5994 `"Roswell, NM"', add
label define city_lbl 6010 `"Roxbury, MA"', add
label define city_lbl 6011 `"Royal Oak, MI"', add
label define city_lbl 6012 `"Rumford Falls, ME"', add
label define city_lbl 6013 `"Rutherford, NJ"', add
label define city_lbl 6014 `"Rutland, VT"', add
label define city_lbl 6030 `"Sacramento, CA"', add
label define city_lbl 6050 `"Saginaw, MI"', add
label define city_lbl 6070 `"Saint Joseph, MO"', add
label define city_lbl 6090 `"Saint Louis, MO"', add
label define city_lbl 6110 `"Saint Paul, MN"', add
label define city_lbl 6130 `"Saint Petersburg, FL"', add
label define city_lbl 6150 `"Salem, MA"', add
label define city_lbl 6170 `"Salem, OR"', add
label define city_lbl 6171 `"Salem, OH"', add
label define city_lbl 6172 `"Salina, KS"', add
label define city_lbl 6190 `"Salinas, CA"', add
label define city_lbl 6191 `"Salisbury, NC"', add
label define city_lbl 6192 `"Salisbury, MD"', add
label define city_lbl 6210 `"Salt Lake City, UT"', add
label define city_lbl 6211 `"San Angelo, TX"', add
label define city_lbl 6220 `"San Angelo, TX"', add
label define city_lbl 6230 `"San Antonio, TX"', add
label define city_lbl 6231 `"San Benito, TX"', add
label define city_lbl 6250 `"San Bernardino, CA"', add
label define city_lbl 6260 `"San Buenaventura (Ventura), CA"', add
label define city_lbl 6270 `"San Diego, CA"', add
label define city_lbl 6280 `"Sandusky, OH"', add
label define city_lbl 6281 `"Sanford, FL"', add
label define city_lbl 6282 `"Sanford, ME"', add
label define city_lbl 6290 `"San Francisco, CA"', add
label define city_lbl 6300 `"San Juan, PR"', add
label define city_lbl 6310 `"San Jose, CA"', add
label define city_lbl 6311 `"San Leandro, CA"', add
label define city_lbl 6312 `"San Mateo, CA"', add
label define city_lbl 6320 `"Santa Barbara, CA"', add
label define city_lbl 6321 `"Santa Cruz, CA"', add
label define city_lbl 6322 `"Santa Fe, NM"', add
label define city_lbl 6330 `"Santa Ana, CA"', add
label define city_lbl 6340 `"Santa Clarita, CA"', add
label define city_lbl 6350 `"Santa Rosa, CA"', add
label define city_lbl 6351 `"Sapulpa, OK"', add
label define city_lbl 6352 `"Saratoga Springs, NY"', add
label define city_lbl 6353 `"Saugus, MA"', add
label define city_lbl 6354 `"Sault Ste. Marie, MI"', add
label define city_lbl 6360 `"Santa Monica, CA"', add
label define city_lbl 6370 `"Savannah, GA"', add
label define city_lbl 6390 `"Schenectedy, NY"', add
label define city_lbl 6410 `"Scranton, PA"', add
label define city_lbl 6430 `"Seattle, WA"', add
label define city_lbl 6431 `"Sedalia, MO"', add
label define city_lbl 6432 `"Selma, AL"', add
label define city_lbl 6433 `"Seminole, OK"', add
label define city_lbl 6434 `"Shaker Heights, OH"', add
label define city_lbl 6435 `"Shamokin, PA"', add
label define city_lbl 6437 `"Sharpsville, PA"', add
label define city_lbl 6438 `"Shawnee, OK"', add
label define city_lbl 6440 `"Sharon, PA"', add
label define city_lbl 6450 `"Sheboygan, WI"', add
label define city_lbl 6451 `"Shelby, NC"', add
label define city_lbl 6452 `"Shelbyville, IN"', add
label define city_lbl 6453 `"Shelton, CT"', add
label define city_lbl 6470 `"Shenandoah Borough, PA"', add
label define city_lbl 6471 `"Sherman, TX"', add
label define city_lbl 6472 `"Shorewood, WI"', add
label define city_lbl 6490 `"Shreveport, LA"', add
label define city_lbl 6500 `"Simi Valley, CA"', add
label define city_lbl 6510 `"Sioux City, IA"', add
label define city_lbl 6530 `"Sioux Falls, SD"', add
label define city_lbl 6550 `"Smithfield, RI (1850)"', add
label define city_lbl 6570 `"Somerville, MA"', add
label define city_lbl 6590 `"South Bend, IN"', add
label define city_lbl 6591 `"South Bethlehem, PA"', add
label define city_lbl 6592 `"South Boise, ID"', add
label define city_lbl 6593 `"South Gate, CA"', add
label define city_lbl 6594 `"South Milwaukee, WI"', add
label define city_lbl 6595 `"South Norwalk, CT"', add
label define city_lbl 6610 `"South Omaha, NE"', add
label define city_lbl 6611 `"South Orange, NJ"', add
label define city_lbl 6612 `"South Pasadena, CA"', add
label define city_lbl 6613 `"South Pittsburgh, PA"', add
label define city_lbl 6614 `"South Portland, ME"', add
label define city_lbl 6615 `"South River, NJ"', add
label define city_lbl 6616 `"South St. Paul, MN"', add
label define city_lbl 6617 `"Southbridge, MA"', add
label define city_lbl 6620 `"Spartanburg, SC"', add
label define city_lbl 6630 `"Spokane, WA"', add
label define city_lbl 6650 `"Springfield, IL"', add
label define city_lbl 6670 `"Springfield, MA"', add
label define city_lbl 6690 `"Springfield, MO"', add
label define city_lbl 6691 `"St. Augustine, FL"', add
label define city_lbl 6692 `"St. Charles, MO"', add
label define city_lbl 6693 `"St. Cloud, MN"', add
label define city_lbl 6710 `"Springfield, OH"', add
label define city_lbl 6730 `"Stamford, CT"', add
label define city_lbl 6731 `"Statesville, NC"', add
label define city_lbl 6732 `"Staunton, VA"', add
label define city_lbl 6733 `"Steelton, PA"', add
label define city_lbl 6734 `"Sterling, IL"', add
label define city_lbl 6750 `"Sterling Heights, MI"', add
label define city_lbl 6770 `"Steubenville, OH"', add
label define city_lbl 6771 `"Stevens Point, WI"', add
label define city_lbl 6772 `"Stillwater, MN"', add
label define city_lbl 6790 `"Stockton, CA"', add
label define city_lbl 6791 `"Stoneham, MA"', add
label define city_lbl 6792 `"Stonington, CT"', add
label define city_lbl 6793 `"Stratford, CT"', add
label define city_lbl 6794 `"Streator, IL"', add
label define city_lbl 6795 `"Struthers, OH"', add
label define city_lbl 6796 `"Suffolk, VA"', add
label define city_lbl 6797 `"Summit, NJ"', add
label define city_lbl 6798 `"Sumter, SC"', add
label define city_lbl 6799 `"Sunbury, PA"', add
label define city_lbl 6810 `"Sunnyvale, CA"', add
label define city_lbl 6830 `"Superior, WI"', add
label define city_lbl 6831 `"Swampscott, MA"', add
label define city_lbl 6832 `"Sweetwater, TX"', add
label define city_lbl 6833 `"Swissvale, PA"', add
label define city_lbl 6850 `"Syracuse, NY"', add
label define city_lbl 6870 `"Tacoma, WA"', add
label define city_lbl 6871 `"Tallahassee, FL"', add
label define city_lbl 6872 `"Tamaqua, PA"', add
label define city_lbl 6890 `"Tampa, FL"', add
label define city_lbl 6910 `"Taunton, MA"', add
label define city_lbl 6911 `"Taylor, PA"', add
label define city_lbl 6912 `"Temple, TX"', add
label define city_lbl 6930 `"Tempe, AZ"', add
label define city_lbl 6950 `"Terre Haute, IN"', add
label define city_lbl 6951 `"Texarkana, TX"', add
label define city_lbl 6952 `"Thomasville, GA"', add
label define city_lbl 6953 `"Thomasville, NC"', add
label define city_lbl 6954 `"Tiffin, OH"', add
label define city_lbl 6960 `"Thousand Oaks, CA"', add
label define city_lbl 6970 `"Toledo, OH"', add
label define city_lbl 6971 `"Tonawanda, NY"', add
label define city_lbl 6990 `"Topeka, KS"', add
label define city_lbl 6991 `"Torrington, CT"', add
label define city_lbl 6992 `"Traverse City, MI"', add
label define city_lbl 7000 `"Torrance, CA"', add
label define city_lbl 7010 `"Trenton, NJ"', add
label define city_lbl 7011 `"Trinidad, CO"', add
label define city_lbl 7030 `"Troy, NY"', add
label define city_lbl 7050 `"Tucson, AZ"', add
label define city_lbl 7070 `"Tulsa, OK"', add
label define city_lbl 7071 `"Turtle Creek, PA"', add
label define city_lbl 7072 `"Tuscaloosa, AL"', add
label define city_lbl 7073 `"Two Rivers, WI"', add
label define city_lbl 7074 `"Tyler, TX"', add
label define city_lbl 7080 `"Union City, NJ"', add
label define city_lbl 7081 `"Uniontown, PA"', add
label define city_lbl 7082 `"University City, MO"', add
label define city_lbl 7083 `"Urbana, IL"', add
label define city_lbl 7090 `"Utica, NY"', add
label define city_lbl 7091 `"Valdosta, GA"', add
label define city_lbl 7092 `"Vallejo, CA"', add
label define city_lbl 7093 `"Valley Stream, NY"', add
label define city_lbl 7100 `"Vancouver, WA"', add
label define city_lbl 7110 `"Vallejo, CA"', add
label define city_lbl 7111 `"Vandergrift, PA"', add
label define city_lbl 7112 `"Venice, CA"', add
label define city_lbl 7120 `"Vicksburg, MS"', add
label define city_lbl 7121 `"Vincennes, IN"', add
label define city_lbl 7122 `"Virginia, MN"', add
label define city_lbl 7123 `"Virginia City, NV"', add
label define city_lbl 7130 `"Virginia Beach, VA"', add
label define city_lbl 7140 `"University City, MO"', add
label define city_lbl 7150 `"Waco, TX"', add
label define city_lbl 7151 `"Wakefield, MA"', add
label define city_lbl 7152 `"Walla Walla, WA"', add
label define city_lbl 7153 `"Wallingford, CT"', add
label define city_lbl 7170 `"Waltham, MA"', add
label define city_lbl 7180 `"Warren, MI"', add
label define city_lbl 7190 `"Warren, OH"', add
label define city_lbl 7191 `"Warren, PA"', add
label define city_lbl 7210 `"Warwick Town, RI"', add
label define city_lbl 7230 `"Washington, DC"', add
label define city_lbl 7231 `"Georgetown, DC"', add
label define city_lbl 7241 `"Washington, PA"', add
label define city_lbl 7242 `"Washington, VA"', add
label define city_lbl 7250 `"Waterbury, CT"', add
label define city_lbl 7270 `"Waterloo, IA"', add
label define city_lbl 7290 `"Waterloo, NY"', add
label define city_lbl 7310 `"Watertown, NY"', add
label define city_lbl 7311 `"Watertown, WI"', add
label define city_lbl 7312 `"Watertown, SD"', add
label define city_lbl 7313 `"Watertown, MA"', add
label define city_lbl 7314 `"Waterville, ME"', add
label define city_lbl 7315 `"Watervliet, NY"', add
label define city_lbl 7316 `"Waukegan, IL"', add
label define city_lbl 7317 `"Waukesha, WI"', add
label define city_lbl 7318 `"Wausau, WI"', add
label define city_lbl 7319 `"Wauwatosa, WI"', add
label define city_lbl 7320 `"West Covina, CA"', add
label define city_lbl 7321 `"Waycross, GA"', add
label define city_lbl 7322 `"Waynesboro, PA"', add
label define city_lbl 7323 `"Webb City, MO"', add
label define city_lbl 7324 `"Webster Groves, MO"', add
label define city_lbl 7325 `"Webster, MA"', add
label define city_lbl 7326 `"Wellesley, MA"', add
label define city_lbl 7327 `"Wenatchee, WA"', add
label define city_lbl 7329 `"West Bay City, MI"', add
label define city_lbl 7330 `"West Hoboken, NJ"', add
label define city_lbl 7331 `"West Bethlehem, PA"', add
label define city_lbl 7332 `"West Chester, PA"', add
label define city_lbl 7333 `"West Frankfort, IL"', add
label define city_lbl 7334 `"West Hartford, CT"', add
label define city_lbl 7335 `"West Haven, CT"', add
label define city_lbl 7340 `"West Allis, WI"', add
label define city_lbl 7350 `"West New York, NJ"', add
label define city_lbl 7351 `"West Orange, NJ"', add
label define city_lbl 7352 `"West Palm Beach, FL"', add
label define city_lbl 7353 `"West Springfield, MA"', add
label define city_lbl 7370 `"West Troy, NY"', add
label define city_lbl 7371 `"West Warwick, RI"', add
label define city_lbl 7372 `"Westbrook, ME"', add
label define city_lbl 7373 `"Westerly, RI"', add
label define city_lbl 7374 `"Westfield, MA"', add
label define city_lbl 7375 `"Westfield, NJ"', add
label define city_lbl 7376 `"Wewoka, OK"', add
label define city_lbl 7377 `"Weymouth, MA"', add
label define city_lbl 7390 `"Wheeling, WV"', add
label define city_lbl 7400 `"White Plains, NY"', add
label define city_lbl 7401 `"Whiting, IN"', add
label define city_lbl 7402 `"Whittier, CA"', add
label define city_lbl 7410 `"Wichita, KS"', add
label define city_lbl 7430 `"Wichita Falls, TX"', add
label define city_lbl 7450 `"Wilkes-Barre, PA"', add
label define city_lbl 7451 `"Wilkinsburg, PA"', add
label define city_lbl 7460 `"Wilkinsburg, PA"', add
label define city_lbl 7470 `"Williamsport, PA"', add
label define city_lbl 7471 `"Willimantic, CT"', add
label define city_lbl 7472 `"Wilmette, IL"', add
label define city_lbl 7490 `"Wilmington, DE"', add
label define city_lbl 7510 `"Wilmington, NC"', add
label define city_lbl 7511 `"Wilson, NC"', add
label define city_lbl 7512 `"Winchester, VA"', add
label define city_lbl 7513 `"Winchester, MA"', add
label define city_lbl 7514 `"Windham, CT"', add
label define city_lbl 7515 `"Winnetka, IL"', add
label define city_lbl 7516 `"Winona, MN"', add
label define city_lbl 7530 `"Winston-Salem, NC"', add
label define city_lbl 7531 `"Winthrop, MA"', add
label define city_lbl 7532 `"Woburn, MA"', add
label define city_lbl 7533 `"Woodlawn, PA"', add
label define city_lbl 7534 `"Woodmont, CT"', add
label define city_lbl 7550 `"Woonsocket, RI"', add
label define city_lbl 7551 `"Wooster, OH"', add
label define city_lbl 7570 `"Worcester, MA"', add
label define city_lbl 7571 `"Wyandotte, MI"', add
label define city_lbl 7572 `"Xenia, OH"', add
label define city_lbl 7573 `"Yakima, WA"', add
label define city_lbl 7590 `"Yonkers, NY"', add
label define city_lbl 7610 `"York, PA"', add
label define city_lbl 7630 `"Youngstown, OH"', add
label define city_lbl 7631 `"Ypsilanti, MI"', add
label define city_lbl 7650 `"Zanesville, OH"', add
label values city city_lbl

label define citypop_lbl 00000 `"0"'
label define citypop_lbl 00347 `"00347"', add
label define citypop_lbl 00354 `"00354"', add
label define citypop_lbl 00367 `"00367"', add
label define citypop_lbl 00372 `"00372"', add
label define citypop_lbl 00432 `"00432"', add
label define citypop_lbl 00434 `"00434"', add
label define citypop_lbl 00436 `"00436"', add
label define citypop_lbl 00469 `"00469"', add
label define citypop_lbl 00490 `"00490"', add
label define citypop_lbl 00497 `"00497"', add
label define citypop_lbl 00524 `"00524"', add
label define citypop_lbl 00534 `"00534"', add
label define citypop_lbl 00547 `"00547"', add
label define citypop_lbl 00585 `"00585"', add
label define citypop_lbl 00607 `"00607"', add
label define citypop_lbl 00624 `"00624"', add
label define citypop_lbl 00631 `"00631"', add
label define citypop_lbl 00654 `"00654"', add
label define citypop_lbl 00657 `"00657"', add
label define citypop_lbl 00674 `"00674"', add
label define citypop_lbl 00679 `"00679"', add
label define citypop_lbl 00685 `"00685"', add
label define citypop_lbl 00709 `"00709"', add
label define citypop_lbl 00735 `"00735"', add
label define citypop_lbl 00751 `"00751"', add
label define citypop_lbl 00787 `"00787"', add
label define citypop_lbl 00788 `"00788"', add
label define citypop_lbl 00798 `"00798"', add
label define citypop_lbl 00802 `"00802"', add
label define citypop_lbl 00805 `"00805"', add
label define citypop_lbl 00824 `"00824"', add
label define citypop_lbl 00835 `"00835"', add
label define citypop_lbl 00842 `"00842"', add
label define citypop_lbl 00843 `"00843"', add
label define citypop_lbl 00846 `"00846"', add
label define citypop_lbl 00852 `"00852"', add
label define citypop_lbl 00869 `"00869"', add
label define citypop_lbl 00878 `"00878"', add
label define citypop_lbl 00879 `"00879"', add
label define citypop_lbl 00881 `"00881"', add
label define citypop_lbl 00913 `"00913"', add
label define citypop_lbl 00917 `"00917"', add
label define citypop_lbl 00921 `"00921"', add
label define citypop_lbl 00929 `"00929"', add
label define citypop_lbl 00953 `"00953"', add
label define citypop_lbl 00961 `"00961"', add
label define citypop_lbl 00967 `"00967"', add
label define citypop_lbl 00968 `"00968"', add
label define citypop_lbl 00972 `"00972"', add
label define citypop_lbl 00976 `"00976"', add
label define citypop_lbl 00977 `"00977"', add
label define citypop_lbl 00984 `"00984"', add
label define citypop_lbl 00993 `"00993"', add
label define citypop_lbl 00996 `"00996"', add
label define citypop_lbl 01002 `"01002"', add
label define citypop_lbl 01003 `"01003"', add
label define citypop_lbl 01004 `"01004"', add
label define citypop_lbl 01005 `"01005"', add
label define citypop_lbl 01007 `"01007"', add
label define citypop_lbl 01008 `"01008"', add
label define citypop_lbl 01009 `"01009"', add
label define citypop_lbl 01011 `"01011"', add
label define citypop_lbl 01012 `"01012"', add
label define citypop_lbl 01013 `"01013"', add
label define citypop_lbl 01014 `"01014"', add
label define citypop_lbl 01015 `"01015"', add
label define citypop_lbl 01020 `"01020"', add
label define citypop_lbl 01021 `"01021"', add
label define citypop_lbl 01023 `"01023"', add
label define citypop_lbl 01025 `"01025"', add
label define citypop_lbl 01027 `"01027"', add
label define citypop_lbl 01028 `"01028"', add
label define citypop_lbl 01032 `"01032"', add
label define citypop_lbl 01033 `"01033"', add
label define citypop_lbl 01034 `"01034"', add
label define citypop_lbl 01037 `"01037"', add
label define citypop_lbl 01038 `"01038"', add
label define citypop_lbl 01040 `"01040"', add
label define citypop_lbl 01041 `"01041"', add
label define citypop_lbl 01042 `"01042"', add
label define citypop_lbl 01043 `"01043"', add
label define citypop_lbl 01045 `"01045"', add
label define citypop_lbl 01046 `"01046"', add
label define citypop_lbl 01048 `"01048"', add
label define citypop_lbl 01049 `"01049"', add
label define citypop_lbl 01051 `"01051"', add
label define citypop_lbl 01052 `"01052"', add
label define citypop_lbl 01055 `"01055"', add
label define citypop_lbl 01056 `"01056"', add
label define citypop_lbl 01059 `"01059"', add
label define citypop_lbl 01060 `"01060"', add
label define citypop_lbl 01062 `"01062"', add
label define citypop_lbl 01064 `"01064"', add
label define citypop_lbl 01066 `"01066"', add
label define citypop_lbl 01067 `"01067"', add
label define citypop_lbl 01068 `"01068"', add
label define citypop_lbl 01069 `"01069"', add
label define citypop_lbl 01070 `"01070"', add
label define citypop_lbl 01073 `"01073"', add
label define citypop_lbl 01075 `"01075"', add
label define citypop_lbl 01077 `"01077"', add
label define citypop_lbl 01078 `"01078"', add
label define citypop_lbl 01080 `"01080"', add
label define citypop_lbl 01081 `"01081"', add
label define citypop_lbl 01082 `"01082"', add
label define citypop_lbl 01083 `"01083"', add
label define citypop_lbl 01084 `"01084"', add
label define citypop_lbl 01086 `"01086"', add
label define citypop_lbl 01087 `"01087"', add
label define citypop_lbl 01088 `"01088"', add
label define citypop_lbl 01090 `"01090"', add
label define citypop_lbl 01092 `"01092"', add
label define citypop_lbl 01094 `"01094"', add
label define citypop_lbl 01095 `"01095"', add
label define citypop_lbl 01096 `"01096"', add
label define citypop_lbl 01097 `"01097"', add
label define citypop_lbl 01098 `"01098"', add
label define citypop_lbl 01099 `"01099"', add
label define citypop_lbl 01100 `"01100"', add
label define citypop_lbl 01103 `"01103"', add
label define citypop_lbl 01112 `"01112"', add
label define citypop_lbl 01114 `"01114"', add
label define citypop_lbl 01115 `"01115"', add
label define citypop_lbl 01116 `"01116"', add
label define citypop_lbl 01118 `"01118"', add
label define citypop_lbl 01119 `"01119"', add
label define citypop_lbl 01125 `"01125"', add
label define citypop_lbl 01126 `"01126"', add
label define citypop_lbl 01127 `"01127"', add
label define citypop_lbl 01129 `"01129"', add
label define citypop_lbl 01131 `"01131"', add
label define citypop_lbl 01132 `"01132"', add
label define citypop_lbl 01133 `"01133"', add
label define citypop_lbl 01134 `"01134"', add
label define citypop_lbl 01135 `"01135"', add
label define citypop_lbl 01139 `"01139"', add
label define citypop_lbl 01140 `"01140"', add
label define citypop_lbl 01141 `"01141"', add
label define citypop_lbl 01142 `"01142"', add
label define citypop_lbl 01143 `"01143"', add
label define citypop_lbl 01145 `"01145"', add
label define citypop_lbl 01148 `"01148"', add
label define citypop_lbl 01149 `"01149"', add
label define citypop_lbl 01151 `"01151"', add
label define citypop_lbl 01153 `"01153"', add
label define citypop_lbl 01154 `"01154"', add
label define citypop_lbl 01155 `"01155"', add
label define citypop_lbl 01156 `"01156"', add
label define citypop_lbl 01157 `"01157"', add
label define citypop_lbl 01159 `"01159"', add
label define citypop_lbl 01160 `"01160"', add
label define citypop_lbl 01163 `"01163"', add
label define citypop_lbl 01165 `"01165"', add
label define citypop_lbl 01166 `"01166"', add
label define citypop_lbl 01167 `"01167"', add
label define citypop_lbl 01169 `"01169"', add
label define citypop_lbl 01170 `"01170"', add
label define citypop_lbl 01171 `"01171"', add
label define citypop_lbl 01172 `"01172"', add
label define citypop_lbl 01174 `"01174"', add
label define citypop_lbl 01175 `"01175"', add
label define citypop_lbl 01178 `"01178"', add
label define citypop_lbl 01181 `"01181"', add
label define citypop_lbl 01182 `"01182"', add
label define citypop_lbl 01184 `"01184"', add
label define citypop_lbl 01187 `"01187"', add
label define citypop_lbl 01188 `"01188"', add
label define citypop_lbl 01191 `"01191"', add
label define citypop_lbl 01193 `"01193"', add
label define citypop_lbl 01194 `"01194"', add
label define citypop_lbl 01196 `"01196"', add
label define citypop_lbl 01198 `"01198"', add
label define citypop_lbl 01200 `"01200"', add
label define citypop_lbl 01206 `"01206"', add
label define citypop_lbl 01210 `"01210"', add
label define citypop_lbl 01213 `"01213"', add
label define citypop_lbl 01216 `"01216"', add
label define citypop_lbl 01217 `"01217"', add
label define citypop_lbl 01220 `"01220"', add
label define citypop_lbl 01226 `"01226"', add
label define citypop_lbl 01232 `"01232"', add
label define citypop_lbl 01233 `"01233"', add
label define citypop_lbl 01236 `"01236"', add
label define citypop_lbl 01237 `"01237"', add
label define citypop_lbl 01240 `"01240"', add
label define citypop_lbl 01242 `"01242"', add
label define citypop_lbl 01243 `"01243"', add
label define citypop_lbl 01244 `"01244"', add
label define citypop_lbl 01245 `"01245"', add
label define citypop_lbl 01247 `"01247"', add
label define citypop_lbl 01248 `"01248"', add
label define citypop_lbl 01249 `"01249"', add
label define citypop_lbl 01250 `"01250"', add
label define citypop_lbl 01256 `"01256"', add
label define citypop_lbl 01260 `"01260"', add
label define citypop_lbl 01261 `"01261"', add
label define citypop_lbl 01262 `"01262"', add
label define citypop_lbl 01264 `"01264"', add
label define citypop_lbl 01265 `"01265"', add
label define citypop_lbl 01273 `"01273"', add
label define citypop_lbl 01277 `"01277"', add
label define citypop_lbl 01279 `"01279"', add
label define citypop_lbl 01280 `"01280"', add
label define citypop_lbl 01282 `"01282"', add
label define citypop_lbl 01283 `"01283"', add
label define citypop_lbl 01284 `"01284"', add
label define citypop_lbl 01287 `"01287"', add
label define citypop_lbl 01288 `"01288"', add
label define citypop_lbl 01289 `"01289"', add
label define citypop_lbl 01290 `"01290"', add
label define citypop_lbl 01295 `"01295"', add
label define citypop_lbl 01300 `"01300"', add
label define citypop_lbl 01304 `"01304"', add
label define citypop_lbl 01305 `"01305"', add
label define citypop_lbl 01308 `"01308"', add
label define citypop_lbl 01310 `"01310"', add
label define citypop_lbl 01315 `"01315"', add
label define citypop_lbl 01316 `"01316"', add
label define citypop_lbl 01317 `"01317"', add
label define citypop_lbl 01319 `"01319"', add
label define citypop_lbl 01325 `"01325"', add
label define citypop_lbl 01329 `"01329"', add
label define citypop_lbl 01330 `"01330"', add
label define citypop_lbl 01332 `"01332"', add
label define citypop_lbl 01336 `"01336"', add
label define citypop_lbl 01338 `"01338"', add
label define citypop_lbl 01339 `"01339"', add
label define citypop_lbl 01340 `"01340"', add
label define citypop_lbl 01341 `"01341"', add
label define citypop_lbl 01346 `"01346"', add
label define citypop_lbl 01351 `"01351"', add
label define citypop_lbl 01352 `"01352"', add
label define citypop_lbl 01353 `"01353"', add
label define citypop_lbl 01364 `"01364"', add
label define citypop_lbl 01370 `"01370"', add
label define citypop_lbl 01376 `"01376"', add
label define citypop_lbl 01379 `"01379"', add
label define citypop_lbl 01380 `"01380"', add
label define citypop_lbl 01382 `"01382"', add
label define citypop_lbl 01388 `"01388"', add
label define citypop_lbl 01389 `"01389"', add
label define citypop_lbl 01394 `"01394"', add
label define citypop_lbl 01395 `"01395"', add
label define citypop_lbl 01397 `"01397"', add
label define citypop_lbl 01400 `"01400"', add
label define citypop_lbl 01405 `"01405"', add
label define citypop_lbl 01407 `"01407"', add
label define citypop_lbl 01408 `"01408"', add
label define citypop_lbl 01409 `"01409"', add
label define citypop_lbl 01417 `"01417"', add
label define citypop_lbl 01422 `"01422"', add
label define citypop_lbl 01423 `"01423"', add
label define citypop_lbl 01424 `"01424"', add
label define citypop_lbl 01425 `"01425"', add
label define citypop_lbl 01426 `"01426"', add
label define citypop_lbl 01431 `"01431"', add
label define citypop_lbl 01435 `"01435"', add
label define citypop_lbl 01436 `"01436"', add
label define citypop_lbl 01437 `"01437"', add
label define citypop_lbl 01438 `"01438"', add
label define citypop_lbl 01441 `"01441"', add
label define citypop_lbl 01442 `"01442"', add
label define citypop_lbl 01448 `"01448"', add
label define citypop_lbl 01449 `"01449"', add
label define citypop_lbl 01450 `"01450"', add
label define citypop_lbl 01452 `"01452"', add
label define citypop_lbl 01453 `"01453"', add
label define citypop_lbl 01464 `"01464"', add
label define citypop_lbl 01465 `"01465"', add
label define citypop_lbl 01471 `"01471"', add
label define citypop_lbl 01473 `"01473"', add
label define citypop_lbl 01487 `"01487"', add
label define citypop_lbl 01492 `"01492"', add
label define citypop_lbl 01494 `"01494"', add
label define citypop_lbl 01495 `"01495"', add
label define citypop_lbl 01498 `"01498"', add
label define citypop_lbl 01500 `"01500"', add
label define citypop_lbl 01501 `"01501"', add
label define citypop_lbl 01503 `"01503"', add
label define citypop_lbl 01508 `"01508"', add
label define citypop_lbl 01511 `"01511"', add
label define citypop_lbl 01512 `"01512"', add
label define citypop_lbl 01514 `"01514"', add
label define citypop_lbl 01515 `"01515"', add
label define citypop_lbl 01516 `"01516"', add
label define citypop_lbl 01520 `"01520"', add
label define citypop_lbl 01521 `"01521"', add
label define citypop_lbl 01523 `"01523"', add
label define citypop_lbl 01525 `"01525"', add
label define citypop_lbl 01526 `"01526"', add
label define citypop_lbl 01533 `"01533"', add
label define citypop_lbl 01536 `"01536"', add
label define citypop_lbl 01539 `"01539"', add
label define citypop_lbl 01543 `"01543"', add
label define citypop_lbl 01548 `"01548"', add
label define citypop_lbl 01550 `"01550"', add
label define citypop_lbl 01551 `"01551"', add
label define citypop_lbl 01552 `"01552"', add
label define citypop_lbl 01556 `"01556"', add
label define citypop_lbl 01557 `"01557"', add
label define citypop_lbl 01568 `"01568"', add
label define citypop_lbl 01570 `"01570"', add
label define citypop_lbl 01572 `"01572"', add
label define citypop_lbl 01576 `"01576"', add
label define citypop_lbl 01580 `"01580"', add
label define citypop_lbl 01582 `"01582"', add
label define citypop_lbl 01585 `"01585"', add
label define citypop_lbl 01586 `"01586"', add
label define citypop_lbl 01587 `"01587"', add
label define citypop_lbl 01589 `"01589"', add
label define citypop_lbl 01596 `"01596"', add
label define citypop_lbl 01598 `"01598"', add
label define citypop_lbl 01601 `"01601"', add
label define citypop_lbl 01606 `"01606"', add
label define citypop_lbl 01607 `"01607"', add
label define citypop_lbl 01617 `"01617"', add
label define citypop_lbl 01618 `"01618"', add
label define citypop_lbl 01630 `"01630"', add
label define citypop_lbl 01639 `"01639"', add
label define citypop_lbl 01642 `"01642"', add
label define citypop_lbl 01643 `"01643"', add
label define citypop_lbl 01644 `"01644"', add
label define citypop_lbl 01645 `"01645"', add
label define citypop_lbl 01647 `"01647"', add
label define citypop_lbl 01651 `"01651"', add
label define citypop_lbl 01652 `"01652"', add
label define citypop_lbl 01662 `"01662"', add
label define citypop_lbl 01663 `"01663"', add
label define citypop_lbl 01674 `"01674"', add
label define citypop_lbl 01677 `"01677"', add
label define citypop_lbl 01680 `"01680"', add
label define citypop_lbl 01681 `"01681"', add
label define citypop_lbl 01683 `"01683"', add
label define citypop_lbl 01684 `"01684"', add
label define citypop_lbl 01693 `"01693"', add
label define citypop_lbl 01694 `"01694"', add
label define citypop_lbl 01696 `"01696"', add
label define citypop_lbl 01698 `"01698"', add
label define citypop_lbl 01700 `"01700"', add
label define citypop_lbl 01701 `"01701"', add
label define citypop_lbl 01704 `"01704"', add
label define citypop_lbl 01705 `"01705"', add
label define citypop_lbl 01707 `"01707"', add
label define citypop_lbl 01706 `"01706"', add
label define citypop_lbl 01709 `"01709"', add
label define citypop_lbl 01713 `"01713"', add
label define citypop_lbl 01722 `"01722"', add
label define citypop_lbl 01724 `"01724"', add
label define citypop_lbl 01726 `"01726"', add
label define citypop_lbl 01727 `"01727"', add
label define citypop_lbl 01731 `"01731"', add
label define citypop_lbl 01734 `"01734"', add
label define citypop_lbl 01736 `"01736"', add
label define citypop_lbl 01739 `"01739"', add
label define citypop_lbl 01743 `"01743"', add
label define citypop_lbl 01744 `"01744"', add
label define citypop_lbl 01748 `"01748"', add
label define citypop_lbl 01750 `"01750"', add
label define citypop_lbl 01753 `"01753"', add
label define citypop_lbl 01755 `"01755"', add
label define citypop_lbl 01756 `"01756"', add
label define citypop_lbl 01765 `"01765"', add
label define citypop_lbl 01766 `"01766"', add
label define citypop_lbl 01770 `"01770"', add
label define citypop_lbl 01771 `"01771"', add
label define citypop_lbl 01772 `"01772"', add
label define citypop_lbl 01774 `"01774"', add
label define citypop_lbl 01776 `"01776"', add
label define citypop_lbl 01777 `"01777"', add
label define citypop_lbl 01780 `"01780"', add
label define citypop_lbl 01781 `"01781"', add
label define citypop_lbl 01783 `"01783"', add
label define citypop_lbl 01789 `"01789"', add
label define citypop_lbl 01800 `"01800"', add
label define citypop_lbl 01802 `"01802"', add
label define citypop_lbl 01805 `"01805"', add
label define citypop_lbl 01806 `"01806"', add
label define citypop_lbl 01807 `"01807"', add
label define citypop_lbl 01815 `"01815"', add
label define citypop_lbl 01817 `"01817"', add
label define citypop_lbl 01818 `"01818"', add
label define citypop_lbl 01821 `"01821"', add
label define citypop_lbl 01823 `"01823"', add
label define citypop_lbl 01827 `"01827"', add
label define citypop_lbl 01831 `"01831"', add
label define citypop_lbl 01835 `"01835"', add
label define citypop_lbl 01836 `"01836"', add
label define citypop_lbl 01839 `"01839"', add
label define citypop_lbl 01843 `"01843"', add
label define citypop_lbl 01844 `"01844"', add
label define citypop_lbl 01845 `"01845"', add
label define citypop_lbl 01854 `"01854"', add
label define citypop_lbl 01858 `"01858"', add
label define citypop_lbl 01863 `"01863"', add
label define citypop_lbl 01871 `"01871"', add
label define citypop_lbl 01875 `"01875"', add
label define citypop_lbl 01880 `"01880"', add
label define citypop_lbl 01881 `"01881"', add
label define citypop_lbl 01887 `"01887"', add
label define citypop_lbl 01889 `"01889"', add
label define citypop_lbl 01891 `"01891"', add
label define citypop_lbl 01895 `"01895"', add
label define citypop_lbl 01896 `"01896"', add
label define citypop_lbl 01908 `"01908"', add
label define citypop_lbl 01909 `"01909"', add
label define citypop_lbl 01910 `"01910"', add
label define citypop_lbl 01916 `"01916"', add
label define citypop_lbl 01928 `"01928"', add
label define citypop_lbl 01929 `"01929"', add
label define citypop_lbl 01930 `"01930"', add
label define citypop_lbl 01931 `"01931"', add
label define citypop_lbl 01932 `"01932"', add
label define citypop_lbl 01936 `"01936"', add
label define citypop_lbl 01937 `"01937"', add
label define citypop_lbl 01938 `"01938"', add
label define citypop_lbl 01939 `"01939"', add
label define citypop_lbl 01940 `"01940"', add
label define citypop_lbl 01944 `"01944"', add
label define citypop_lbl 01947 `"01947"', add
label define citypop_lbl 01950 `"01950"', add
label define citypop_lbl 01954 `"01954"', add
label define citypop_lbl 01956 `"01956"', add
label define citypop_lbl 01957 `"01957"', add
label define citypop_lbl 01961 `"01961"', add
label define citypop_lbl 01963 `"01963"', add
label define citypop_lbl 01965 `"01965"', add
label define citypop_lbl 01966 `"01966"', add
label define citypop_lbl 01970 `"01970"', add
label define citypop_lbl 01976 `"01976"', add
label define citypop_lbl 01978 `"01978"', add
label define citypop_lbl 01979 `"01979"', add
label define citypop_lbl 01981 `"01981"', add
label define citypop_lbl 01982 `"01982"', add
label define citypop_lbl 01985 `"01985"', add
label define citypop_lbl 01986 `"01986"', add
label define citypop_lbl 01987 `"01987"', add
label define citypop_lbl 01989 `"01989"', add
label define citypop_lbl 01900 `"01900"', add
label define citypop_lbl 01992 `"01992"', add
label define citypop_lbl 01993 `"01993"', add
label define citypop_lbl 01994 `"01994"', add
label define citypop_lbl 01995 `"01995"', add
label define citypop_lbl 01998 `"01998"', add
label define citypop_lbl 02002 `"02002"', add
label define citypop_lbl 02004 `"02004"', add
label define citypop_lbl 02016 `"02016"', add
label define citypop_lbl 02020 `"02020"', add
label define citypop_lbl 02029 `"02029"', add
label define citypop_lbl 02031 `"02031"', add
label define citypop_lbl 02033 `"02033"', add
label define citypop_lbl 02034 `"02034"', add
label define citypop_lbl 02035 `"02035"', add
label define citypop_lbl 02037 `"02037"', add
label define citypop_lbl 02042 `"02042"', add
label define citypop_lbl 02044 `"02044"', add
label define citypop_lbl 02045 `"02045"', add
label define citypop_lbl 02047 `"02047"', add
label define citypop_lbl 02057 `"02057"', add
label define citypop_lbl 02058 `"02058"', add
label define citypop_lbl 02060 `"02060"', add
label define citypop_lbl 02069 `"02069"', add
label define citypop_lbl 02081 `"02081"', add
label define citypop_lbl 02097 `"02097"', add
label define citypop_lbl 02103 `"02103"', add
label define citypop_lbl 02106 `"02106"', add
label define citypop_lbl 02107 `"02107"', add
label define citypop_lbl 02108 `"02108"', add
label define citypop_lbl 02109 `"02109"', add
label define citypop_lbl 02152 `"02152"', add
label define citypop_lbl 02171 `"02171"', add
label define citypop_lbl 02182 `"02182"', add
label define citypop_lbl 02192 `"02192"', add
label define citypop_lbl 02193 `"02193"', add
label define citypop_lbl 02194 `"02194"', add
label define citypop_lbl 02195 `"02195"', add
label define citypop_lbl 02198 `"02198"', add
label define citypop_lbl 02206 `"02206"', add
label define citypop_lbl 02211 `"02211"', add
label define citypop_lbl 02215 `"02215"', add
label define citypop_lbl 02220 `"02220"', add
label define citypop_lbl 02221 `"02221"', add
label define citypop_lbl 02230 `"02230"', add
label define citypop_lbl 02234 `"02234"', add
label define citypop_lbl 02235 `"02235"', add
label define citypop_lbl 02238 `"02238"', add
label define citypop_lbl 02239 `"02239"', add
label define citypop_lbl 02254 `"02254"', add
label define citypop_lbl 02263 `"02263"', add
label define citypop_lbl 02265 `"02265"', add
label define citypop_lbl 02278 `"02278"', add
label define citypop_lbl 02285 `"02285"', add
label define citypop_lbl 02291 `"02291"', add
label define citypop_lbl 02296 `"02296"', add
label define citypop_lbl 02303 `"02303"', add
label define citypop_lbl 02308 `"02308"', add
label define citypop_lbl 02316 `"02316"', add
label define citypop_lbl 02335 `"02335"', add
label define citypop_lbl 02344 `"02344"', add
label define citypop_lbl 02369 `"02369"', add
label define citypop_lbl 02372 `"02372"', add
label define citypop_lbl 02386 `"02386"', add
label define citypop_lbl 02394 `"02394"', add
label define citypop_lbl 02401 `"02401"', add
label define citypop_lbl 02417 `"02417"', add
label define citypop_lbl 02418 `"02418"', add
label define citypop_lbl 02435 `"02435"', add
label define citypop_lbl 02438 `"02438"', add
label define citypop_lbl 02439 `"02439"', add
label define citypop_lbl 02448 `"02448"', add
label define citypop_lbl 02470 `"02470"', add
label define citypop_lbl 02471 `"02471"', add
label define citypop_lbl 02486 `"02486"', add
label define citypop_lbl 02487 `"02487"', add
label define citypop_lbl 02493 `"02493"', add
label define citypop_lbl 02507 `"02507"', add
label define citypop_lbl 02511 `"02511"', add
label define citypop_lbl 02535 `"02535"', add
label define citypop_lbl 02539 `"02539"', add
label define citypop_lbl 02550 `"02550"', add
label define citypop_lbl 02552 `"02552"', add
label define citypop_lbl 02583 `"02583"', add
label define citypop_lbl 02603 `"02603"', add
label define citypop_lbl 02605 `"02605"', add
label define citypop_lbl 02612 `"02612"', add
label define citypop_lbl 02617 `"02617"', add
label define citypop_lbl 02622 `"02622"', add
label define citypop_lbl 02664 `"02664"', add
label define citypop_lbl 02670 `"02670"', add
label define citypop_lbl 02676 `"02676"', add
label define citypop_lbl 02679 `"02679"', add
label define citypop_lbl 02691 `"02691"', add
label define citypop_lbl 02702 `"02702"', add
label define citypop_lbl 02708 `"02708"', add
label define citypop_lbl 02715 `"02715"', add
label define citypop_lbl 02722 `"02722"', add
label define citypop_lbl 02735 `"02735"', add
label define citypop_lbl 02751 `"02751"', add
label define citypop_lbl 02757 `"02757"', add
label define citypop_lbl 02761 `"02761"', add
label define citypop_lbl 02775 `"02775"', add
label define citypop_lbl 02787 `"02787"', add
label define citypop_lbl 02788 `"02788"', add
label define citypop_lbl 02791 `"02791"', add
label define citypop_lbl 02793 `"02793"', add
label define citypop_lbl 02800 `"02800"', add
label define citypop_lbl 02811 `"02811"', add
label define citypop_lbl 02814 `"02814"', add
label define citypop_lbl 02823 `"02823"', add
label define citypop_lbl 02830 `"02830"', add
label define citypop_lbl 02844 `"02844"', add
label define citypop_lbl 02853 `"02853"', add
label define citypop_lbl 02863 `"02863"', add
label define citypop_lbl 02871 `"02871"', add
label define citypop_lbl 02872 `"02872"', add
label define citypop_lbl 02901 `"02901"', add
label define citypop_lbl 02904 `"02904"', add
label define citypop_lbl 02924 `"02924"', add
label define citypop_lbl 02926 `"02926"', add
label define citypop_lbl 02937 `"02937"', add
label define citypop_lbl 02938 `"02938"', add
label define citypop_lbl 02947 `"02947"', add
label define citypop_lbl 02958 `"02958"', add
label define citypop_lbl 02969 `"02969"', add
label define citypop_lbl 02984 `"02984"', add
label define citypop_lbl 02985 `"02985"', add
label define citypop_lbl 03017 `"03017"', add
label define citypop_lbl 03023 `"03023"', add
label define citypop_lbl 03036 `"03036"', add
label define citypop_lbl 03040 `"03040"', add
label define citypop_lbl 03054 `"03054"', add
label define citypop_lbl 03061 `"03061"', add
label define citypop_lbl 03084 `"03084"', add
label define citypop_lbl 03128 `"03128"', add
label define citypop_lbl 03136 `"03136"', add
label define citypop_lbl 03143 `"03143"', add
label define citypop_lbl 03144 `"03144"', add
label define citypop_lbl 03163 `"03163"', add
label define citypop_lbl 03191 `"03191"', add
label define citypop_lbl 03224 `"03224"', add
label define citypop_lbl 03250 `"03250"', add
label define citypop_lbl 03260 `"03260"', add
label define citypop_lbl 03280 `"03280"', add
label define citypop_lbl 03281 `"03281"', add
label define citypop_lbl 03292 `"03292"', add
label define citypop_lbl 03305 `"03305"', add
label define citypop_lbl 03313 `"03313"', add
label define citypop_lbl 03318 `"03318"', add
label define citypop_lbl 03323 `"03323"', add
label define citypop_lbl 03325 `"03325"', add
label define citypop_lbl 03332 `"03332"', add
label define citypop_lbl 03344 `"03344"', add
label define citypop_lbl 03346 `"03346"', add
label define citypop_lbl 03380 `"03380"', add
label define citypop_lbl 03393 `"03393"', add
label define citypop_lbl 03409 `"03409"', add
label define citypop_lbl 03428 `"03428"', add
label define citypop_lbl 03440 `"03440"', add
label define citypop_lbl 03443 `"03443"', add
label define citypop_lbl 03456 `"03456"', add
label define citypop_lbl 03469 `"03469"', add
label define citypop_lbl 03472 `"03472"', add
label define citypop_lbl 03482 `"03482"', add
label define citypop_lbl 03526 `"03526"', add
label define citypop_lbl 03542 `"03542"', add
label define citypop_lbl 03546 `"03546"', add
label define citypop_lbl 03550 `"03550"', add
label define citypop_lbl 03563 `"03563"', add
label define citypop_lbl 03577 `"03577"', add
label define citypop_lbl 03579 `"03579"', add
label define citypop_lbl 03585 `"03585"', add
label define citypop_lbl 03609 `"03609"', add
label define citypop_lbl 03613 `"03613"', add
label define citypop_lbl 03650 `"03650"', add
label define citypop_lbl 03653 `"03653"', add
label define citypop_lbl 03664 `"03664"', add
label define citypop_lbl 03673 `"03673"', add
label define citypop_lbl 03683 `"03683"', add
label define citypop_lbl 03684 `"03684"', add
label define citypop_lbl 03691 `"03691"', add
label define citypop_lbl 03694 `"03694"', add
label define citypop_lbl 03699 `"03699"', add
label define citypop_lbl 03710 `"03710"', add
label define citypop_lbl 03722 `"03722"', add
label define citypop_lbl 03727 `"03727"', add
label define citypop_lbl 03728 `"03728"', add
label define citypop_lbl 03736 `"03736"', add
label define citypop_lbl 03759 `"03759"', add
label define citypop_lbl 03815 `"03815"', add
label define citypop_lbl 03826 `"03826"', add
label define citypop_lbl 03829 `"03829"', add
label define citypop_lbl 03845 `"03845"', add
label define citypop_lbl 03852 `"03852"', add
label define citypop_lbl 03855 `"03855"', add
label define citypop_lbl 03870 `"03870"', add
label define citypop_lbl 03890 `"03890"', add
label define citypop_lbl 03930 `"03930"', add
label define citypop_lbl 03931 `"03931"', add
label define citypop_lbl 03940 `"03940"', add
label define citypop_lbl 03960 `"03960"', add
label define citypop_lbl 03967 `"03967"', add
label define citypop_lbl 03992 `"03992"', add
label define citypop_lbl 04032 `"04032"', add
label define citypop_lbl 04070 `"04070"', add
label define citypop_lbl 04084 `"04084"', add
label define citypop_lbl 04169 `"04169"', add
label define citypop_lbl 04239 `"04239"', add
label define citypop_lbl 04250 `"04250"', add
label define citypop_lbl 04253 `"04253"', add
label define citypop_lbl 04258 `"04258"', add
label define citypop_lbl 04272 `"04272"', add
label define citypop_lbl 04277 `"04277"', add
label define citypop_lbl 04294 `"04294"', add
label define citypop_lbl 04345 `"04345"', add
label define citypop_lbl 04356 `"04356"', add
label define citypop_lbl 04370 `"04370"', add
label define citypop_lbl 04373 `"04373"', add
label define citypop_lbl 04415 `"04415"', add
label define citypop_lbl 04424 `"04424"', add
label define citypop_lbl 04431 `"04431"', add
label define citypop_lbl 04443 `"04443"', add
label define citypop_lbl 04447 `"04447"', add
label define citypop_lbl 04473 `"04473"', add
label define citypop_lbl 04476 `"04476"', add
label define citypop_lbl 04482 `"04482"', add
label define citypop_lbl 04506 `"04506"', add
label define citypop_lbl 04517 `"04517"', add
label define citypop_lbl 04522 `"04522"', add
label define citypop_lbl 04531 `"04531"', add
label define citypop_lbl 04538 `"04538"', add
label define citypop_lbl 04556 `"04556"', add
label define citypop_lbl 04557 `"04557"', add
label define citypop_lbl 04566 `"04566"', add
label define citypop_lbl 04615 `"04615"', add
label define citypop_lbl 04656 `"04656"', add
label define citypop_lbl 04667 `"04667"', add
label define citypop_lbl 04676 `"04676"', add
label define citypop_lbl 04725 `"04725"', add
label define citypop_lbl 04743 `"04743"', add
label define citypop_lbl 04784 `"04784"', add
label define citypop_lbl 04847 `"04847"', add
label define citypop_lbl 04927 `"04927"', add
label define citypop_lbl 04938 `"04938"', add
label define citypop_lbl 04945 `"04945"', add
label define citypop_lbl 04969 `"04969"', add
label define citypop_lbl 05040 `"05040"', add
label define citypop_lbl 05056 `"05056"', add
label define citypop_lbl 05061 `"05061"', add
label define citypop_lbl 05163 `"05163"', add
label define citypop_lbl 05310 `"05310"', add
label define citypop_lbl 05377 `"05377"', add
label define citypop_lbl 05402 `"05402"', add
label define citypop_lbl 05408 `"05408"', add
label define citypop_lbl 05409 `"05409"', add
label define citypop_lbl 05559 `"05559"', add
label define citypop_lbl 05575 `"05575"', add
label define citypop_lbl 05630 `"05630"', add
label define citypop_lbl 05634 `"05634"', add
label define citypop_lbl 05649 `"05649"', add
label define citypop_lbl 05704 `"05704"', add
label define citypop_lbl 05721 `"05721"', add
label define citypop_lbl 05734 `"05734"', add
label define citypop_lbl 05738 `"05738"', add
label define citypop_lbl 05743 `"05743"', add
label define citypop_lbl 05759 `"05759"', add
label define citypop_lbl 05763 `"05763"', add
label define citypop_lbl 05801 `"05801"', add
label define citypop_lbl 05815 `"05815"', add
label define citypop_lbl 05820 `"05820"', add
label define citypop_lbl 05825 `"05825"', add
label define citypop_lbl 05875 `"05875"', add
label define citypop_lbl 05891 `"05891"', add
label define citypop_lbl 05908 `"05908"', add
label define citypop_lbl 05962 `"05962"', add
label define citypop_lbl 05966 `"05966"', add
label define citypop_lbl 05970 `"05970"', add
label define citypop_lbl 06069 `"06069"', add
label define citypop_lbl 06103 `"06103"', add
label define citypop_lbl 06161 `"06161"', add
label define citypop_lbl 06281 `"06281"', add
label define citypop_lbl 06294 `"06294"', add
label define citypop_lbl 06305 `"06305"', add
label define citypop_lbl 06314 `"06314"', add
label define citypop_lbl 06362 `"06362"', add
label define citypop_lbl 06364 `"06364"', add
label define citypop_lbl 06374 `"06374"', add
label define citypop_lbl 06383 `"06383"', add
label define citypop_lbl 06458 `"06458"', add
label define citypop_lbl 06464 `"06464"', add
label define citypop_lbl 06501 `"06501"', add
label define citypop_lbl 06512 `"06512"', add
label define citypop_lbl 06631 `"06631"', add
label define citypop_lbl 06699 `"06699"', add
label define citypop_lbl 06709 `"06709"', add
label define citypop_lbl 06738 `"06738"', add
label define citypop_lbl 06768 `"06768"', add
label define citypop_lbl 06790 `"06790"', add
label define citypop_lbl 07007 `"07007"', add
label define citypop_lbl 07240 `"07240"', add
label define citypop_lbl 07360 `"07360"', add
label define citypop_lbl 07410 `"07410"', add
label define citypop_lbl 07440 `"07440"', add
label define citypop_lbl 07708 `"07708"', add
label define citypop_lbl 07767 `"07767"', add
label define citypop_lbl 07822 `"07822"', add
label define citypop_lbl 07858 `"07858"', add
label define citypop_lbl 07868 `"07868"', add
label define citypop_lbl 07897 `"07897"', add
label define citypop_lbl 08014 `"08014"', add
label define citypop_lbl 08022 `"08022"', add
label define citypop_lbl 08160 `"08160"', add
label define citypop_lbl 08568 `"08568"', add
label define citypop_lbl 08591 `"08591"', add
label define citypop_lbl 08711 `"08711"', add
label define citypop_lbl 08755 `"08755"', add
label define citypop_lbl 08783 `"08783"', add
label define citypop_lbl 08835 `"08835"', add
label define citypop_lbl 09041 `"09041"', add
label define citypop_lbl 09148 `"09148"', add
label define citypop_lbl 09359 `"09359"', add
label define citypop_lbl 09497 `"09497"', add
label define citypop_lbl 09513 `"09513"', add
label define citypop_lbl 09533 `"09533"', add
label define citypop_lbl 09834 `"09834"', add
label define citypop_lbl 10069 `"10069"', add
label define citypop_lbl 10280 `"10280"', add
label define citypop_lbl 11105 `"11105"', add
label define citypop_lbl 12033 `"12033"', add
label define citypop_lbl 14484 `"14484"', add
label define citypop_lbl 14564 `"14564"', add
label define citypop_lbl 14698 `"14698"', add
label define citypop_lbl 15043 `"15043"', add
label define citypop_lbl 15176 `"15176"', add
label define citypop_lbl 15856 `"15856"', add
label define citypop_lbl 15952 `"15952"', add
label define citypop_lbl 16235 `"16235"', add
label define citypop_lbl 16306 `"16306"', add
label define citypop_lbl 16882 `"16882"', add
label define citypop_lbl 18496 `"18496"', add
label define citypop_lbl 19313 `"19313"', add
label define citypop_lbl 19704 `"19704"', add
label define citypop_lbl 20716 `"20716"', add
label define citypop_lbl 27837 `"27837"', add
label define citypop_lbl 28333 `"28333"', add
label define citypop_lbl 28428 `"28428"', add
label define citypop_lbl 28960 `"28960"', add
label define citypop_lbl 29669 `"29669"', add
label define citypop_lbl 30051 `"30051"', add
label define citypop_lbl 33968 `"33968"', add
label define citypop_lbl 34853 `"34853"', add
label define citypop_lbl 36210 `"36210"', add
label define citypop_lbl 36948 `"36948"', add
label define citypop_lbl 38471 `"38471"', add
label define citypop_lbl 38494 `"38494"', add
label define citypop_lbl 67109 `"67109"', add
label define citypop_lbl 70716 `"70716"', add
label define citypop_lbl 71367 `"71367"', add
label define citypop_lbl 73226 `"73226"', add
label define citypop_lbl 80083 `"80083"', add
label define citypop_lbl 82138 `"82138"', add
label define citypop_lbl 82144 `"82144"', add
label values citypop citypop_lbl

label define pumasupr_lbl 01100 `"01100"'
label define pumasupr_lbl 01200 `"01200"', add
label define pumasupr_lbl 01300 `"01300"', add
label define pumasupr_lbl 01400 `"01400"', add
label define pumasupr_lbl 01500 `"01500"', add
label define pumasupr_lbl 01600 `"01600"', add
label define pumasupr_lbl 01701 `"01701"', add
label define pumasupr_lbl 01702 `"01702"', add
label define pumasupr_lbl 02100 `"02100"', add
label define pumasupr_lbl 04100 `"04100"', add
label define pumasupr_lbl 04200 `"04200"', add
label define pumasupr_lbl 04301 `"04301"', add
label define pumasupr_lbl 04302 `"04302"', add
label define pumasupr_lbl 04303 `"04303"', add
label define pumasupr_lbl 04304 `"04304"', add
label define pumasupr_lbl 04305 `"04305"', add
label define pumasupr_lbl 04306 `"04306"', add
label define pumasupr_lbl 04401 `"04401"', add
label define pumasupr_lbl 04402 `"04402"', add
label define pumasupr_lbl 05100 `"05100"', add
label define pumasupr_lbl 05200 `"05200"', add
label define pumasupr_lbl 05300 `"05300"', add
label define pumasupr_lbl 05400 `"05400"', add
label define pumasupr_lbl 05500 `"05500"', add
label define pumasupr_lbl 06010 `"06010"', add
label define pumasupr_lbl 06020 `"06020"', add
label define pumasupr_lbl 06030 `"06030"', add
label define pumasupr_lbl 06040 `"06040"', add
label define pumasupr_lbl 06050 `"06050"', add
label define pumasupr_lbl 06060 `"06060"', add
label define pumasupr_lbl 06071 `"06071"', add
label define pumasupr_lbl 06072 `"06072"', add
label define pumasupr_lbl 06080 `"06080"', add
label define pumasupr_lbl 06090 `"06090"', add
label define pumasupr_lbl 06100 `"06100"', add
label define pumasupr_lbl 06110 `"06110"', add
label define pumasupr_lbl 06121 `"06121"', add
label define pumasupr_lbl 06122 `"06122"', add
label define pumasupr_lbl 06130 `"06130"', add
label define pumasupr_lbl 06140 `"06140"', add
label define pumasupr_lbl 06151 `"06151"', add
label define pumasupr_lbl 06152 `"06152"', add
label define pumasupr_lbl 06153 `"06153"', add
label define pumasupr_lbl 06161 `"06161"', add
label define pumasupr_lbl 06162 `"06162"', add
label define pumasupr_lbl 06163 `"06163"', add
label define pumasupr_lbl 06170 `"06170"', add
label define pumasupr_lbl 06180 `"06180"', add
label define pumasupr_lbl 06190 `"06190"', add
label define pumasupr_lbl 06201 `"06201"', add
label define pumasupr_lbl 06202 `"06202"', add
label define pumasupr_lbl 06203 `"06203"', add
label define pumasupr_lbl 06210 `"06210"', add
label define pumasupr_lbl 06220 `"06220"', add
label define pumasupr_lbl 06230 `"06230"', add
label define pumasupr_lbl 06301 `"06301"', add
label define pumasupr_lbl 06302 `"06302"', add
label define pumasupr_lbl 06303 `"06303"', add
label define pumasupr_lbl 06304 `"06304"', add
label define pumasupr_lbl 06305 `"06305"', add
label define pumasupr_lbl 06306 `"06306"', add
label define pumasupr_lbl 06307 `"06307"', add
label define pumasupr_lbl 06401 `"06401"', add
label define pumasupr_lbl 06402 `"06402"', add
label define pumasupr_lbl 06403 `"06403"', add
label define pumasupr_lbl 06404 `"06404"', add
label define pumasupr_lbl 06405 `"06405"', add
label define pumasupr_lbl 06406 `"06406"', add
label define pumasupr_lbl 06407 `"06407"', add
label define pumasupr_lbl 06408 `"06408"', add
label define pumasupr_lbl 06409 `"06409"', add
label define pumasupr_lbl 06410 `"06410"', add
label define pumasupr_lbl 06411 `"06411"', add
label define pumasupr_lbl 06501 `"06501"', add
label define pumasupr_lbl 06502 `"06502"', add
label define pumasupr_lbl 06503 `"06503"', add
label define pumasupr_lbl 06504 `"06504"', add
label define pumasupr_lbl 06505 `"06505"', add
label define pumasupr_lbl 06601 `"06601"', add
label define pumasupr_lbl 06602 `"06602"', add
label define pumasupr_lbl 06603 `"06603"', add
label define pumasupr_lbl 06701 `"06701"', add
label define pumasupr_lbl 06702 `"06702"', add
label define pumasupr_lbl 06703 `"06703"', add
label define pumasupr_lbl 06704 `"06704"', add
label define pumasupr_lbl 06705 `"06705"', add
label define pumasupr_lbl 08101 `"08101"', add
label define pumasupr_lbl 08102 `"08102"', add
label define pumasupr_lbl 08103 `"08103"', add
label define pumasupr_lbl 08104 `"08104"', add
label define pumasupr_lbl 08201 `"08201"', add
label define pumasupr_lbl 08202 `"08202"', add
label define pumasupr_lbl 08203 `"08203"', add
label define pumasupr_lbl 08204 `"08204"', add
label define pumasupr_lbl 08205 `"08205"', add
label define pumasupr_lbl 09100 `"09100"', add
label define pumasupr_lbl 09200 `"09200"', add
label define pumasupr_lbl 09300 `"09300"', add
label define pumasupr_lbl 09400 `"09400"', add
label define pumasupr_lbl 09500 `"09500"', add
label define pumasupr_lbl 09600 `"09600"', add
label define pumasupr_lbl 10100 `"10100"', add
label define pumasupr_lbl 11100 `"11100"', add
label define pumasupr_lbl 12010 `"12010"', add
label define pumasupr_lbl 12020 `"12020"', add
label define pumasupr_lbl 12030 `"12030"', add
label define pumasupr_lbl 12040 `"12040"', add
label define pumasupr_lbl 12051 `"12051"', add
label define pumasupr_lbl 12052 `"12052"', add
label define pumasupr_lbl 12060 `"12060"', add
label define pumasupr_lbl 12070 `"12070"', add
label define pumasupr_lbl 12081 `"12081"', add
label define pumasupr_lbl 12082 `"12082"', add
label define pumasupr_lbl 12083 `"12083"', add
label define pumasupr_lbl 12084 `"12084"', add
label define pumasupr_lbl 12085 `"12085"', add
label define pumasupr_lbl 12091 `"12091"', add
label define pumasupr_lbl 12092 `"12092"', add
label define pumasupr_lbl 12093 `"12093"', add
label define pumasupr_lbl 12100 `"12100"', add
label define pumasupr_lbl 12110 `"12110"', add
label define pumasupr_lbl 12120 `"12120"', add
label define pumasupr_lbl 12130 `"12130"', add
label define pumasupr_lbl 12140 `"12140"', add
label define pumasupr_lbl 12150 `"12150"', add
label define pumasupr_lbl 12161 `"12161"', add
label define pumasupr_lbl 12162 `"12162"', add
label define pumasupr_lbl 12171 `"12171"', add
label define pumasupr_lbl 12172 `"12172"', add
label define pumasupr_lbl 12173 `"12173"', add
label define pumasupr_lbl 12181 `"12181"', add
label define pumasupr_lbl 12182 `"12182"', add
label define pumasupr_lbl 12183 `"12183"', add
label define pumasupr_lbl 12184 `"12184"', add
label define pumasupr_lbl 12185 `"12185"', add
label define pumasupr_lbl 13010 `"13010"', add
label define pumasupr_lbl 13020 `"13020"', add
label define pumasupr_lbl 13030 `"13030"', add
label define pumasupr_lbl 13040 `"13040"', add
label define pumasupr_lbl 13050 `"13050"', add
label define pumasupr_lbl 13060 `"13060"', add
label define pumasupr_lbl 13070 `"13070"', add
label define pumasupr_lbl 13080 `"13080"', add
label define pumasupr_lbl 13090 `"13090"', add
label define pumasupr_lbl 13100 `"13100"', add
label define pumasupr_lbl 13110 `"13110"', add
label define pumasupr_lbl 13120 `"13120"', add
label define pumasupr_lbl 13130 `"13130"', add
label define pumasupr_lbl 13140 `"13140"', add
label define pumasupr_lbl 13150 `"13150"', add
label define pumasupr_lbl 15101 `"15101"', add
label define pumasupr_lbl 15102 `"15102"', add
label define pumasupr_lbl 16100 `"16100"', add
label define pumasupr_lbl 16200 `"16200"', add
label define pumasupr_lbl 16300 `"16300"', add
label define pumasupr_lbl 17010 `"17010"', add
label define pumasupr_lbl 17020 `"17020"', add
label define pumasupr_lbl 17030 `"17030"', add
label define pumasupr_lbl 17040 `"17040"', add
label define pumasupr_lbl 17050 `"17050"', add
label define pumasupr_lbl 17060 `"17060"', add
label define pumasupr_lbl 17070 `"17070"', add
label define pumasupr_lbl 17080 `"17080"', add
label define pumasupr_lbl 17090 `"17090"', add
label define pumasupr_lbl 17100 `"17100"', add
label define pumasupr_lbl 17201 `"17201"', add
label define pumasupr_lbl 17202 `"17202"', add
label define pumasupr_lbl 17300 `"17300"', add
label define pumasupr_lbl 17401 `"17401"', add
label define pumasupr_lbl 17402 `"17402"', add
label define pumasupr_lbl 17403 `"17403"', add
label define pumasupr_lbl 17404 `"17404"', add
label define pumasupr_lbl 17405 `"17405"', add
label define pumasupr_lbl 17501 `"17501"', add
label define pumasupr_lbl 17502 `"17502"', add
label define pumasupr_lbl 17503 `"17503"', add
label define pumasupr_lbl 17504 `"17504"', add
label define pumasupr_lbl 17505 `"17505"', add
label define pumasupr_lbl 18010 `"18010"', add
label define pumasupr_lbl 18020 `"18020"', add
label define pumasupr_lbl 18030 `"18030"', add
label define pumasupr_lbl 18040 `"18040"', add
label define pumasupr_lbl 18050 `"18050"', add
label define pumasupr_lbl 18060 `"18060"', add
label define pumasupr_lbl 18070 `"18070"', add
label define pumasupr_lbl 18080 `"18080"', add
label define pumasupr_lbl 18091 `"18091"', add
label define pumasupr_lbl 18092 `"18092"', add
label define pumasupr_lbl 18100 `"18100"', add
label define pumasupr_lbl 18110 `"18110"', add
label define pumasupr_lbl 19100 `"19100"', add
label define pumasupr_lbl 19200 `"19200"', add
label define pumasupr_lbl 19300 `"19300"', add
label define pumasupr_lbl 19400 `"19400"', add
label define pumasupr_lbl 19500 `"19500"', add
label define pumasupr_lbl 20100 `"20100"', add
label define pumasupr_lbl 20200 `"20200"', add
label define pumasupr_lbl 20300 `"20300"', add
label define pumasupr_lbl 20400 `"20400"', add
label define pumasupr_lbl 20500 `"20500"', add
label define pumasupr_lbl 21100 `"21100"', add
label define pumasupr_lbl 21200 `"21200"', add
label define pumasupr_lbl 21300 `"21300"', add
label define pumasupr_lbl 21400 `"21400"', add
label define pumasupr_lbl 21500 `"21500"', add
label define pumasupr_lbl 21600 `"21600"', add
label define pumasupr_lbl 21700 `"21700"', add
label define pumasupr_lbl 22100 `"22100"', add
label define pumasupr_lbl 22200 `"22200"', add
label define pumasupr_lbl 22300 `"22300"', add
label define pumasupr_lbl 22400 `"22400"', add
label define pumasupr_lbl 22500 `"22500"', add
label define pumasupr_lbl 22600 `"22600"', add
label define pumasupr_lbl 22701 `"22701"', add
label define pumasupr_lbl 22702 `"22702"', add
label define pumasupr_lbl 22800 `"22800"', add
label define pumasupr_lbl 22999 `"22999"', add
label define pumasupr_lbl 23100 `"23100"', add
label define pumasupr_lbl 23200 `"23200"', add
label define pumasupr_lbl 24100 `"24100"', add
label define pumasupr_lbl 24201 `"24201"', add
label define pumasupr_lbl 24202 `"24202"', add
label define pumasupr_lbl 24300 `"24300"', add
label define pumasupr_lbl 24401 `"24401"', add
label define pumasupr_lbl 24402 `"24402"', add
label define pumasupr_lbl 24403 `"24403"', add
label define pumasupr_lbl 24404 `"24404"', add
label define pumasupr_lbl 24501 `"24501"', add
label define pumasupr_lbl 24502 `"24502"', add
label define pumasupr_lbl 25010 `"25010"', add
label define pumasupr_lbl 25020 `"25020"', add
label define pumasupr_lbl 25030 `"25030"', add
label define pumasupr_lbl 25040 `"25040"', add
label define pumasupr_lbl 25050 `"25050"', add
label define pumasupr_lbl 25060 `"25060"', add
label define pumasupr_lbl 25070 `"25070"', add
label define pumasupr_lbl 25080 `"25080"', add
label define pumasupr_lbl 25090 `"25090"', add
label define pumasupr_lbl 25100 `"25100"', add
label define pumasupr_lbl 25110 `"25110"', add
label define pumasupr_lbl 25120 `"25120"', add
label define pumasupr_lbl 25130 `"25130"', add
label define pumasupr_lbl 26010 `"26010"', add
label define pumasupr_lbl 26020 `"26020"', add
label define pumasupr_lbl 26030 `"26030"', add
label define pumasupr_lbl 26040 `"26040"', add
label define pumasupr_lbl 26051 `"26051"', add
label define pumasupr_lbl 26052 `"26052"', add
label define pumasupr_lbl 26060 `"26060"', add
label define pumasupr_lbl 26070 `"26070"', add
label define pumasupr_lbl 26080 `"26080"', add
label define pumasupr_lbl 26090 `"26090"', add
label define pumasupr_lbl 26100 `"26100"', add
label define pumasupr_lbl 26110 `"26110"', add
label define pumasupr_lbl 26121 `"26121"', add
label define pumasupr_lbl 26122 `"26122"', add
label define pumasupr_lbl 26123 `"26123"', add
label define pumasupr_lbl 26124 `"26124"', add
label define pumasupr_lbl 26131 `"26131"', add
label define pumasupr_lbl 26132 `"26132"', add
label define pumasupr_lbl 26133 `"26133"', add
label define pumasupr_lbl 26134 `"26134"', add
label define pumasupr_lbl 27100 `"27100"', add
label define pumasupr_lbl 27200 `"27200"', add
label define pumasupr_lbl 27300 `"27300"', add
label define pumasupr_lbl 27400 `"27400"', add
label define pumasupr_lbl 27500 `"27500"', add
label define pumasupr_lbl 27600 `"27600"', add
label define pumasupr_lbl 27710 `"27710"', add
label define pumasupr_lbl 27720 `"27720"', add
label define pumasupr_lbl 27800 `"27800"', add
label define pumasupr_lbl 27900 `"27900"', add
label define pumasupr_lbl 28100 `"28100"', add
label define pumasupr_lbl 28200 `"28200"', add
label define pumasupr_lbl 28300 `"28300"', add
label define pumasupr_lbl 28400 `"28400"', add
label define pumasupr_lbl 28500 `"28500"', add
label define pumasupr_lbl 28600 `"28600"', add
label define pumasupr_lbl 29100 `"29100"', add
label define pumasupr_lbl 29200 `"29200"', add
label define pumasupr_lbl 29300 `"29300"', add
label define pumasupr_lbl 29400 `"29400"', add
label define pumasupr_lbl 29500 `"29500"', add
label define pumasupr_lbl 29600 `"29600"', add
label define pumasupr_lbl 29701 `"29701"', add
label define pumasupr_lbl 29702 `"29702"', add
label define pumasupr_lbl 29800 `"29800"', add
label define pumasupr_lbl 29900 `"29900"', add
label define pumasupr_lbl 30100 `"30100"', add
label define pumasupr_lbl 30200 `"30200"', add
label define pumasupr_lbl 31100 `"31100"', add
label define pumasupr_lbl 31201 `"31201"', add
label define pumasupr_lbl 31202 `"31202"', add
label define pumasupr_lbl 32100 `"32100"', add
label define pumasupr_lbl 32201 `"32201"', add
label define pumasupr_lbl 32202 `"32202"', add
label define pumasupr_lbl 32203 `"32203"', add
label define pumasupr_lbl 33100 `"33100"', add
label define pumasupr_lbl 33200 `"33200"', add
label define pumasupr_lbl 34011 `"34011"', add
label define pumasupr_lbl 34012 `"34012"', add
label define pumasupr_lbl 34020 `"34020"', add
label define pumasupr_lbl 34030 `"34030"', add
label define pumasupr_lbl 34041 `"34041"', add
label define pumasupr_lbl 34042 `"34042"', add
label define pumasupr_lbl 34050 `"34050"', add
label define pumasupr_lbl 34060 `"34060"', add
label define pumasupr_lbl 34070 `"34070"', add
label define pumasupr_lbl 34080 `"34080"', add
label define pumasupr_lbl 34090 `"34090"', add
label define pumasupr_lbl 34101 `"34101"', add
label define pumasupr_lbl 34102 `"34102"', add
label define pumasupr_lbl 34110 `"34110"', add
label define pumasupr_lbl 34120 `"34120"', add
label define pumasupr_lbl 35100 `"35100"', add
label define pumasupr_lbl 35200 `"35200"', add
label define pumasupr_lbl 35300 `"35300"', add
label define pumasupr_lbl 35400 `"35400"', add
label define pumasupr_lbl 36010 `"36010"', add
label define pumasupr_lbl 36021 `"36021"', add
label define pumasupr_lbl 36022 `"36022"', add
label define pumasupr_lbl 36030 `"36030"', add
label define pumasupr_lbl 36041 `"36041"', add
label define pumasupr_lbl 36042 `"36042"', add
label define pumasupr_lbl 36051 `"36051"', add
label define pumasupr_lbl 36052 `"36052"', add
label define pumasupr_lbl 36060 `"36060"', add
label define pumasupr_lbl 36070 `"36070"', add
label define pumasupr_lbl 36081 `"36081"', add
label define pumasupr_lbl 36082 `"36082"', add
label define pumasupr_lbl 36083 `"36083"', add
label define pumasupr_lbl 36084 `"36084"', add
label define pumasupr_lbl 36085 `"36085"', add
label define pumasupr_lbl 36091 `"36091"', add
label define pumasupr_lbl 36092 `"36092"', add
label define pumasupr_lbl 36101 `"36101"', add
label define pumasupr_lbl 36102 `"36102"', add
label define pumasupr_lbl 36103 `"36103"', add
label define pumasupr_lbl 36111 `"36111"', add
label define pumasupr_lbl 36112 `"36112"', add
label define pumasupr_lbl 36113 `"36113"', add
label define pumasupr_lbl 36114 `"36114"', add
label define pumasupr_lbl 36121 `"36121"', add
label define pumasupr_lbl 36122 `"36122"', add
label define pumasupr_lbl 36123 `"36123"', add
label define pumasupr_lbl 36124 `"36124"', add
label define pumasupr_lbl 36125 `"36125"', add
label define pumasupr_lbl 36130 `"36130"', add
label define pumasupr_lbl 36141 `"36141"', add
label define pumasupr_lbl 36142 `"36142"', add
label define pumasupr_lbl 36143 `"36143"', add
label define pumasupr_lbl 36151 `"36151"', add
label define pumasupr_lbl 36152 `"36152"', add
label define pumasupr_lbl 36153 `"36153"', add
label define pumasupr_lbl 37010 `"37010"', add
label define pumasupr_lbl 37020 `"37020"', add
label define pumasupr_lbl 37030 `"37030"', add
label define pumasupr_lbl 37040 `"37040"', add
label define pumasupr_lbl 37050 `"37050"', add
label define pumasupr_lbl 37060 `"37060"', add
label define pumasupr_lbl 37070 `"37070"', add
label define pumasupr_lbl 37080 `"37080"', add
label define pumasupr_lbl 37090 `"37090"', add
label define pumasupr_lbl 37100 `"37100"', add
label define pumasupr_lbl 37110 `"37110"', add
label define pumasupr_lbl 37120 `"37120"', add
label define pumasupr_lbl 37130 `"37130"', add
label define pumasupr_lbl 37140 `"37140"', add
label define pumasupr_lbl 38100 `"38100"', add
label define pumasupr_lbl 39010 `"39010"', add
label define pumasupr_lbl 39020 `"39020"', add
label define pumasupr_lbl 39030 `"39030"', add
label define pumasupr_lbl 39040 `"39040"', add
label define pumasupr_lbl 39050 `"39050"', add
label define pumasupr_lbl 39061 `"39061"', add
label define pumasupr_lbl 39062 `"39062"', add
label define pumasupr_lbl 39063 `"39063"', add
label define pumasupr_lbl 39070 `"39070"', add
label define pumasupr_lbl 39080 `"39080"', add
label define pumasupr_lbl 39090 `"39090"', add
label define pumasupr_lbl 39100 `"39100"', add
label define pumasupr_lbl 39110 `"39110"', add
label define pumasupr_lbl 39120 `"39120"', add
label define pumasupr_lbl 39130 `"39130"', add
label define pumasupr_lbl 39141 `"39141"', add
label define pumasupr_lbl 39142 `"39142"', add
label define pumasupr_lbl 39150 `"39150"', add
label define pumasupr_lbl 39160 `"39160"', add
label define pumasupr_lbl 39171 `"39171"', add
label define pumasupr_lbl 39172 `"39172"', add
label define pumasupr_lbl 39180 `"39180"', add
label define pumasupr_lbl 40100 `"40100"', add
label define pumasupr_lbl 40201 `"40201"', add
label define pumasupr_lbl 40202 `"40202"', add
label define pumasupr_lbl 40300 `"40300"', add
label define pumasupr_lbl 40400 `"40400"', add
label define pumasupr_lbl 40500 `"40500"', add
label define pumasupr_lbl 41100 `"41100"', add
label define pumasupr_lbl 41200 `"41200"', add
label define pumasupr_lbl 41300 `"41300"', add
label define pumasupr_lbl 41400 `"41400"', add
label define pumasupr_lbl 41501 `"41501"', add
label define pumasupr_lbl 41502 `"41502"', add
label define pumasupr_lbl 41503 `"41503"', add
label define pumasupr_lbl 42010 `"42010"', add
label define pumasupr_lbl 42020 `"42020"', add
label define pumasupr_lbl 42030 `"42030"', add
label define pumasupr_lbl 42040 `"42040"', add
label define pumasupr_lbl 42050 `"42050"', add
label define pumasupr_lbl 42060 `"42060"', add
label define pumasupr_lbl 42071 `"42071"', add
label define pumasupr_lbl 42072 `"42072"', add
label define pumasupr_lbl 42073 `"42073"', add
label define pumasupr_lbl 42080 `"42080"', add
label define pumasupr_lbl 42090 `"42090"', add
label define pumasupr_lbl 42100 `"42100"', add
label define pumasupr_lbl 42110 `"42110"', add
label define pumasupr_lbl 42120 `"42120"', add
label define pumasupr_lbl 42130 `"42130"', add
label define pumasupr_lbl 42140 `"42140"', add
label define pumasupr_lbl 42151 `"42151"', add
label define pumasupr_lbl 42152 `"42152"', add
label define pumasupr_lbl 42153 `"42153"', add
label define pumasupr_lbl 42160 `"42160"', add
label define pumasupr_lbl 42170 `"42170"', add
label define pumasupr_lbl 42180 `"42180"', add
label define pumasupr_lbl 42190 `"42190"', add
label define pumasupr_lbl 44100 `"44100"', add
label define pumasupr_lbl 44200 `"44200"', add
label define pumasupr_lbl 45100 `"45100"', add
label define pumasupr_lbl 45200 `"45200"', add
label define pumasupr_lbl 45300 `"45300"', add
label define pumasupr_lbl 45400 `"45400"', add
label define pumasupr_lbl 45500 `"45500"', add
label define pumasupr_lbl 45600 `"45600"', add
label define pumasupr_lbl 45700 `"45700"', add
label define pumasupr_lbl 45800 `"45800"', add
label define pumasupr_lbl 46100 `"46100"', add
label define pumasupr_lbl 47010 `"47010"', add
label define pumasupr_lbl 47020 `"47020"', add
label define pumasupr_lbl 47030 `"47030"', add
label define pumasupr_lbl 47040 `"47040"', add
label define pumasupr_lbl 47050 `"47050"', add
label define pumasupr_lbl 47060 `"47060"', add
label define pumasupr_lbl 47070 `"47070"', add
label define pumasupr_lbl 47081 `"47081"', add
label define pumasupr_lbl 47082 `"47082"', add
label define pumasupr_lbl 47090 `"47090"', add
label define pumasupr_lbl 47101 `"47101"', add
label define pumasupr_lbl 47102 `"47102"', add
label define pumasupr_lbl 48010 `"48010"', add
label define pumasupr_lbl 48020 `"48020"', add
label define pumasupr_lbl 48030 `"48030"', add
label define pumasupr_lbl 48040 `"48040"', add
label define pumasupr_lbl 48050 `"48050"', add
label define pumasupr_lbl 48060 `"48060"', add
label define pumasupr_lbl 48070 `"48070"', add
label define pumasupr_lbl 48080 `"48080"', add
label define pumasupr_lbl 48090 `"48090"', add
label define pumasupr_lbl 48101 `"48101"', add
label define pumasupr_lbl 48102 `"48102"', add
label define pumasupr_lbl 48103 `"48103"', add
label define pumasupr_lbl 48104 `"48104"', add
label define pumasupr_lbl 48111 `"48111"', add
label define pumasupr_lbl 48112 `"48112"', add
label define pumasupr_lbl 48113 `"48113"', add
label define pumasupr_lbl 48120 `"48120"', add
label define pumasupr_lbl 48130 `"48130"', add
label define pumasupr_lbl 48140 `"48140"', add
label define pumasupr_lbl 48150 `"48150"', add
label define pumasupr_lbl 48160 `"48160"', add
label define pumasupr_lbl 48170 `"48170"', add
label define pumasupr_lbl 48181 `"48181"', add
label define pumasupr_lbl 48182 `"48182"', add
label define pumasupr_lbl 48183 `"48183"', add
label define pumasupr_lbl 48184 `"48184"', add
label define pumasupr_lbl 48185 `"48185"', add
label define pumasupr_lbl 48186 `"48186"', add
label define pumasupr_lbl 48187 `"48187"', add
label define pumasupr_lbl 48190 `"48190"', add
label define pumasupr_lbl 48200 `"48200"', add
label define pumasupr_lbl 48210 `"48210"', add
label define pumasupr_lbl 48221 `"48221"', add
label define pumasupr_lbl 48222 `"48222"', add
label define pumasupr_lbl 48231 `"48231"', add
label define pumasupr_lbl 48232 `"48232"', add
label define pumasupr_lbl 48233 `"48233"', add
label define pumasupr_lbl 48240 `"48240"', add
label define pumasupr_lbl 48250 `"48250"', add
label define pumasupr_lbl 48260 `"48260"', add
label define pumasupr_lbl 48270 `"48270"', add
label define pumasupr_lbl 49100 `"49100"', add
label define pumasupr_lbl 49200 `"49200"', add
label define pumasupr_lbl 49301 `"49301"', add
label define pumasupr_lbl 49302 `"49302"', add
label define pumasupr_lbl 50100 `"50100"', add
label define pumasupr_lbl 51011 `"51011"', add
label define pumasupr_lbl 51012 `"51012"', add
label define pumasupr_lbl 51020 `"51020"', add
label define pumasupr_lbl 51030 `"51030"', add
label define pumasupr_lbl 51040 `"51040"', add
label define pumasupr_lbl 51050 `"51050"', add
label define pumasupr_lbl 51060 `"51060"', add
label define pumasupr_lbl 51070 `"51070"', add
label define pumasupr_lbl 51080 `"51080"', add
label define pumasupr_lbl 51090 `"51090"', add
label define pumasupr_lbl 51100 `"51100"', add
label define pumasupr_lbl 51110 `"51110"', add
label define pumasupr_lbl 51120 `"51120"', add
label define pumasupr_lbl 53010 `"53010"', add
label define pumasupr_lbl 53020 `"53020"', add
label define pumasupr_lbl 53030 `"53030"', add
label define pumasupr_lbl 53040 `"53040"', add
label define pumasupr_lbl 53050 `"53050"', add
label define pumasupr_lbl 53060 `"53060"', add
label define pumasupr_lbl 53070 `"53070"', add
label define pumasupr_lbl 53081 `"53081"', add
label define pumasupr_lbl 53082 `"53082"', add
label define pumasupr_lbl 53090 `"53090"', add
label define pumasupr_lbl 53100 `"53100"', add
label define pumasupr_lbl 54100 `"54100"', add
label define pumasupr_lbl 54200 `"54200"', add
label define pumasupr_lbl 54300 `"54300"', add
label define pumasupr_lbl 55100 `"55100"', add
label define pumasupr_lbl 55200 `"55200"', add
label define pumasupr_lbl 55300 `"55300"', add
label define pumasupr_lbl 55400 `"55400"', add
label define pumasupr_lbl 55500 `"55500"', add
label define pumasupr_lbl 55600 `"55600"', add
label define pumasupr_lbl 55700 `"55700"', add
label define pumasupr_lbl 55800 `"55800"', add
label define pumasupr_lbl 55900 `"55900"', add
label define pumasupr_lbl 56100 `"56100"', add
label define pumasupr_lbl 72100 `"72100"', add
label define pumasupr_lbl 72200 `"72200"', add
label define pumasupr_lbl 72300 `"72300"', add
label define pumasupr_lbl 72400 `"72400"', add
label define pumasupr_lbl 72500 `"72500"', add
label define pumasupr_lbl 72600 `"72600"', add
label define pumasupr_lbl 72700 `"72700"', add
label define pumasupr_lbl 72800 `"72800"', add
label values pumasupr pumasupr_lbl

label define conspuma_lbl 001 `"001"'
label define conspuma_lbl 002 `"002"', add
label define conspuma_lbl 003 `"003"', add
label define conspuma_lbl 004 `"004"', add
label define conspuma_lbl 005 `"005"', add
label define conspuma_lbl 006 `"006"', add
label define conspuma_lbl 007 `"007"', add
label define conspuma_lbl 008 `"008"', add
label define conspuma_lbl 009 `"009"', add
label define conspuma_lbl 010 `"010"', add
label define conspuma_lbl 011 `"011"', add
label define conspuma_lbl 012 `"012"', add
label define conspuma_lbl 013 `"013"', add
label define conspuma_lbl 014 `"014"', add
label define conspuma_lbl 015 `"015"', add
label define conspuma_lbl 016 `"016"', add
label define conspuma_lbl 017 `"017"', add
label define conspuma_lbl 018 `"018"', add
label define conspuma_lbl 019 `"019"', add
label define conspuma_lbl 020 `"020"', add
label define conspuma_lbl 021 `"021"', add
label define conspuma_lbl 022 `"022"', add
label define conspuma_lbl 023 `"023"', add
label define conspuma_lbl 024 `"024"', add
label define conspuma_lbl 025 `"025"', add
label define conspuma_lbl 026 `"026"', add
label define conspuma_lbl 027 `"027"', add
label define conspuma_lbl 028 `"028"', add
label define conspuma_lbl 029 `"029"', add
label define conspuma_lbl 030 `"030"', add
label define conspuma_lbl 031 `"031"', add
label define conspuma_lbl 032 `"032"', add
label define conspuma_lbl 033 `"033"', add
label define conspuma_lbl 034 `"034"', add
label define conspuma_lbl 035 `"035"', add
label define conspuma_lbl 036 `"036"', add
label define conspuma_lbl 037 `"037"', add
label define conspuma_lbl 038 `"038"', add
label define conspuma_lbl 039 `"039"', add
label define conspuma_lbl 040 `"040"', add
label define conspuma_lbl 041 `"041"', add
label define conspuma_lbl 042 `"042"', add
label define conspuma_lbl 043 `"043"', add
label define conspuma_lbl 044 `"044"', add
label define conspuma_lbl 045 `"045"', add
label define conspuma_lbl 046 `"046"', add
label define conspuma_lbl 047 `"047"', add
label define conspuma_lbl 048 `"048"', add
label define conspuma_lbl 049 `"049"', add
label define conspuma_lbl 050 `"050"', add
label define conspuma_lbl 051 `"051"', add
label define conspuma_lbl 052 `"052"', add
label define conspuma_lbl 053 `"053"', add
label define conspuma_lbl 054 `"054"', add
label define conspuma_lbl 055 `"055"', add
label define conspuma_lbl 056 `"056"', add
label define conspuma_lbl 057 `"057"', add
label define conspuma_lbl 058 `"058"', add
label define conspuma_lbl 059 `"059"', add
label define conspuma_lbl 060 `"060"', add
label define conspuma_lbl 061 `"061"', add
label define conspuma_lbl 062 `"062"', add
label define conspuma_lbl 063 `"063"', add
label define conspuma_lbl 064 `"064"', add
label define conspuma_lbl 065 `"065"', add
label define conspuma_lbl 066 `"066"', add
label define conspuma_lbl 067 `"067"', add
label define conspuma_lbl 068 `"068"', add
label define conspuma_lbl 069 `"069"', add
label define conspuma_lbl 070 `"070"', add
label define conspuma_lbl 071 `"071"', add
label define conspuma_lbl 072 `"072"', add
label define conspuma_lbl 073 `"073"', add
label define conspuma_lbl 074 `"074"', add
label define conspuma_lbl 075 `"075"', add
label define conspuma_lbl 076 `"076"', add
label define conspuma_lbl 077 `"077"', add
label define conspuma_lbl 078 `"078"', add
label define conspuma_lbl 079 `"079"', add
label define conspuma_lbl 080 `"080"', add
label define conspuma_lbl 081 `"081"', add
label define conspuma_lbl 082 `"082"', add
label define conspuma_lbl 083 `"083"', add
label define conspuma_lbl 084 `"084"', add
label define conspuma_lbl 085 `"085"', add
label define conspuma_lbl 086 `"086"', add
label define conspuma_lbl 087 `"087"', add
label define conspuma_lbl 088 `"088"', add
label define conspuma_lbl 089 `"089"', add
label define conspuma_lbl 090 `"090"', add
label define conspuma_lbl 091 `"091"', add
label define conspuma_lbl 092 `"092"', add
label define conspuma_lbl 093 `"093"', add
label define conspuma_lbl 094 `"094"', add
label define conspuma_lbl 095 `"095"', add
label define conspuma_lbl 096 `"096"', add
label define conspuma_lbl 097 `"097"', add
label define conspuma_lbl 098 `"098"', add
label define conspuma_lbl 099 `"099"', add
label define conspuma_lbl 100 `"100"', add
label define conspuma_lbl 101 `"101"', add
label define conspuma_lbl 102 `"102"', add
label define conspuma_lbl 103 `"103"', add
label define conspuma_lbl 104 `"104"', add
label define conspuma_lbl 105 `"105"', add
label define conspuma_lbl 106 `"106"', add
label define conspuma_lbl 107 `"107"', add
label define conspuma_lbl 108 `"108"', add
label define conspuma_lbl 109 `"109"', add
label define conspuma_lbl 110 `"110"', add
label define conspuma_lbl 111 `"111"', add
label define conspuma_lbl 112 `"112"', add
label define conspuma_lbl 113 `"113"', add
label define conspuma_lbl 114 `"114"', add
label define conspuma_lbl 115 `"115"', add
label define conspuma_lbl 116 `"116"', add
label define conspuma_lbl 117 `"117"', add
label define conspuma_lbl 118 `"118"', add
label define conspuma_lbl 119 `"119"', add
label define conspuma_lbl 120 `"120"', add
label define conspuma_lbl 121 `"121"', add
label define conspuma_lbl 122 `"122"', add
label define conspuma_lbl 123 `"123"', add
label define conspuma_lbl 124 `"124"', add
label define conspuma_lbl 125 `"125"', add
label define conspuma_lbl 126 `"126"', add
label define conspuma_lbl 127 `"127"', add
label define conspuma_lbl 128 `"128"', add
label define conspuma_lbl 129 `"129"', add
label define conspuma_lbl 130 `"130"', add
label define conspuma_lbl 131 `"131"', add
label define conspuma_lbl 132 `"132"', add
label define conspuma_lbl 133 `"133"', add
label define conspuma_lbl 134 `"134"', add
label define conspuma_lbl 135 `"135"', add
label define conspuma_lbl 136 `"136"', add
label define conspuma_lbl 137 `"137"', add
label define conspuma_lbl 138 `"138"', add
label define conspuma_lbl 139 `"139"', add
label define conspuma_lbl 140 `"140"', add
label define conspuma_lbl 141 `"141"', add
label define conspuma_lbl 142 `"142"', add
label define conspuma_lbl 143 `"143"', add
label define conspuma_lbl 144 `"144"', add
label define conspuma_lbl 145 `"145"', add
label define conspuma_lbl 146 `"146"', add
label define conspuma_lbl 147 `"147"', add
label define conspuma_lbl 148 `"148"', add
label define conspuma_lbl 149 `"149"', add
label define conspuma_lbl 150 `"150"', add
label define conspuma_lbl 151 `"151"', add
label define conspuma_lbl 152 `"152"', add
label define conspuma_lbl 153 `"153"', add
label define conspuma_lbl 154 `"154"', add
label define conspuma_lbl 155 `"155"', add
label define conspuma_lbl 156 `"156"', add
label define conspuma_lbl 157 `"157"', add
label define conspuma_lbl 158 `"158"', add
label define conspuma_lbl 159 `"159"', add
label define conspuma_lbl 160 `"160"', add
label define conspuma_lbl 161 `"161"', add
label define conspuma_lbl 162 `"162"', add
label define conspuma_lbl 163 `"163"', add
label define conspuma_lbl 164 `"164"', add
label define conspuma_lbl 165 `"165"', add
label define conspuma_lbl 166 `"166"', add
label define conspuma_lbl 167 `"167"', add
label define conspuma_lbl 168 `"168"', add
label define conspuma_lbl 169 `"169"', add
label define conspuma_lbl 170 `"170"', add
label define conspuma_lbl 171 `"171"', add
label define conspuma_lbl 172 `"172"', add
label define conspuma_lbl 173 `"173"', add
label define conspuma_lbl 174 `"174"', add
label define conspuma_lbl 175 `"175"', add
label define conspuma_lbl 176 `"176"', add
label define conspuma_lbl 177 `"177"', add
label define conspuma_lbl 178 `"178"', add
label define conspuma_lbl 179 `"179"', add
label define conspuma_lbl 180 `"180"', add
label define conspuma_lbl 181 `"181"', add
label define conspuma_lbl 182 `"182"', add
label define conspuma_lbl 183 `"183"', add
label define conspuma_lbl 184 `"184"', add
label define conspuma_lbl 185 `"185"', add
label define conspuma_lbl 186 `"186"', add
label define conspuma_lbl 187 `"187"', add
label define conspuma_lbl 188 `"188"', add
label define conspuma_lbl 189 `"189"', add
label define conspuma_lbl 190 `"190"', add
label define conspuma_lbl 191 `"191"', add
label define conspuma_lbl 192 `"192"', add
label define conspuma_lbl 193 `"193"', add
label define conspuma_lbl 194 `"194"', add
label define conspuma_lbl 195 `"195"', add
label define conspuma_lbl 196 `"196"', add
label define conspuma_lbl 197 `"197"', add
label define conspuma_lbl 198 `"198"', add
label define conspuma_lbl 199 `"199"', add
label define conspuma_lbl 200 `"200"', add
label define conspuma_lbl 201 `"201"', add
label define conspuma_lbl 202 `"202"', add
label define conspuma_lbl 203 `"203"', add
label define conspuma_lbl 204 `"204"', add
label define conspuma_lbl 205 `"205"', add
label define conspuma_lbl 206 `"206"', add
label define conspuma_lbl 207 `"207"', add
label define conspuma_lbl 208 `"208"', add
label define conspuma_lbl 209 `"209"', add
label define conspuma_lbl 210 `"210"', add
label define conspuma_lbl 211 `"211"', add
label define conspuma_lbl 212 `"212"', add
label define conspuma_lbl 213 `"213"', add
label define conspuma_lbl 214 `"214"', add
label define conspuma_lbl 215 `"215"', add
label define conspuma_lbl 216 `"216"', add
label define conspuma_lbl 217 `"217"', add
label define conspuma_lbl 218 `"218"', add
label define conspuma_lbl 219 `"219"', add
label define conspuma_lbl 220 `"220"', add
label define conspuma_lbl 221 `"221"', add
label define conspuma_lbl 222 `"222"', add
label define conspuma_lbl 223 `"223"', add
label define conspuma_lbl 224 `"224"', add
label define conspuma_lbl 225 `"225"', add
label define conspuma_lbl 226 `"226"', add
label define conspuma_lbl 227 `"227"', add
label define conspuma_lbl 228 `"228"', add
label define conspuma_lbl 229 `"229"', add
label define conspuma_lbl 230 `"230"', add
label define conspuma_lbl 231 `"231"', add
label define conspuma_lbl 232 `"232"', add
label define conspuma_lbl 233 `"233"', add
label define conspuma_lbl 234 `"234"', add
label define conspuma_lbl 235 `"235"', add
label define conspuma_lbl 236 `"236"', add
label define conspuma_lbl 237 `"237"', add
label define conspuma_lbl 238 `"238"', add
label define conspuma_lbl 239 `"239"', add
label define conspuma_lbl 240 `"240"', add
label define conspuma_lbl 241 `"241"', add
label define conspuma_lbl 242 `"242"', add
label define conspuma_lbl 243 `"243"', add
label define conspuma_lbl 244 `"244"', add
label define conspuma_lbl 245 `"245"', add
label define conspuma_lbl 246 `"246"', add
label define conspuma_lbl 247 `"247"', add
label define conspuma_lbl 248 `"248"', add
label define conspuma_lbl 249 `"249"', add
label define conspuma_lbl 250 `"250"', add
label define conspuma_lbl 251 `"251"', add
label define conspuma_lbl 252 `"252"', add
label define conspuma_lbl 253 `"253"', add
label define conspuma_lbl 254 `"254"', add
label define conspuma_lbl 255 `"255"', add
label define conspuma_lbl 256 `"256"', add
label define conspuma_lbl 257 `"257"', add
label define conspuma_lbl 258 `"258"', add
label define conspuma_lbl 259 `"259"', add
label define conspuma_lbl 260 `"260"', add
label define conspuma_lbl 261 `"261"', add
label define conspuma_lbl 262 `"262"', add
label define conspuma_lbl 263 `"263"', add
label define conspuma_lbl 264 `"264"', add
label define conspuma_lbl 265 `"265"', add
label define conspuma_lbl 266 `"266"', add
label define conspuma_lbl 267 `"267"', add
label define conspuma_lbl 268 `"268"', add
label define conspuma_lbl 269 `"269"', add
label define conspuma_lbl 270 `"270"', add
label define conspuma_lbl 271 `"271"', add
label define conspuma_lbl 272 `"272"', add
label define conspuma_lbl 273 `"273"', add
label define conspuma_lbl 274 `"274"', add
label define conspuma_lbl 275 `"275"', add
label define conspuma_lbl 276 `"276"', add
label define conspuma_lbl 277 `"277"', add
label define conspuma_lbl 278 `"278"', add
label define conspuma_lbl 279 `"279"', add
label define conspuma_lbl 280 `"280"', add
label define conspuma_lbl 281 `"281"', add
label define conspuma_lbl 282 `"282"', add
label define conspuma_lbl 283 `"283"', add
label define conspuma_lbl 284 `"284"', add
label define conspuma_lbl 285 `"285"', add
label define conspuma_lbl 286 `"286"', add
label define conspuma_lbl 287 `"287"', add
label define conspuma_lbl 288 `"288"', add
label define conspuma_lbl 289 `"289"', add
label define conspuma_lbl 290 `"290"', add
label define conspuma_lbl 291 `"291"', add
label define conspuma_lbl 292 `"292"', add
label define conspuma_lbl 293 `"293"', add
label define conspuma_lbl 294 `"294"', add
label define conspuma_lbl 295 `"295"', add
label define conspuma_lbl 296 `"296"', add
label define conspuma_lbl 297 `"297"', add
label define conspuma_lbl 298 `"298"', add
label define conspuma_lbl 299 `"299"', add
label define conspuma_lbl 300 `"300"', add
label define conspuma_lbl 301 `"301"', add
label define conspuma_lbl 302 `"302"', add
label define conspuma_lbl 303 `"303"', add
label define conspuma_lbl 304 `"304"', add
label define conspuma_lbl 305 `"305"', add
label define conspuma_lbl 306 `"306"', add
label define conspuma_lbl 307 `"307"', add
label define conspuma_lbl 308 `"308"', add
label define conspuma_lbl 309 `"309"', add
label define conspuma_lbl 310 `"310"', add
label define conspuma_lbl 311 `"311"', add
label define conspuma_lbl 312 `"312"', add
label define conspuma_lbl 313 `"313"', add
label define conspuma_lbl 314 `"314"', add
label define conspuma_lbl 315 `"315"', add
label define conspuma_lbl 316 `"316"', add
label define conspuma_lbl 317 `"317"', add
label define conspuma_lbl 318 `"318"', add
label define conspuma_lbl 319 `"319"', add
label define conspuma_lbl 320 `"320"', add
label define conspuma_lbl 321 `"321"', add
label define conspuma_lbl 322 `"322"', add
label define conspuma_lbl 323 `"323"', add
label define conspuma_lbl 324 `"324"', add
label define conspuma_lbl 325 `"325"', add
label define conspuma_lbl 326 `"326"', add
label define conspuma_lbl 327 `"327"', add
label define conspuma_lbl 328 `"328"', add
label define conspuma_lbl 329 `"329"', add
label define conspuma_lbl 330 `"330"', add
label define conspuma_lbl 331 `"331"', add
label define conspuma_lbl 332 `"332"', add
label define conspuma_lbl 333 `"333"', add
label define conspuma_lbl 334 `"334"', add
label define conspuma_lbl 335 `"335"', add
label define conspuma_lbl 336 `"336"', add
label define conspuma_lbl 337 `"337"', add
label define conspuma_lbl 338 `"338"', add
label define conspuma_lbl 339 `"339"', add
label define conspuma_lbl 340 `"340"', add
label define conspuma_lbl 341 `"341"', add
label define conspuma_lbl 342 `"342"', add
label define conspuma_lbl 343 `"343"', add
label define conspuma_lbl 344 `"344"', add
label define conspuma_lbl 345 `"345"', add
label define conspuma_lbl 346 `"346"', add
label define conspuma_lbl 347 `"347"', add
label define conspuma_lbl 348 `"348"', add
label define conspuma_lbl 349 `"349"', add
label define conspuma_lbl 350 `"350"', add
label define conspuma_lbl 351 `"351"', add
label define conspuma_lbl 352 `"352"', add
label define conspuma_lbl 353 `"353"', add
label define conspuma_lbl 354 `"354"', add
label define conspuma_lbl 355 `"355"', add
label define conspuma_lbl 356 `"356"', add
label define conspuma_lbl 357 `"357"', add
label define conspuma_lbl 358 `"358"', add
label define conspuma_lbl 359 `"359"', add
label define conspuma_lbl 360 `"360"', add
label define conspuma_lbl 361 `"361"', add
label define conspuma_lbl 362 `"362"', add
label define conspuma_lbl 363 `"363"', add
label define conspuma_lbl 364 `"364"', add
label define conspuma_lbl 365 `"365"', add
label define conspuma_lbl 366 `"366"', add
label define conspuma_lbl 367 `"367"', add
label define conspuma_lbl 368 `"368"', add
label define conspuma_lbl 369 `"369"', add
label define conspuma_lbl 370 `"370"', add
label define conspuma_lbl 371 `"371"', add
label define conspuma_lbl 372 `"372"', add
label define conspuma_lbl 373 `"373"', add
label define conspuma_lbl 374 `"374"', add
label define conspuma_lbl 375 `"375"', add
label define conspuma_lbl 376 `"376"', add
label define conspuma_lbl 377 `"377"', add
label define conspuma_lbl 378 `"378"', add
label define conspuma_lbl 379 `"379"', add
label define conspuma_lbl 380 `"380"', add
label define conspuma_lbl 381 `"381"', add
label define conspuma_lbl 382 `"382"', add
label define conspuma_lbl 383 `"383"', add
label define conspuma_lbl 384 `"384"', add
label define conspuma_lbl 385 `"385"', add
label define conspuma_lbl 386 `"386"', add
label define conspuma_lbl 387 `"387"', add
label define conspuma_lbl 388 `"388"', add
label define conspuma_lbl 389 `"389"', add
label define conspuma_lbl 390 `"390"', add
label define conspuma_lbl 391 `"391"', add
label define conspuma_lbl 392 `"392"', add
label define conspuma_lbl 393 `"393"', add
label define conspuma_lbl 394 `"394"', add
label define conspuma_lbl 395 `"395"', add
label define conspuma_lbl 396 `"396"', add
label define conspuma_lbl 397 `"397"', add
label define conspuma_lbl 398 `"398"', add
label define conspuma_lbl 399 `"399"', add
label define conspuma_lbl 400 `"400"', add
label define conspuma_lbl 401 `"401"', add
label define conspuma_lbl 402 `"402"', add
label define conspuma_lbl 403 `"403"', add
label define conspuma_lbl 404 `"404"', add
label define conspuma_lbl 405 `"405"', add
label define conspuma_lbl 406 `"406"', add
label define conspuma_lbl 407 `"407"', add
label define conspuma_lbl 408 `"408"', add
label define conspuma_lbl 409 `"409"', add
label define conspuma_lbl 410 `"410"', add
label define conspuma_lbl 411 `"411"', add
label define conspuma_lbl 412 `"412"', add
label define conspuma_lbl 413 `"413"', add
label define conspuma_lbl 414 `"414"', add
label define conspuma_lbl 415 `"415"', add
label define conspuma_lbl 416 `"416"', add
label define conspuma_lbl 417 `"417"', add
label define conspuma_lbl 418 `"418"', add
label define conspuma_lbl 419 `"419"', add
label define conspuma_lbl 420 `"420"', add
label define conspuma_lbl 421 `"421"', add
label define conspuma_lbl 422 `"422"', add
label define conspuma_lbl 423 `"423"', add
label define conspuma_lbl 424 `"424"', add
label define conspuma_lbl 425 `"425"', add
label define conspuma_lbl 426 `"426"', add
label define conspuma_lbl 427 `"427"', add
label define conspuma_lbl 428 `"428"', add
label define conspuma_lbl 429 `"429"', add
label define conspuma_lbl 430 `"430"', add
label define conspuma_lbl 431 `"431"', add
label define conspuma_lbl 432 `"432"', add
label define conspuma_lbl 433 `"433"', add
label define conspuma_lbl 434 `"434"', add
label define conspuma_lbl 435 `"435"', add
label define conspuma_lbl 436 `"436"', add
label define conspuma_lbl 437 `"437"', add
label define conspuma_lbl 438 `"438"', add
label define conspuma_lbl 439 `"439"', add
label define conspuma_lbl 440 `"440"', add
label define conspuma_lbl 441 `"441"', add
label define conspuma_lbl 442 `"442"', add
label define conspuma_lbl 443 `"443"', add
label define conspuma_lbl 444 `"444"', add
label define conspuma_lbl 445 `"445"', add
label define conspuma_lbl 446 `"446"', add
label define conspuma_lbl 447 `"447"', add
label define conspuma_lbl 448 `"448"', add
label define conspuma_lbl 449 `"449"', add
label define conspuma_lbl 450 `"450"', add
label define conspuma_lbl 451 `"451"', add
label define conspuma_lbl 452 `"452"', add
label define conspuma_lbl 453 `"453"', add
label define conspuma_lbl 454 `"454"', add
label define conspuma_lbl 455 `"455"', add
label define conspuma_lbl 456 `"456"', add
label define conspuma_lbl 457 `"457"', add
label define conspuma_lbl 458 `"458"', add
label define conspuma_lbl 459 `"459"', add
label define conspuma_lbl 460 `"460"', add
label define conspuma_lbl 461 `"461"', add
label define conspuma_lbl 462 `"462"', add
label define conspuma_lbl 463 `"463"', add
label define conspuma_lbl 464 `"464"', add
label define conspuma_lbl 465 `"465"', add
label define conspuma_lbl 466 `"466"', add
label define conspuma_lbl 467 `"467"', add
label define conspuma_lbl 468 `"468"', add
label define conspuma_lbl 469 `"469"', add
label define conspuma_lbl 470 `"470"', add
label define conspuma_lbl 471 `"471"', add
label define conspuma_lbl 472 `"472"', add
label define conspuma_lbl 473 `"473"', add
label define conspuma_lbl 474 `"474"', add
label define conspuma_lbl 475 `"475"', add
label define conspuma_lbl 476 `"476"', add
label define conspuma_lbl 477 `"477"', add
label define conspuma_lbl 478 `"478"', add
label define conspuma_lbl 479 `"479"', add
label define conspuma_lbl 480 `"480"', add
label define conspuma_lbl 481 `"481"', add
label define conspuma_lbl 482 `"482"', add
label define conspuma_lbl 483 `"483"', add
label define conspuma_lbl 484 `"484"', add
label define conspuma_lbl 485 `"485"', add
label define conspuma_lbl 486 `"486"', add
label define conspuma_lbl 487 `"487"', add
label define conspuma_lbl 488 `"488"', add
label define conspuma_lbl 489 `"489"', add
label define conspuma_lbl 490 `"490"', add
label define conspuma_lbl 491 `"491"', add
label define conspuma_lbl 492 `"492"', add
label define conspuma_lbl 493 `"493"', add
label define conspuma_lbl 494 `"494"', add
label define conspuma_lbl 495 `"495"', add
label define conspuma_lbl 496 `"496"', add
label define conspuma_lbl 497 `"497"', add
label define conspuma_lbl 498 `"498"', add
label define conspuma_lbl 499 `"499"', add
label define conspuma_lbl 500 `"500"', add
label define conspuma_lbl 501 `"501"', add
label define conspuma_lbl 502 `"502"', add
label define conspuma_lbl 503 `"503"', add
label define conspuma_lbl 504 `"504"', add
label define conspuma_lbl 505 `"505"', add
label define conspuma_lbl 506 `"506"', add
label define conspuma_lbl 507 `"507"', add
label define conspuma_lbl 508 `"508"', add
label define conspuma_lbl 509 `"509"', add
label define conspuma_lbl 510 `"510"', add
label define conspuma_lbl 511 `"511"', add
label define conspuma_lbl 512 `"512"', add
label define conspuma_lbl 513 `"513"', add
label define conspuma_lbl 514 `"514"', add
label define conspuma_lbl 515 `"515"', add
label define conspuma_lbl 516 `"516"', add
label define conspuma_lbl 517 `"517"', add
label define conspuma_lbl 518 `"518"', add
label define conspuma_lbl 519 `"519"', add
label define conspuma_lbl 520 `"520"', add
label define conspuma_lbl 521 `"521"', add
label define conspuma_lbl 522 `"522"', add
label define conspuma_lbl 523 `"523"', add
label define conspuma_lbl 524 `"524"', add
label define conspuma_lbl 525 `"525"', add
label define conspuma_lbl 526 `"526"', add
label define conspuma_lbl 527 `"527"', add
label define conspuma_lbl 528 `"528"', add
label define conspuma_lbl 529 `"529"', add
label define conspuma_lbl 530 `"530"', add
label define conspuma_lbl 531 `"531"', add
label define conspuma_lbl 532 `"532"', add
label define conspuma_lbl 533 `"533"', add
label define conspuma_lbl 534 `"534"', add
label define conspuma_lbl 535 `"535"', add
label define conspuma_lbl 536 `"536"', add
label define conspuma_lbl 537 `"537"', add
label define conspuma_lbl 538 `"538"', add
label define conspuma_lbl 539 `"539"', add
label define conspuma_lbl 540 `"540"', add
label define conspuma_lbl 541 `"541"', add
label define conspuma_lbl 542 `"542"', add
label define conspuma_lbl 543 `"543"', add
label values conspuma conspuma_lbl

label define appal_lbl 0 `"Not in Appalachia"'
label define appal_lbl 1 `"Northern Appalachia"', add
label define appal_lbl 2 `"Central Appalachia"', add
label define appal_lbl 3 `"Southern Appalachia"', add
label values appal appal_lbl

label define appald_lbl 00 `"Not in Appalachia"'
label define appald_lbl 10 `"Northern Appalachia"', add
label define appald_lbl 11 `"Northern Applachia"', add
label define appald_lbl 12 `"North Central Appalachia"', add
label define appald_lbl 20 `"Central Appalachia"', add
label define appald_lbl 30 `"Southern Appalachia"', add
label define appald_lbl 31 `"South Central Appalachia"', add
label define appald_lbl 32 `"Southern Appalachia"', add
label values appald appald_lbl

label define homeland_lbl 1 `"PUMA does not include a homeland area"'
label define homeland_lbl 2 `"PUMA includes a homeland area"', add
label values homeland homeland_lbl

label define cntry_lbl 630 `"Puerto Rico"'
label define cntry_lbl 840 `"United States"', add
label values cntry cntry_lbl

label define gq_lbl 0 `"Vacant unit"'
label define gq_lbl 1 `"Households under 1970 definition"', add
label define gq_lbl 2 `"Additional households under 1990 definition"', add
label define gq_lbl 3 `"Group quarters--Institutions"', add
label define gq_lbl 4 `"Other group quarters"', add
label define gq_lbl 5 `"Additional households under 2000 definition"', add
label define gq_lbl 6 `"Fragment"', add
label values gq gq_lbl

label define acrehous_lbl 0 `"N/A"'
label define acrehous_lbl 1 `"House on less than 10 acres"', add
label define acrehous_lbl 2 `"House on 10 acres or more"', add
label define acrehous_lbl 3 `"House on less than 3 cuerdas (1980-1990)"', add
label define acrehous_lbl 4 `"House on 3+ cuerdas (1980-1990)"', add
label define acrehous_lbl 5 `"House on less than 10 cuerdas (2000 and PRCS)"', add
label define acrehous_lbl 6 `"House on 10 or more cuerdas (2000 and PRCS)"', add
label values acrehous acrehous_lbl

label define rentgrs_lbl 0000 `"N/A"'
label define rentgrs_lbl 0010 `"$1-19"', add
label define rentgrs_lbl 0025 `"$20-29"', add
label define rentgrs_lbl 0035 `"$30-39"', add
label define rentgrs_lbl 0045 `"$40-49"', add
label define rentgrs_lbl 0055 `"$50-59"', add
label define rentgrs_lbl 0065 `"$60-69"', add
label define rentgrs_lbl 0075 `"$70-79"', add
label define rentgrs_lbl 0090 `"$80-99"', add
label define rentgrs_lbl 0110 `"$100-119"', add
label define rentgrs_lbl 0135 `"$120-149"', add
label define rentgrs_lbl 0175 `"$150-199"', add
label define rentgrs_lbl 0200 `"$200+"', add
label values rentgrs rentgrs_lbl

label define kitchen_lbl 0 `"N/A"'
label define kitchen_lbl 1 `"No"', add
label define kitchen_lbl 2 `"No, or shared use"', add
label define kitchen_lbl 3 `"Yes, shared use"', add
label define kitchen_lbl 4 `"Yes (shared or exclusive use)"', add
label define kitchen_lbl 5 `"Yes, exclusive use"', add
label values kitchen kitchen_lbl

label define rooms_lbl 00 `"N/A"'
label define rooms_lbl 01 `"1 room"', add
label define rooms_lbl 02 `"2"', add
label define rooms_lbl 03 `"3"', add
label define rooms_lbl 04 `"4"', add
label define rooms_lbl 05 `"5"', add
label define rooms_lbl 06 `"6"', add
label define rooms_lbl 07 `"7"', add
label define rooms_lbl 08 `"8"', add
label define rooms_lbl 09 `"9 (9+, 1960-2007)"', add
label define rooms_lbl 10 `"10"', add
label define rooms_lbl 11 `"11"', add
label define rooms_lbl 12 `"12"', add
label define rooms_lbl 13 `"13"', add
label define rooms_lbl 14 `"14"', add
label define rooms_lbl 15 `"15"', add
label define rooms_lbl 16 `"16"', add
label define rooms_lbl 17 `"17"', add
label define rooms_lbl 18 `"18"', add
label define rooms_lbl 19 `"19"', add
label define rooms_lbl 20 `"20"', add
label define rooms_lbl 21 `"21"', add
label define rooms_lbl 22 `"22"', add
label define rooms_lbl 23 `"23"', add
label define rooms_lbl 24 `"24"', add
label define rooms_lbl 25 `"25"', add
label define rooms_lbl 26 `"26"', add
label define rooms_lbl 27 `"27"', add
label values rooms rooms_lbl

label define plumbing_lbl 00 `"N/A"'
label define plumbing_lbl 10 `"Without complete plumbing"', add
label define plumbing_lbl 11 `"Lacking only hot water"', add
label define plumbing_lbl 12 `"Lacking other or all plumbing facilities"', add
label define plumbing_lbl 13 `"Has some facilities"', add
label define plumbing_lbl 14 `"Has no facilities"', add
label define plumbing_lbl 20 `"With complete plumbing"', add
label define plumbing_lbl 21 `"Used only by household"', add
label define plumbing_lbl 22 `"Shared with others"', add
label values plumbing plumbing_lbl

label define builtyr2_lbl 00 `"N/A"'
label define builtyr2_lbl 01 `"1939 or earlier"', add
label define builtyr2_lbl 02 `"1940-1949"', add
label define builtyr2_lbl 03 `"1950-1959"', add
label define builtyr2_lbl 04 `"1960-1969"', add
label define builtyr2_lbl 05 `"1970-1979"', add
label define builtyr2_lbl 06 `"1980-1989"', add
label define builtyr2_lbl 07 `"1990-1994 (1990-1999 in the 2005-2011 ACS and the PRCS)"', add
label define builtyr2_lbl 08 `"1995-1999 (1995-1998 in the 2000-2002 ACS)"', add
label define builtyr2_lbl 09 `"2000-2004 (1999-2002 in the 2000-2002 ACS)"', add
label define builtyr2_lbl 10 `"2005 (2005 or later in the 2005-2007 and 2006-2011 ACS/PRCS)"', add
label define builtyr2_lbl 11 `"2006"', add
label define builtyr2_lbl 12 `"2007"', add
label define builtyr2_lbl 13 `"2008"', add
label define builtyr2_lbl 14 `"2009"', add
label define builtyr2_lbl 15 `"2010"', add
label define builtyr2_lbl 16 `"2011"', add
label values builtyr2 builtyr2_lbl

label define unitsstr_lbl 00 `"N/A"'
label define unitsstr_lbl 01 `"Mobile home or trailer"', add
label define unitsstr_lbl 02 `"Boat, tent, van, other"', add
label define unitsstr_lbl 03 `"1-family house, detached"', add
label define unitsstr_lbl 04 `"1-family house, attached"', add
label define unitsstr_lbl 05 `"2-family building"', add
label define unitsstr_lbl 06 `"3-4 family building"', add
label define unitsstr_lbl 07 `"5-9 family building"', add
label define unitsstr_lbl 08 `"10-19 family building"', add
label define unitsstr_lbl 09 `"20-49 family building"', add
label define unitsstr_lbl 10 `"50+ family building"', add
label values unitsstr unitsstr_lbl

label define bedrooms_lbl 00 `"N/A"'
label define bedrooms_lbl 01 `"No bedrooms"', add
label define bedrooms_lbl 02 `"1"', add
label define bedrooms_lbl 03 `"2"', add
label define bedrooms_lbl 04 `"3"', add
label define bedrooms_lbl 05 `"4 (4+ in 1960)"', add
label define bedrooms_lbl 06 `"5+ (1970-2000, ACS, PRCS)"', add
label define bedrooms_lbl 07 `"6"', add
label define bedrooms_lbl 08 `"7"', add
label define bedrooms_lbl 09 `"8"', add
label define bedrooms_lbl 10 `"9"', add
label define bedrooms_lbl 11 `"10"', add
label define bedrooms_lbl 12 `"11"', add
label define bedrooms_lbl 13 `"12"', add
label define bedrooms_lbl 14 `"13"', add
label define bedrooms_lbl 15 `"14"', add
label define bedrooms_lbl 16 `"15"', add
label define bedrooms_lbl 17 `"16"', add
label define bedrooms_lbl 18 `"17"', add
label define bedrooms_lbl 19 `"18"', add
label define bedrooms_lbl 20 `"19"', add
label define bedrooms_lbl 21 `"20"', add
label define bedrooms_lbl 22 `"21"', add
label values bedrooms bedrooms_lbl

label define qacrehou_lbl 0 `"Not allocated or N/A"'
label define qacrehou_lbl 3 `"Logical edit"', add
label define qacrehou_lbl 4 `"Allocated, hot deck"', add
label define qacrehou_lbl 5 `"Cold deck allocation (select variables)"', add
label values qacrehou qacrehou_lbl

label define qbedroom_lbl 0 `"Not allocated"'
label define qbedroom_lbl 3 `"Allocated, direct"', add
label define qbedroom_lbl 4 `"Allocated, hot deck"', add
label define qbedroom_lbl 5 `"Cold deck allocation (select variables)"', add
label define qbedroom_lbl 9 `"Allocated, direct/indirect"', add
label values qbedroom qbedroom_lbl

label define qbuilty2_lbl 0 `"Not allocated"'
label define qbuilty2_lbl 3 `"Allocated, direct"', add
label define qbuilty2_lbl 4 `"Allocated"', add
label define qbuilty2_lbl 5 `"Allocated, indirect"', add
label define qbuilty2_lbl 6 `""Don't Know""', add
label define qbuilty2_lbl 9 `"Allocated, direct/indirect"', add
label values qbuilty2 qbuilty2_lbl

label define qkitchen_lbl 0 `"Not allocated"'
label define qkitchen_lbl 3 `"Allocated, direct"', add
label define qkitchen_lbl 4 `"Allocated"', add
label define qkitchen_lbl 5 `"Cold deck allocation (select variables)"', add
label define qkitchen_lbl 9 `"Allocated, direct/indirect"', add
label values qkitchen qkitchen_lbl

label define qrooms_lbl 0 `"Not allocated"'
label define qrooms_lbl 4 `"Allocated"', add
label values qrooms qrooms_lbl

label define qunitsst_lbl 0 `"Not allocated"'
label define qunitsst_lbl 3 `"Allocated, direct"', add
label define qunitsst_lbl 4 `"Allocated"', add
label define qunitsst_lbl 5 `"Cold deck allocation-select variables"', add
label define qunitsst_lbl 9 `"Allocated, direct/indirect"', add
label values qunitsst qunitsst_lbl

! gzip -f ACS2005raw.dat

save ACS2005raw, replace

log close

