All data of observations on and around the Juneau Icefield divide, where we are simulating firn processes using the Meyer-Hewitt continuum model for meltwater flow through compacting snow (2017).

# juneauIcefield_massBalance_v1.0

*From JIF_GlacierWide_README*

Direct field measurements of point glaciological data are combined with weather and geodetic data to estimate the seasonal and annual mass balance at each glacier in both a conventional and reference surface format (Cogley and others, 2011). The basic analysis framework (O'Neel, 2019, in prep; McNeil et. al, 2019) is the same at each glacier to enable cross-comparison between output time series. However, in this data release for Taku and Lemon Creek glaciers temperature lapse rates are optimized using on-icefield weather data. This changes the degree day factor in the melt model, giving small post-geodetic calibration differences on the order of 2-3 cm. Details are described in McNeil (2019). Vocabulary used follows Cogley and others (2011) Glossary of Glacier Mass Balance.

**Reference:** O'Neel and Others (2019)

**CRS:** UTM Zone 8N (EPSG 26908) ref to WGS84

**Units:** m.w.e. 

**Data of interest:** (from JuneauICefield_UTM8N.csv) 
	- MG6 (some yearly data 1962-2019 in Input_Taku_Glaciological_Data) and TSQG1 (none in Input_Taku_Glaciological_Data); other close datapoints include MG2 and LLG1
	- Input_Taku_Glaciological_Data (MG6)
	- Input_Taku_Daily_Weather (from airport, elev 5m, [58.3566, -134.5640])

**JuneauIcefieldDivide QGIS:** JuneauIcefield_UTM8N_massbalance

# juneauIcfield_weather_v1.0

*From README_WeatherData_JIF.pdf*

**Reference:** Baker, E. H., McGee, S., Campbell, S. W., Pierce, J. L., McNeil, C. J., 2019, Weather Station Data on the Juneau Icefield (ver. 1.0, June 2019): U.S. Geological Survey data release, https://doi.org/10.5066/P90RCN51.

**Units:** Celsius (temperature), incremental - mm (precipitation, mulitply by site specific factor), cumulative - m, W m^2 (radiation)

**CRS:** WGS 84

**Data of Interest:** C26 is closest; Precip and Temp; Lvl1 and Lvl2 data are both fine for this analysis

