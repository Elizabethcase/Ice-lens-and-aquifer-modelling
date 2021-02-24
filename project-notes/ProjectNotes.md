# Modelling the Juneau Icefield with M&H's Firn Aquifer Model

The Juneau Icefield is a set of temperate glaciers. Much circumstantial evidence points to the existence of a firn aquifer, and potentially, a configuration of stacked aquifers that form each season and are underlain by a thick ice lens.

## Main working questions

1. What aquifer structures/configurations does the existing data support?
2. Can the Meyer/Hewitt Model produce a stacked firn aquifer?
	- Secondary question: can it produce big ice slabs?
3. What kinds of aquifer configurations have been found?

## To do 
### Feb 22-26
- [ ] create Simplified (no compaction, constant water input, etc) version of RunCode
- [x] create gaus, exp, uniform, & ice lens initial phi
### Feb 1-5
- [[20210202-Tue]]: run systematically varying Q + A 
### Jan 25-29
### Jan 18-22
- [ ] transfer and clean up MH paper onto Obsidian
### Jan 4-8
- [x] create [[firn aquifer file structure | file structure]] for project
- [x] create initial literature review for project
- [ ] create to do for MH code
- [x] transfer notes on MH paper from Notion to Obsidian
- [x] develop [[firn aquifer timeline | timeline]] for project

## Notes on data
### Juneau Icefield
#### Weather
- [description of data](/Users/elizabeth/Documents/home_research/projects/JIRP Firn Aquifer/data/readme_data.md)
- as of 2/2/2021
> *From README_WeatherData_JIF.pdf*
> **Reference:** Baker, E. H., McGee, S., Campbell, S. W., Pierce, J. L., McNeil, C. J., 2019, Weather Station Data on the Juneau Icefield (ver. 1.0, June 2019): U.S. Geological Survey data release, https://doi.org/10.5066/P90RCN51.
> **Units:** Celsius (temperature), incremental - mm (precipitation, mulitply by site specific factor), cumulative - m, W m^2 (radiation)
> **CRS:** WGS 84
> **Data of Interest:** C26 is closest; Precip and Temp; Lvl1 and Lvl2 data are both fine for this analysis


#### Mass balance
- [description of data](/Users/elizabeth/Documents/home_research/projects/JIRP Firn Aquifer/data/readme_data.md)
- as of 2/2/2021

> *From JIF_GlacierWide_README*
> Direct field measurements of point glaciological data are combined with weather and geodetic data to estimate the seasonal and annual mass balance at each glacier in both a conventional and reference surface format (Cogley and others, 2011). The basic analysis framework (O'Neel, 2019, in prep; McNeil et. al, 2019) is the same at each glacier to enable cross-comparison between output time series. However, in this data release for Taku and Lemon Creek glaciers temperature lapse rates are optimized using on-icefield weather data. This changes the degree day factor in the melt model, giving small post-geodetic calibration differences on the order of 2-3 cm. Details are described in McNeil (2019). Vocabulary used follows Cogley and others (2011) Glossary of Glacier Mass Balance.
> **Reference:** O'Neel and Others (2019)
> **CRS:** UTM Zone 8N (EPSG 26908) ref to WGS84
> **Units:** m.w.e. 
> **Data of interest:** (from JuneauICefield_UTM8N.csv)
> - MG6 (some yearly data 1962-2019 in Input_Taku_Glaciological_Data) and TSQG1 (none in Input_Taku_Glaciological_Data); other close datapoints include MG2 and LLG1
> - Input_Taku_Glaciological_Data (MG6)	
> - Input_Taku_Daily_Weather (from airport, elev 5m, [58.3566, -134.5640])

> **Related files**
> - **JuneauIcefieldDivide QGIS:** JuneauIcefield_UTM8N_massbalance
## Notes on code
- MH code description in [[@meyerContinuumModelMeltwater2017#^c629cb]]

## Literature
1. [Meyer CR and Hewitt IJ (2017)](zotero://open-pdf/library/items/XM76EABX) A continuum model for meltwater flow through compacting snow. _The Cryosphere_ **11**(6), 2799–2813 (doi:[10.5194/tc-11-2799-2017](https://doi.org/10.5194/tc-11-2799-2017)) [[@meyerContinuumModelMeltwater2017]]
2. [Fountain AG and Walder JS (1998)](zotero://open-pdf/library/items/DXJE9DM5) Water flow through temperate glaciers. _Reviews of Geophysics_ **36**(3), 299–328 (doi:[https://doi.org/10.1029/97RG03579](https://doi.org/10.1029/97RG03579))
<!-- <iframe width=900 height = 700 src = "https://www.connectedpapers.com/main/dced78217b94e98f2b621df8eafbfeeaf3582546/Water-flow-through-temperate-glaciers/graph"></iframe> -->
3. Angelen JH van, Lenaerts JTM, Broeke MR van den, Fettweis X and Meijgaard E van (2013) Rapid loss of firn pore space accelerates 21st century Greenland mass loss. _Geophysical Research Letters_ **40**(10), 2109–2113 (doi:[10.1002/grl.50490](https://doi.org/10.1002/grl.50490))
4. Ettema J, Broeke MR van den, Meijgaard E van, Berg WJ van de, Bamber JL, Box JE and Bales RC (2009) Higher surface mass balance of the Greenland ice sheet revealed by high-resolution climate modeling. _Geophysical Research Letters_ **36**(12) (doi:[10.1029/2009GL038110](https://doi.org/10.1029/2009GL038110))
5. Harper J, Humphrey N, Pfeffer WT, Brown J and Fettweis X (2012) Greenland ice-sheet contribution to sea-level rise buffered by meltwater storage in firn. _Nature_ **491**(7423), 240–243 (doi:[10.1038/nature11566](https://doi.org/10.1038/nature11566))
6. Jansson P, Hock R and Schneider T (2003) The concept of glacier storage: a review. _Journal of Hydrology_ **282**(1), 116–129 (doi:[10.1016/S0022-1694(03)00258-0](https://doi.org/10.1016/S0022-1694(03)00258-0))
7. Kawashima K and Yamada T (1997) Experimental studies on the transformation from firn to ice in the wet-snow zone of temperate glaciers. _Annals of Glaciology_ **24**, 181–185 (doi:[10.3189/S0260305500012143](https://doi.org/10.3189/S0260305500012143))
8. Machguth H, MacFerrin M, As D van, Box JE, Charalampidis C, Colgan W, Fausto RS, Meijer HAJ, Mosley-Thompson E and Wal RSW van de (2016) Greenland meltwater storage in firn limited by near-surface ice formation. _Nature Clim Change_ **6**(4), 390–393 (doi:[10.1038/nclimate2899](https://doi.org/10.1038/nclimate2899))
9. McNeil C, O’Neel S, Loso M, Pelto M, Sass L, Baker EH and Campbell S (2020) Explaining mass balance and retreat dichotomies at Taku and Lemon Creek Glaciers, Alaska. _Journal of Glaciology_ **66**(258), 530–542 (doi:[10.1017/jog.2020.22](https://doi.org/10.1017/jog.2020.22))
10. Miller O, Solomon DK, Miège C, Koenig L, Forster R, Schmerr N, Ligtenberg SRM and Montgomery L (2018) Direct Evidence of Meltwater Flow Within a Firn Aquifer in Southeast Greenland. _Geophysical Research Letters_ **45**(1), 207–215 (doi:[10.1002/2017GL075707](https://doi.org/10.1002/2017GL075707))
11. Ochwat NE, Marshall SJ, Moorman BJ, Criscitiello AS and Copland L (2020) Meltwater Storage in the firn of Kaskawulsh Glacier, Yukon Territory, Canada. _The Cryosphere Discussions_, 1–21 (doi:[https://doi.org/10.5194/tc-2020-119](https://doi.org/10.5194/tc-2020-119))
12. Peña S de la, Howat IM, Nienow PW, Broeke MR van den, Mosley-Thompson E, Price SF, Mair D, Noël B and Sole AJ (2015) Changes in the firn structure of the western Greenland Ice Sheet caused by recent warming. _The Cryosphere_ **9**(3), 1203–1211 (doi:[https://doi.org/10.5194/tc-9-1203-2015](https://doi.org/10.5194/tc-9-1203-2015))
13. Schneider T (2001) Hydrological processes in firn on Storglaciären, Sweden. [http://urn.kb.se/resolve?urn=urn:nbn:se:su:diva-144224](http://urn.kb.se/resolve?urn=urn:nbn:se:su:diva-144224)
14. Vandecrux B, Mottram R, Langen PL, Fausto RS, Olesen M, Stevens CM, Verjans V, Leeson A, Ligtenberg S, Kuipers Munneke P, Marchenko S, van Pelt W, Meyer CR, Simonsen SB, Heilig A, Samimi S, Marshall S, Machguth H, MacFerrin M, Niwano M, Miller O, Voss CI and Box JE (2020) The firn meltwater Retention Model Intercomparison Project (RetMIP): evaluation of nine firn models at four weather station sites on the Greenland ice sheet. _The Cryosphere_ **14**(11), 3785–3810 (doi:[https://doi.org/10.5194/tc-14-3785-2020](https://doi.org/10.5194/tc-14-3785-2020))
15. Williamson AG (2014) The Hydrological System of Storglaciären, Sweden: Integrating Modelling with Observations. Thesis, Scott Polar Research Institute, University of Cambridge. [https://www.repository.cam.ac.uk/handle/1810/264249](https://www.repository.cam.ac.uk/handle/1810/264249)