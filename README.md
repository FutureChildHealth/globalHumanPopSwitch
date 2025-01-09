# Global human population ended self-facilitation in the 1950s
<a href="https://zenodo.org/doi/10.5281/zenodo.10951676"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.10951676.svg" alt="DOI"></a>
<img align="right" src="www/pop.png" alt="population icon" width="180" style="margin-top: 20px">

The global human population switched from depensation to compensation in the 1950s

<br>
Prof <a href="https://globalecologyflinders.com/people/#DIRECTOR">Corey J. A. Bradshaw</a> <br>
<a href="http://globalecologyflinders.com" target="_blank">Global Ecology</a> | <em><a href="https://globalecologyflinders.com/partuyarta-ngadluku-wardli-kuu/" target="_blank">Partuyarta Ngadluku Wardli Kuu</a></em>, <a href="http://flinders.edu.au" target="_blank">Flinders University</a>, Adelaide, Australia <br>
February 2024 <br>
<a href=mailto:corey.bradshaw@flinders.edu.au>e-mail</a> <br>
<br>

Accompanies preprint:<br>
<br>
Bradshaw, CJA, MA Judge, DT Blumstein, PR Ehrlich, ANZ Dasgupta, M Wackernagel, LJZ Weeda, PN Le Souëf. <a href="https://doi.org/10.2139/ssrn.4788003">Global human population ended self-facilitation in the 1950s</a>. <em>SSRN</em> doi:10.2139/ssrn.4788003

## Abstract
Applied to human populations, the ecological concept of carrying capacity is necessarily complicated because human beings are the “ultimate ecosystem engineers” who deliberately and successfully moderate their environment for their benefit. For at least that last few hundred years, human ingenuity, access to massive stocks of fossil fuels, and technological development have driven facilitation whereby increasing human abundance promoted higher population growth rates. However, this positive relationship broke down during the 1950s, and by 1962, the global human population entered a negative phase where the growth rate consistently declined as population increased. The onset of the negative phase was approximately contemporaneous with a global biocapacity deficit that began in 1970 and has increased ever since. The onset of the negative phase varies regionally, with the lowest-income regions entering this phase later than higher-income regions. A Ricker logistic model fitted to the negative phase predicts that the global population could reach 11.6 to 12.3 billion people between 2065 and 2074. The same model fitted to the facilitation phase predicts a maximum sustainable population of 2.5 billion people that Earth’s carrying capacity might be able to maintain over a longer time period. This is in line with economics-based estimates of equitable wealth distribution (3.3 billion people) and an ecological footprint of 0.5 planets (2.35 billion people). The negative phase also correlates strongly with the trend in global temperature anomaly. The Earth cannot sustain the future human population, or even today’s, without a major overhaul of socio-cultural practices for using land, water, energy, biodiversity, and other resources. Sustainable development therefore continues to face the challenge of meeting the needs and aspirations of current societies while simultaneously ensuring that future generations can meet their own.

## Scripts
- <code>humanpoptransition.R</code> (main code)
- <code>new_lmer_AIC_tables3.R</code> (source functions)
- <code>r.squared.R</code> (source functions)

## Data
- <em>agestructure1950_2021.csv</em>: <a href="https://data.un.org/Data.aspx?d=POP&f=tableCode%3a1">United Nations Population Division age structure data from 1950-2021</a>
- <em>AGO.csv</em>: <a href="https://data.un.org/Data.aspx?d=POP&f=tableCode%3a1">United Nations Population Division data for Angola</a>
- <em>BDI.csv</em>: <a href="https://data.un.org/Data.aspx?d=POP&f=tableCode%3a1">United Nations Population Division data for Burundi</a>
- <em>BFA.csv</em>: <a href="https://data.un.org/Data.aspx?d=POP&f=tableCode%3a1">United Nations Population Division data for Burkina Faso</a>
- <em>china1950-2021.csv</em>: <a href="https://data.un.org/Data.aspx?d=POP&f=tableCode%3a1">United Nations Population Division data for China</a>
- <em>COD.csv</em>: <a href="https://data.un.org/Data.aspx?d=POP&f=tableCode%3a1">United Nations Population Division data for Republic of Congo</a>
- <em>consump.csv</em>: <a href="[https://data.un.org/Data.aspx?d=POP&f=tableCode%3a1](https://www.energyinst.org/statistical-review/resources-and-data-downloads)">global consumption data</a>
- <em>GMB.csv</em>: <a href="https://data.un.org/Data.aspx?d=POP&f=tableCode%3a1">United Nations Population Division data for Gambia</a>
- <em>MLI.csv</em>: <a href="https://data.un.org/Data.aspx?d=POP&f=tableCode%3a1">United Nations Population Division data for Mali</a>
- <em>Nchild-01.csv</em>: number of children aged 0-1 years from 1950-2021 (United Nations Population Division)
- <em>NER.csv</em>: <a href="https://data.un.org/Data.aspx?d=POP&f=tableCode%3a1">United Nations Population Division data for Niger</a>
- <em>NGA.csv</em>: <a href="https://data.un.org/Data.aspx?d=POP&f=tableCode%3a1">United Nations Population Division data for Nigeria</a>
- <em>Npre1950.csv</em>: <a href="http://www.census.gov/data/tables/time-series/demo/international-programs/historical-est-worldpop.html">historical estimates of global population size</a>
- <em>popregions.csv</em>: population data by major global region (United Nations Population Division)
- <em>popXtempanom.csv</em>: world population size relative to global temperature anomaly (<a href="http://www.metoffice.gov.uk/hadobs/hadcrut5/data/HadCRUT.5.0.2.0/download.html">HadCRUT.05.0.3.0</a>)
- <em>TCD.csv</em>: <a href="https://data.un.org/Data.aspx?d=POP&f=tableCode%3a1">United Nations Population Division data for Chad</a>
- <em>UGA.csv</em>: <a href="https://data.un.org/Data.aspx?d=POP&f=tableCode%3a1">United Nations Population Division data for Uganda</a>
- <em>UNpop.csv</em>: <a href="https://data.un.org/Data.aspx?d=POP&f=tableCode%3a1">United Nations Population Division data</a>
- <em>worldpophist.csv</em>: <a href="https://data.un.org/Data.aspx?d=POP&f=tableCode%3a1">historical United Nations Population Division data</a>

## R libraries
<code>plotrix</code>, <code>boot</code>, <code>tmvnsim</code>, <code>wCorr</code>, <code>truncnorm</code>, <code>orcutt</code>, <code>lmtest</code>, <code>performance</code>, <code>sjPlot</code>, <code>dismo</code>, <code>gbm</code>
<br>
<br>

<p><a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Horizontal_RGB_Master.png" alt="Flinders University" width="150" style="margin-top: 20px"></a> &nbsp; <a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="GEL" width="85" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://www.uwa.edu.au/"><img align="bottom-left" src="www/uwa2.png" alt="UWA" width="100" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://www.telethonkids.org.au"><img align="bottom-left" src="www/tkilogo.png" alt="TKI" width="90" style="margin-top: 20px"></a> &nbsp; &nbsp; &nbsp; <a href="https://www.ucla.edu/"><img align="bottom-left" src="www/uclalogo.png" alt="UCLA logo" width="90" style="margin-top: 20px"></a> <a href="https://www.stanford.edu/"><img align="bottom-left" src="www/stanfordlogo.webp" alt="Stanford logo" width="90" style="margin-top: 20px"></a> &nbsp; <a href="https://www.footprintnetwork.org/"><img align="bottom-left" src="www/GFNlogo.png" alt="Global Footprint Network logo" width="90" style="margin-top: 20px"></a>
