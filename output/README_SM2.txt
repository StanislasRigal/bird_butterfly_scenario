###
Explanation of columns:
###

Species: species scientific name
Name: species English name
Habitat.indicator: is the species is listed in FaBI, FoBI, GBI or WBI
All.bird/butterfly.index:  is the species is listed in BiI or BuI
Habitat.indicator.used: is the species included in the FaBI, FoBI, GBI or WBI computed in the article
trend_past: additive trend (natural logarithm of the multiplicative trend) between 2000 and 2021 at the European scale

##### Covariates effect on local abundance

open_vs_forest: effect of open land category versus forest category on local abundance
urban_vs_forest: effect of urban category versus forest category on local abundance
spring_mean_temperature: effect of mean temperature between March and June on local abundance
spring_rainfall: effect of rainfall amount between March and June on local abundance
landscape_diversity: effect of Shannon landscape diversity (based on CORINE classes) on local abundance
primary_productivity: effect of drymatter production (i.e. local primary production) on local abundance

##### Effect of covariate trends on species trend

urbanisation_increase_effect_on_trend: change in yearly slope due to increase in urbanisation
treedensity_increase_effect_on_trend: change in yearly slope due to increase in tree density
non_intensively_managed_forest_effect_on_trend: change in yearly slope when amount of non intensively managed forests is high
intensively_managed_forest_effect_on_trend: change in yearly slope when amount of intensively managed forests is high
farmland_increase_effect_on_trend: change in yearly slope due to increase in farmland area
low_intensively_managed_farmland_effect_on_trend: change in yearly slope when amount of low intensively managed farmland is high
high_intensively_managed_farmland_effect_on_trend: change in yearly slope when amount of high intensively managed farmland is high
temperature_increase_effect_on_trend: change in yearly slope due to increase in spring temperature
temperature_variation_increase_effect_on_trend: change in yearly slope due to increase in variance of spring temperature
rainfall_increase_effect_on_trend: change in yearly slope due to increase in spring rainfall
landscape_diversity_increase_effect_on_trend: change in yearly slope due to increase in shannon landscape diversity
protected_area_effect_on_trend: change in yearly slope when amount of protected area is high

##### Standard deviation of the effect of covariate trends on species trend

Same column names with “_sd” suffix.

dev_exp: deviance explain by the model
