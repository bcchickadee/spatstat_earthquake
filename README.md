# Point Pattern Analyis of Earthquakes in Indonesia

By Sumin Yun [(LinkedIn)](https://www.linkedin.com/in/bcchickadee/)

## Introduction

This GitHub repository contains the source code for "Point Pattern Analysis of Earthquakes in Indonesia". This project was intially created as the final report of the course "Statistical Analysis for Spatial Data", a graduate-level course in the Department of Statistics in Seoul National University.

## Data Usage

This project utilizes mainly 3 data sources:

1. The ["Earthquakes in Indonesia" dataset from Kaggle](https://www.kaggle.com/datasets/kekavigi/earthquakes-in-indonesia). Contains data of earthquakes in Indonesia from 2008 to 2023. Also contains numerous covariates regarding earthquake information, such as the magnitude, magnitude type, depth (in km), azimuth gap, and a short description of the location.
2. The [Global Earthquake Model (GEM) Active Faults Database](https://github.com/GEMScienceTools/gem-global-active-faults). Contains spatial distribution of faults worldwide, including Indonesia. The dataset provides vector data of crustal faults and subduction trenches.
3. [Slab2](https://www.usgs.gov/data/slab2-a-comprehensive-subduction-zone-geometry-model), a 3D model of subduction zone geometry produced by USGS. Slab2 contains NetCDF grids containing the depth of the slab surface, dip, and strike at regular intervals. We use multiple mosaic grids for several subduction regions in Indonesia: Sumatra-Java (`sum`), Sulawesi (`sul`), Halmahera (`hal`), Philippines (`phi`), and New Guinea (`png`).

Also, `R` libraries such as `rnaturalearth` were utilized for map visualization.

## Content

The repository contains a single `R` file (`script.r`) that contains the code used in this project. The script contains the code needed to output the figures and tables in the report, as well as the code that constructs the models and necessary spatial objects.