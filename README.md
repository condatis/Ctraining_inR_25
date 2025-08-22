---
title: "Condatis in R 2025 (Ctraining_inR_25)"
author: "Jenny Hodgson and Claudia Gutierrez-Arellano" 

output:
  html_document: 
  
editor_options: 
  markdown: 
    wrap: 72
---

*Training first delivered in September 2025, by Jenny Hodgson and
Claudia Gutierrez-Arellano*

This github repository contains most of the materials you will need to
go through the training examples for yourself.

If you're unfamiliar with github, don't worry! Just click on the green
button that says "\<\> Code" to reveal a short menu, and then choose
"Download ZIP". You can unzip the download anywhere that's convenient on
your own computer.

## Session summary - Condatis connectivity analysis to plan resilient habitat networks

Condatis is decision support software to identify the best locations for
habitat restoration to increase connectivity across landscapes. It can
help planners to find multiple-win solutions for land management, so
that biodiversity can recover and natural ecosystems can co-exist with
human activities. Condatis’s metric of long-distance connectivity is
particularly relevant under climate change, ensuring that networks of
habitat facilitate range shifts. Condatis:

1.  Highlights pathways across a landscape that allow both dispersal and
    reproduction of species;

2.  Pinpoints bottlenecks in the habitat network, where there are
    restricted opportunities for colonisation, and where restoration
    would be most impactful;

3.  Ranks the feasible sites for habitat restoration, to efficiently
    enhance the existing habitat network.

Condatis has been developed since 2015 and has been used by conservation
practitioners around the world, via its web application
([www.condatis.org.uk](www.condatis.org.uk)). Through collaboration with
users, we have designed it to be fast and user-friendly, and to work
with the kind of data that organisations typically have. This workshop
will provide a brief overview of the software and a series of graded
exercises, showcasing the different functions of Condatis.

## Background requirements

There will be examples using the Condatis web app as well as the R
environment. Attendees will need some familiarity with GIS, e.g.
reading, writing and querying raster data. You may use your own favoured
GIS software or spatial package to view our input and output data,
either R, ArcGIS Pro or QGIS. It may be advisable to update your GIS
package to the latest version. We will not be offering in-depth GIS help
during the session. R is essential for one of our exercises (for the
analysis and not just for its GIS capability). Please install R version
4.3 or later before you start, and install the packages terra, sf,
sfheaders, tidyverse and viridis. Before using the Condatis web app, you
need to create a user account for yourself at
<https://webapp.condatis.org.uk>.

## Programme

In the [Workshop
programme](https://condatis.github.io/Ctraining_inR_25/), you will find
a schedule and a short description of the activities planned

## Copyright

Please note that we are making this code and documentation available
under a GPL license; for more details see the LICENSE file. Copyright
Jenny Hodgson and Claudia Gutierrez-Arellano 2025 unless otherwise
stated (e.g. datasets). Please cite this archive as Jenny Hodgson and
Claudia Gutierrez-Arellano (2025) Condatis connectivity analysis to plan
resilient habitat networks - Condatis training in R. Available at:
<https://github.com/condatis/Ctraining_inR_25>.
