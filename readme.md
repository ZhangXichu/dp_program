## 🌦️ dp_program — Spatio-Temporal Interpolation in R

This project was created for my [Master's thesis](https://is.muni.cz/th/rtyw9/) at Masaryk University.

It focuses on interpolation of temperature and precipitation over Czechia using:

- 🌐 Generalized Additive Models (GAM) with splines

- 🗺️ Spatio-temporal kriging

- 📊 Exploratory spatial data analysis and variogram fitting

The core implementation is in dp_impl.R, with supporting theory and simulations in dp_theory.R. memory_management.R is used to list all R objects in memory, sorted by size, so large ones can be removed.