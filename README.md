# LEO-CubeSat-Constellation-Design-for-Remote-Monitoring

## Abstract

Nepal's glacial lakes, floods and landslides pose constant threats, yet ground sensors are costly and often absent in the most vulnerable regions. LEO CubeSats offer frequent passes at low cost, but a constellation with a good average revisit can still leave long stretches of the day uncovered because passes cluster. This project uses Bayesian Optimization to choose the orbits of a CubeSat constellation for monitoring Nepal's glacial lakes, rather than fixing a Walker or Sun-synchronous pattern in advance. The optimisation targets the shape of the daily revisit schedule, not only its average.

## Process

1. **Simulation:** a single MATLAB script on the Satellite Communications Toolbox models a 24-hour scenario with the Nepal border, glacial lake targets and a ground station. Each candidate orbit is scored by full SGP4 propagation and coverage simulation.
2. **Baseline orbit:** Latin Hypercube Sampling seeds a Gaussian Process surrogate, and Expected Improvement proposes each next orbit, balancing coverage, revisit, ground-station contact, latency and lake coverage.
3. **Incremental growth:** satellites are added one at a time. An inner search places a satellite in an existing plane, and an outer search adds a new plane. A hard veto rejects any step that worsens worst-lake revisit.
4. **Temporal objective:** Wasserstein distance, largest-gap and gap-uniformity terms favour satellites that fly at times the constellation leaves uncovered.
5. **Verification:** independent geometric checks and per-lake analysis confirm coverage and revisit.
6. **Remote monitoring:** a proposed SAR payload gives day-and-night, cloud-independent imaging, with Sentinel-1 data used to demonstrate change detection.

## Result

A grown multi-plane constellation clearly improves on the optimised single satellite in coverage, average revisit and downlink latency, while keeping every target glacial lake covered. The remaining weakness is a large overnight gap in passes, which is geometric rather than a failure of the search: satellites in the same plane can only shift passes within that plane's own time window.

## Conclusion

Optimisation-driven design produces an effective small constellation for glacial lake monitoring in Nepal. The daily gap is structural, so the next steps are more orbital planes with RAAN spaced around the day, multi-day and end-of-life analysis, validation against real TLE data, and completing the SAR monitoring chain with Google Earth Engine.

## Requirements

- MATLAB with the Satellite Communications Toolbox
- Nepal border shapefile

