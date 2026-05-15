# COVID-19 Survey Data from Statistics Austria (November 2020)

Survey sample of \\n = 2287\\ participants tested for COVID-19 in
November 2020 by Statistics Austria, used as the case study in Guerrier,
Kuzmics & Victoria-Feser (2024).

## Usage

``` r
covid19_austria
```

## Format

A `data.frame` with 2287 rows and 3 variables:

- Y:

  Binary; 1 if participant \\i\\ tested positive in the survey sample, 0
  otherwise.

- Z:

  Binary; 1 if participant \\i\\ was declared positive by the official
  procedure, 0 otherwise.

- weights:

  Sampling weights.

## Source

Statistics Austria. 2020. "Prävalenz von SARS-CoV-2-Infektionen liegt
bei 3,1%."

## References

Guerrier, S., Kuzmics, C., Victoria-Feser, M.-P. (2024). Assessing
COVID-19 Prevalence in Austria with Infection Surveys and Case Count
Data as Auxiliary Information. *Journal of the American Statistical
Association*, 119(547), 1722-1735.
[doi:10.1080/01621459.2024.2313790](https://doi.org/10.1080/01621459.2024.2313790)
