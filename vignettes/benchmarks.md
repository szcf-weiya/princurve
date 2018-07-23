Benchmarking princurve
================
Robrecht Cannoodt

<!-- github markdown using
rmarkdown::render("vignettes/benchmarks.Rmd", output_format = "github_document")
-->
princurve 2.1 contains major optimisations if the `approx_points` parameter is used. This is showcased on a toy example, where the number of points was varied between 10<sup>2</sup> and 10<sup>6</sup>.

``` r
# warning: running this will take ± 2 hours
# princurvelegacy can be installed using devtools::install_github("dynverse/princurve@legacy")
num_points <- round(10^seq(log10(100), log10(100000), length.out = 25))
lambda <- rnorm(max(num_points), 0, .2)
x <- cbind(lambda, lambda^2) + rnorm(length(lambda) * 2, 0, .02)

benchmarks <- map_df(
  rev(num_points),
  function(np) {
    xsel <- x[seq_len(np),]
    microbenchmark::microbenchmark(
      "princurve1.1" = princurvelegacy::principal.curve(xsel),
      "princurve2.1" = princurve::principal_curve(xsel, approx_points = 100),
      times = 10L,
      unit = "ms"
    ) %>%
      summary() %>%
      as_data_frame() %>%
      mutate(num_points = np)
  }
)
```

We can see princurve 2.1 scales linearly w.r.t. the number of rows in the dataset, whereas princurve 1.1 scales quadratically. This is due to the addition of the approximation step added in between the smoothing and the projection steps (explained in more detail in the [algorithm](algorithm.md) vignette).

``` r
ggplot(benchmarks, aes(num_points, median / 1000)) +
  geom_point() +
  geom_line() +
  facet_wrap(~expr, ncol = 1, scales = "free") +
  theme_bw() +
  labs(x = "Number of rows in dataset", y = "Time (s)") +
  scale_colour_brewer(palette = "Set1")
```

![](benchmarks_files/figure-markdown_github/compare-1.png)