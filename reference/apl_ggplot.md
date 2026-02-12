# Plot Association Plot with ggplot

Uses ggplot to plot an Association Plot

## Usage

``` r
apl_ggplot(
  rows,
  rows_group = NULL,
  cols,
  cols_group = NULL,
  rows_scored = NULL,
  rows_color = "#0066FF",
  rows_high_color = "#FF0000",
  cols_color = "#601A4A",
  cols_high_color = "#EE442F",
  score_color = "rainbow",
  row_labs = FALSE,
  col_labs = FALSE,
  show_score = FALSE,
  show_cols = FALSE,
  show_rows = TRUE
)
```

## Arguments

- rows:

  Row APL-coordinates

- rows_group:

  Row AP-coordinates to highlight

- cols:

  Column AP-coordinates

- cols_group:

  Column AP-coordinates for the group to be highlighted.

- rows_scored:

  Row AP-coordinates of rows above a score cutoff.

- rows_color:

  Color for rows

- rows_high_color:

  Color for rows to be highlighted.

- cols_color:

  Column points color.

- cols_high_color:

  Color for column points to be highlighted..

- score_color:

  Color scheme for row points with a score.

- row_labs:

  Logical. Whether labels for rows indicated by rows_idx should be
  labeled with text. Default TRUE.

- col_labs:

  Logical. Whether labels for columns indicated by cols_idx shouls be
  labeled with text. Default FALSE.

- show_score:

  Logical. Whether the S-alpha score should be shown in the plot.

- show_cols:

  Logical. Whether column points should be plotted.

- show_rows:

  Logical. Whether row points should be plotted.

## Value

ggplot Association Plot
