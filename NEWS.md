# switchbox 2.1.0.9001

## Bug fixes

- Fixed an issue where `swap_train_1tsp()` could store a gene pair in the
  wrong order when trained with an unsigned scoring function (i.e.
  `score_fn = swap_calculate_basic_tsp_scores`). This could produce a
  classifier whose predictions did not match its reported score, and in some
  cases a much less accurate classifier than intended. Classifiers trained
  with the default (signed) scoring function were not affected.

- Fixed the ROC curve returned by `swap_ktsp_result()` and
  `swap_train_test_results()` (and therefore plotted by `swap_plot_roc()`),
  which was mirrored across the diagonal. The curve — and the AUC shown in
  its legend — previously appeared worse than the classifier actually
  performed. `swap_prediction_stats()`'s reported `auc` value was always
  correct; only the plotted curve was affected.

- Fixed `swap_calculate_signed_tsp_scores()` with `handle_ties = TRUE`: when
  called without `restricted_pairs`, the returned scores were exactly double
  their correct magnitude. This did not change which feature pairs were
  selected or how classifiers built from them predicted — only the numeric
  score values reported to the user (e.g. in `classifier$score` or
  `scores$score`) were affected.

## Documentation

- Clarified the "Getting Started" vignette's explanation of what the
  `gene1`/`gene2` order in `classifier$TSPs` means for prediction.
