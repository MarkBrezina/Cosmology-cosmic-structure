# Run readout

Stage: development

Reserved test opened: False

Lower RMS log-speed error is better. Every galaxy has equal aggregate weight.

```
           model subset  galaxies  points  rmse_galaxy_equal_km_s  rms_log10_speed_dex  standardized_mse_galaxy_equal  invalid_points  eligible_for_comparison
             rar    all       117    2212               20.092914             0.076218                       7.369391               0                     True
      shared_nfw    all       117    2212               29.828982             0.131306                      16.312418               0                     True
 radius_gradient    all       117    2212               33.716775             0.135619                      22.145293               0                     True
density_gradient    all       117    2212               33.917573             0.186282                      20.553949               1                    False
 screened_lambda    all       117    2212               41.645556             0.201261                      32.242274               0                     True
         baryons    all       117    2212               41.649377             0.201276                      32.247939               0                     True
```

See paired_bootstrap.csv for paired galaxy-resampling comparisons. Negative candidate-minus-reference values favor density_gradient.

Inspect fold_fits.json (or fixed_development_calibration.json) for optimizer failures and boundary hits.

A fit can reject or support this particular phenomenological force law; it cannot establish a unified dark-energy explanation.

The density proxies are not independent measurements of local gas+stellar volume density or external cosmic void density.

The RAR benchmark was developed from SPARC itself, so its family has historical dataset exposure. Its scale parameter here is fitted only to each training fold.

The shared NFW scaling is a restricted two-parameter population rule plus mass calibration, not a full LambdaCDM model comparison.

Void observations are not fitted by this galaxy run. See void_test.py and data/VOID_INPUT.md.
