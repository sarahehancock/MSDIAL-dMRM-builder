# MS-DIAL to Agilent dMRM builder

Repository: [https://github.com/sarahehancock/MSDIAL-dMRM-builder](https://github.com/sarahehancock/MSDIAL-dMRM-builder)

Github page: [https://sarahehancock.github.io/MSDIAL-dMRM-builder/](https://sarahehancock.github.io/MSDIAL-dMRM-builder/)

A Shiny app that turns MS-DIAL lipidomics alignment results into a positive-mode
dynamic MRM (dMRM) transition list for an Agilent triple quadrupole, using
negative-mode data to identify acyl chains.

## Run it

```r
install.packages(c("shiny", "DT"))
shiny::runApp("app")
```

Needs R 4.1 or later. No other packages.

## Inputs

| File | Required | Used for |
|---|---|---|
| Positive MS-DIAL alignment result export (`.txt`, tab-delimited) | yes | Targets, precursors, retention times, product ions, abundance |
| Negative MS-DIAL alignment result export | no | Acyl-chain identities only |

The exports must include the `MS/MS spectrum` column. Sample intensity columns
and the `File type` header row are used to work out abundance.

## How transitions are chosen

| Decision | Source |
|---|---|
| Which features are targeted | `Metabolite name`: every annotated positive row (unknowns skipped) |
| Precursor m/z and RT | `Reference m/z` and `Average Rt(min)` |
| Class and sum composition | `Metabolite name` only |
| Molecular species | The name after `\|` in positive and the matched negative feature, plus chains mined from the `MS/MS spectrum` column of both files |
| Product ion | The most intense expected fragment present in the positive MS2 spectrum |
| Product m/z in the export | Calculated from the class rule, not the measured peak |

1. **Quantifier**: the most intense class-specific product (head group or
   neutral loss). TG and CL have none, so a chain fragment is used.
2. **One transition is kept** unless another target shares the Q1 window and
   co-elutes (a discriminating second product is added, or the feature is
   flagged `RT separation only`), or more than one isomeric species was found
   (each isomer gets a transition on a chain the others do not carry).
3. **PC and ether PC** use m/z 184 only by default; species come from negative
   mode.
4. **Capacity**: where more MRMs are concurrent than
   `cycle time / (minimum dwell + inter-scan delay)` allows, the lowest-priority
   TG transitions are dropped. Priority is the TG's mean abundance in the
   samples in use, scaled for isomer-specific ions by relative fragment
   intensity. The `Dropped TG` tab gives the reason for each.

Positive-mode fragment rules are listed in the app's `Rules` tab.

## Outputs

- **Agilent dMRM list (.csv)**: Compound Group, Compound Name, ISTD?, Precursor
  Ion, MS1 Res, Product Ion, MS2 Res, Primary, Trigger, Ret Time (min), Delta
  Ret Time, Fragmentor, Collision Energy, Cell Accelerator Voltage, Polarity.
- **Annotated transitions (.csv)**: the same rows with fragment annotation,
  evidence, abundance, worst-case concurrency and dwell.
- **Dropped TG list (.csv)**.

## Check before acquiring

- **Collision energies are starting values**, one per class, not optimised.
  Edit them in the `Collision energies` tab.
- **Fragmentor (380 V), inter-scan delay (3.5 ms) and minimum dwell (3 ms) are
  defaults**, not instrument-verified. Set them for your model.
- **The CSV header has not been tested against MassHunter.** Compare it with a
  method exported from your instrument. Compound names contain
  `: _ / ; ( ) [ ] - .` from lipid nomenclature.
- **Retention times come from the MS-DIAL run.** Use RT scale and offset if the
  QQQ gradient differs.
- **Peak width (12.5 s), window (1.0 min) and cycle time (1000 ms) defaults**
  were measured on one 25-minute reversed-phase method.
- Transitions flagged `predicted` were not seen in the positive MS2 spectrum.

## Hosting

Shiny needs R to run, so the app cannot be served as plain files. Two options:

- **GitHub Pages via shinylive** (`.github/workflows/deploy-pages.yml`): the
  workflow compiles the app to a static site where R runs in the visitor's
  browser. Uploaded files never leave their computer. Enable under
  Settings > Pages > Source: GitHub Actions. First load is slow (tens of MB).
- **A Shiny server** such as shinyapps.io or Posit Connect Cloud:
  deploy the `app` folder. Uploaded files are processed on that server.

## License

GNU General Public License v3.0, see `LICENSE`. Copyright 2026 Sarah Hancock.
