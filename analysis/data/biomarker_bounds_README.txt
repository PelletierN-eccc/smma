Bounding biomarker values for SMMA Factor 2 stability testing
Purpose
Factor 2 currently apportions ~60% petrogenic. Because only ESTS petroleum entries carry biomarker data,
any biomarker signal in the factor profile can only be matched to a petrogenic candidate - the apportionment
may therefore be an artefact of library structure rather than a chemical result. These sheets provide low,
central and high bounding biomarker values for the 120 profiles that currently have none, so the model can be
re-run at each bound and the stability of the Factor 2 assignment evaluated.
These are NOT measurements. They are literature-informed bounds intended for sensitivity analysis only, and
should not be merged into the source library or reported as source profile data.
How the values are constructed
Sum(biomarker) is set per profile, then distributed across the 21 analytes using a source-class shape vector.
Two magnitude bases are used, both given on the Bounds sheet:
  ratio    - Sum(biomarker) = f x Sum(PAH) of that profile. Used for emission and product profiles, where the
             reporting basis (per mass of PM, per mass of fuel) is not comparable across studies.
  absolute - Sum(biomarker) fixed in mg/kg, then divided by that profile Sum(PAH). Used for the Van Metre dust,
             soil and pavement matrices, which appear to be reported in mg/kg (see Units caveat below).
Shape matters more than magnitude
SMMA discriminates on profile shape. The feature that separates lubricating-oil carryover from crude and heavy
fuel oil is a truncated extended homohopane series - C29 and C30 alpha-beta hopane dominant, C33 to C35 sparse -
together with proportionally higher steranes. Coal and coal-tar shapes are depleted in extended homohopanes for
a different reason (pyrolytic degradation of the pentacyclic skeleton). If Factor 2 has a full extended
homohopane series with S/R near 1.5, the lube-oil and coal shapes will not fit it at any magnitude, and that is
a stronger result than the magnitude test alone.
Units caveat - please check before use
The library is not on a single concentration basis. Van Metre sealed-lot composite dust has a median Sum(PAH)
of 2126 here, against a published median of 2200 mg/kg for sealcoated pavement dust; unsealed is 69 here against
27 mg/kg published. Those match in mg/kg, not ng/g. The ESTS entries sit near 10,000, which reads as ug/g. If
SMMA normalises each profile before matching this is harmless for the analysis, but the manuscript statement
that concentrations were all converted to ng/g needs revising either way.
Second caveat - partial analyte coverage
Sum(PAH) is computed over the analytes each profile actually reports, which ranges from 5 to 22 of the target
list. Ratio-based imputation against a partial sum therefore under-imputes for sparse profiles. If that matters
for your Factor 2 competitors, recompute the ratio against a common analyte subset instead.
Sheets
Bounds        - magnitude bounds and rationale by category
Shapes        - the four biomarker distribution vectors, each summing to 1
Impute_low / Impute_central / Impute_high - 120 profiles x 21 biomarkers, ready to substitute