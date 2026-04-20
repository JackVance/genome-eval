"""Trait configurations for the generic PRS runner.

Each entry provides:
  - pgs_id: PGS Catalog ID for the weight file
  - label: plain-language trait name
  - topic: ledger topic key
  - r_squared: variance explained in the source paper's European validation
  - anchors (optional): population mean + SD per sex, in the trait's units
  - units: the trait's natural units ("cm", "years of schooling", "IQ points", "SD", ...)
  - claim_template + notes_template: format strings for the ledger
  - citation: short citation string
  - replication_count + inference_confidence: tiering inputs
  - self_report_key (optional): the profile's self_reported_phenotypes key (e.g. "iq")
      that corresponds to this trait. If set, run_prs.py ALWAYS emits a
      self-report cross-check block in the finding notes (MATCH / MISMATCH /
      "no self-report on file"). This is a Rule-10/11 guardrail: you cannot
      accidentally ship a PRS finding without the cross-check.
  - self_report_value_field (optional): name of the numeric field to read inside
      the self-report dict. Common values: "value_cm" for height, "value_approx"
      for IQ, "value_years" for educational attainment.
  - self_report_unit_conversion (optional): factor to convert self-report to the
      trait's ledger unit if they differ (rare — defaults to identity).

Anthropometric / cognitive anchors are rough European figures; update per
population if a different subject ancestry is used.
"""
from __future__ import annotations

TRAITS = {
    # --- Height --- (existing, using Lello-style smaller panel)
    'height_xie2020': {
        'pgs_id': 'PGS000297',
        'label': 'Height',
        'topic': 'prs_height_pgs000297',
        'r_squared': 0.20,
        'units': 'cm',
        'anchors': {
            'male': {'mean': 177.0, 'sd': 7.0},
            'female': {'mean': 164.0, 'sd': 6.5},
        },
        'self_report_key': 'height',
        'self_report_value_field': 'value_cm',
        'citation': 'Xie T et al. Circ Genom Precis Med (2020). DOI: 10.1161/circgen.119.002775',
        'replication_count': 10,
        'inference_confidence': 'moderate',
        'claim_template': (
            '{label}: PRS z-score {z} (percentile ~{percentile}) -> predicted {predicted} '
            '(95% CI {ci_low}-{ci_high} cm); population SD is {pheno_sd} cm per 1 SD; r^2={r2}.'
        ),
        'notes_template': (
            'Height PRS from {pgs_id} ({label}). Panel: {panel} variants; {contributing} contributing '
            '(coverage {coverage}). r^2={r2} in European validation; remaining '
            'variance is unexplained by this score (other genetic + environmental factors).'
        ),
    },

    # --- Height, Yengo 2022 (larger, more predictive) ---
    'height_yengo2022': {
        'pgs_id': 'PGS002804',
        'label': 'Height (Yengo 2022 European, 1.1M variants, SBayesC)',
        'topic': 'prs_height_pgs002804',
        'r_squared': 0.40,
        'units': 'cm',
        'anchors': {
            'male': {'mean': 177.0, 'sd': 7.0},
            'female': {'mean': 164.0, 'sd': 6.5},
        },
        'self_report_key': 'height',
        'self_report_value_field': 'value_cm',
        'citation': 'Yengo L et al. Nature (2022). A saturated map of common genetic variants associated with human height. N=5.4M.',
        'replication_count': 15,
        'inference_confidence': 'high',
        'claim_template': (
            '{label} (polygenic prediction): {predicted} '
            '(95% CI {ci_low}-{ci_high} cm); percentile ~{percentile}; '
            'PRS z-score {z}; r^2={r2} — state-of-the-art common-variant coverage.'
        ),
        'notes_template': (
            'Height PRS from {pgs_id} (Yengo 2022 saturated map, N=5.4M). Panel: {panel} variants; '
            '{contributing} contributing (coverage {coverage}). r^2={r2} in European validation. '
            'This is the best-predicting common-variant height PRS currently published. '
            'Residual variance reflects rare variants, non-additive effects, and environmental factors.'
        ),
    },

    # --- Educational attainment (Privé LDpred2 on UK Biobank years of education) ---
    'educational_attainment': {
        'pgs_id': 'PGS002231',
        'label': 'Educational attainment (years of schooling proxy)',
        'topic': 'prs_educational_attainment_pgs002231',
        'r_squared': 0.12,
        'units': 'years of schooling',
        'anchors': {
            'any': {'mean': 14.0, 'sd': 3.0},  # US adult mean ~14 yrs; SD ~3
        },
        'self_report_key': 'educational_attainment',
        'self_report_value_field': 'value_years',
        'citation': 'Privé F et al. Am J Hum Genet (2022). Portability of 245 polygenic scores when derived from the UK Biobank and applied to 9 ancestry groups. Qualifications / years of education, ~950k variants, LDpred2.',
        'replication_count': 8,
        'inference_confidence': 'moderate',
        'claim_template': (
            'Educational attainment: PRS z-score {z} (percentile ~{percentile}) -> predicted {predicted} '
            '(95% CI {ci_low}-{ci_high} years); population SD is {pheno_sd} years per 1 SD; r^2={r2}. '
            'Individual prediction is weak at this r^2. Within-family effects are ~half of between-family '
            '(gene-environment correlation). Not a judgment of ability.'
        ),
        'notes_template': (
            'Educational-attainment PRS from {pgs_id} (Okbay 2022 EA4, N=3M). Panel: {panel} variants; '
            '{contributing} contributing (coverage {coverage}). r^2={r2} in European validation — '
            'population-level association is genuine but individual prediction is weak. '
            'Within-family PRS effects are consistently found to be about half the size of between-family '
            'effects, meaning much of the cross-family signal is gene-environment correlation '
            '(including assortative mating and parenting) rather than direct genetic effect. '
            'This score should not be interpreted as a measure of individual ability or predetermined outcome.'
        ),
    },

    # --- Coronary artery disease (Patel 2023 GPS_Mult) ---
    # Binary outcome — no continuous anchor. Percentile + OR/SD is the
    # meaningful read; empirical 1000G EUR calibration produces the z-score.
    'coronary_artery_disease': {
        'pgs_id': 'PGS003725',
        'label': 'Coronary artery disease (Patel 2023 GPS_Mult)',
        'topic': 'prs_cad_pgs003725',
        'r_squared': 0.10,        # liability-scale approximation from OR/SD ≈ 2.14
        'units': '(binary)',       # no continuous trait, predicted is None
        'anchors': {},             # intentionally empty
        'self_report_key': None,   # no self-report on file; would be yes/no diagnosis
        'citation': 'Patel AP et al. Nat Med (2023). A multi-ancestry polygenic risk score improves risk prediction for CAD. N=1.4M training.',
        'replication_count': 12,
        'inference_confidence': 'moderate',
        'claim_template': (
            '{label} (polygenic risk): PRS z-score {z}; percentile ~{percentile}. '
            'Source-paper OR per SD ≈ 2.14 (HR/SD ≈ 1.75). r^2 (liability-scale est.) ≈ {r2}. '
            'CAD is a common-outcome binary trait; interpret as percentile, not a predicted value.'
        ),
        'notes_template': (
            'Coronary-artery-disease PRS from {pgs_id} (Patel 2023 GPS_Mult, N=1.4M; 269k cases). '
            'Panel: {panel} variants; {contributing} contributing (coverage {coverage}). '
            'Benchmark: outperforms Inouye 2018 metaGRS and Khera 2018 in independent validation. '
            'Multi-ancestry design; primary validation cohorts are EUR. r^2 is liability-scale '
            'estimate (~{r2}) from reported OR/SD ≈ 2.14. Not a diagnosis — risk is probabilistic '
            'and only one factor among age, sex, lipids, BP, smoking, diabetes, family history.'
        ),
    },

    # --- Type 2 diabetes (Ge 2022 PRS-CSx on DIAMANTE) ---
    'type_2_diabetes': {
        'pgs_id': 'PGS002308',
        'label': 'Type 2 diabetes (Ge 2022, PRS-CSx on DIAMANTE)',
        'topic': 'prs_t2d_pgs002308',
        'r_squared': 0.092,
        'units': '(binary)',
        'anchors': {},
        'self_report_key': None,
        'citation': 'Ge T et al. Am J Hum Genet (2022). Development and validation of a trans-ancestry PRS for T2D. DIAMANTE N=1.1M.',
        'replication_count': 10,
        'inference_confidence': 'moderate',
        'claim_template': (
            '{label} (polygenic risk): PRS z-score {z}; percentile ~{percentile}. '
            'Source-paper AUROC = 0.793; OR/SD ≈ 1.96; R^2 = {r2}. Interpret as percentile.'
        ),
        'notes_template': (
            'T2D PRS from {pgs_id} (Ge 2022 PRS-CSx on DIAMANTE, N=1.1M; 81.7% EUR + 16.1% EAS + 2.2% AFR). '
            'Panel: {panel} variants; {contributing} contributing (coverage {coverage}). '
            'EUR-validation AUROC = 0.793; OR/SD = 1.96; R^2 = {r2}. '
            'Risk modulated by BMI, diet, physical activity — PRS is one factor.'
        ),
    },

    # --- LDL cholesterol (Graham 2021 multi-ancestry GLGC) ---
    'ldl_cholesterol': {
        'pgs_id': 'PGS000889',
        'label': 'LDL cholesterol (Graham 2021 GLGC multi-ancestry)',
        'topic': 'prs_ldl_pgs000889',
        'r_squared': 0.158,
        'units': 'mg/dL',
        'anchors': {
            # Adult US reference: NHANES median LDL-C ~112 mg/dL, SD ~35. Untreated.
            'any': {'mean': 112.0, 'sd': 35.0},
        },
        'self_report_key': None,
        'citation': 'Graham SE et al. Nature (2021). The power of genetic diversity in genome-wide association studies of lipids. GLGC N=1.65M.',
        'replication_count': 12,
        'inference_confidence': 'moderate',
        'claim_template': (
            '{label} (polygenic prediction): PRS z-score {z} (percentile ~{percentile}) '
            '→ predicted {predicted} (95% CI {ci_low}-{ci_high} mg/dL); population SD = '
            '{pheno_sd} mg/dL per 1 SD; r^2 = {r2}. Untreated, fasting adult baseline. '
            'Current LDL depends heavily on diet + statins + age.'
        ),
        'notes_template': (
            'LDL-C PRS from {pgs_id} (Graham 2021, GLGC multi-ancestry, N=1.65M). '
            'Panel: {panel} variants; {contributing} contributing (coverage {coverage}). '
            'R^2 = {r2} in MVP EUR validation (n=68,381). Predicted value is untreated, '
            'fasting adult baseline — real LDL depends on diet, statin use, age, weight.'
        ),
    },

    # --- Systolic blood pressure (Shetty 2023 All-of-Us PRS-CS) ---
    'systolic_blood_pressure': {
        'pgs_id': 'PGS003971',
        'label': 'Systolic blood pressure (Shetty 2023 AoU PRS-CS)',
        'topic': 'prs_sbp_pgs003971',
        'r_squared': 0.18,
        'units': 'mmHg',
        'anchors': {
            # US adult mean SBP ~122 mmHg, SD ~15. Varies with age; use mid-adult.
            'any': {'mean': 122.0, 'sd': 15.0},
        },
        'self_report_key': None,
        'citation': 'Shetty PB et al. (2023). PRS-CS blood-pressure scores from All-of-Us. N=127k training.',
        'replication_count': 8,
        'inference_confidence': 'moderate',
        'claim_template': (
            '{label} (polygenic prediction): PRS z-score {z} (percentile ~{percentile}) '
            '→ predicted {predicted} (95% CI {ci_low}-{ci_high} mmHg); population SD = '
            '{pheno_sd} mmHg per 1 SD; r^2 = {r2}. Genetic component only — measured BP '
            'depends heavily on age, BMI, sodium intake, stress, measurement technique.'
        ),
        'notes_template': (
            'SBP PRS from {pgs_id} (Shetty 2023 All-of-Us PRS-CS, HapMap3 variants). '
            'Panel: {panel} variants; {contributing} contributing (coverage {coverage}). '
            'Multi-ancestry training; EUR validation R^2 = {r2}. Predicted value is a '
            'genetic-propensity baseline; actual SBP drifts substantially with age + '
            'modifiable factors.'
        ),
    },

    # --- IQ / cognitive ability (Privé LDpred2 on UK Biobank intelligence phenotype) ---
    'cognitive_ability': {
        'pgs_id': 'PGS002135',
        'label': 'General cognitive ability (g-like)',
        'topic': 'prs_cognitive_ability_pgs002135',
        'r_squared': 0.07,
        'units': 'IQ points',
        'anchors': {
            'any': {'mean': 100.0, 'sd': 15.0},
        },
        'self_report_key': 'iq',
        'self_report_value_field': 'value_approx',
        'citation': 'Privé F et al. Am J Hum Genet (2022). Portability of 245 polygenic scores when derived from the UK Biobank and applied to 9 ancestry groups. Intelligence phenotype, ~903k variants, LDpred2.',
        'replication_count': 5,
        'inference_confidence': 'low',
        'claim_template': (
            'Cognitive ability (g-like): PRS z-score {z} (percentile ~{percentile}) -> predicted {predicted} '
            '(95% CI {ci_low}-{ci_high} IQ points); population SD is {pheno_sd} IQ points per 1 SD; r^2={r2}. '
            'Individual prediction is weak at this r^2. Not a measure of intelligence.'
        ),
        'notes_template': (
            'Cognitive ability PRS from {pgs_id} (Savage 2018, N=269,867). Panel: {panel} variants; '
            '{contributing} contributing (coverage {coverage}). r^2={r2} in European validation — '
            'substantially weaker than educational-attainment PRS and far weaker than the heritability '
            'estimates from twin/family studies (50-80%). The gap reflects missing heritability: rare '
            'variants, non-additive effects, phenotype measurement noise, and the modest sample size '
            'relative to the polygenic architecture. This score should NOT be interpreted as an '
            'intelligence measurement; it is a weak statistical signal useful at population level, '
            'noisy at individual level. Tail predictions (very high or very low) are particularly '
            'unreliable due to regression-to-the-mean in partial-coverage PRS.'
        ),
    },
}
