# Antigen Selection Justification

This document provides detailed justification for the selection of specific antigens for the mixture model GWAS analysis, based on three key criteria:

1. **Distribution Characteristics:** Evidence of bimodality with overlap from diagnostic plots
2. **Genetic Signal:** Genome-wide significant associations from Butler-Laporte et al. (2020)
3. **Power Gain Potential:** Estimated inflation factor (1/λ²) indicating meaningful power loss with hard cutoffs

## Selected Antigens

### CMV (Cytomegalovirus)

#### cmv_pp150
- **Threshold:** 100
- **Butler-Laporte Evidence:** pp28 MFI shows rs12698418 (chr7, P=4.4e-08) near EN2 gene
- **Rationale:** CMV is highly prevalent (~50-80% seropositivity) with clear bimodal distributions
- **Expected Power Gain:** High - CMV typically shows excellent component separation

#### cmv_pp28  
- **Threshold:** 200
- **Butler-Laporte Evidence:** rs12698418 (chr7, P=4.4e-08) near EN2 gene
- **Rationale:** Strong genetic signal outside MHC region, indicating non-HLA genetic control
- **Expected Power Gain:** Moderate - good separation but higher threshold may reduce misclassification

### HSV (Herpes Simplex Virus)

#### hsv1
- **Threshold:** 150
- **Butler-Laporte Evidence:** No genome-wide significant hits in Table 2
- **Rationale:** High prevalence (~60-90% seropositivity), clear bimodal distribution, important for herpesvirus family representation
- **Expected Power Gain:** High - typically excellent component separation

#### hsv2
- **Threshold:** 150  
- **Butler-Laporte Evidence:** 
  - Seropositivity: rs538162817 (chr3, P=2.9e-08) near GRK7
  - mgG-2 MFI: rs144232229 (chr3, P=1.9e-09) near TAMM41
- **Rationale:** Multiple genetic associations outside MHC, strong evidence for genetic control
- **Expected Power Gain:** High - clear genetic signal with good distribution separation

### EBV (Epstein-Barr Virus)

#### ebv_vca
- **Threshold:** 250
- **Butler-Laporte Evidence:**
  - rs9379862 (chr6, P=1.1e-08) near BTN3A2
  - rs9271536 (chr6, P=7.3e-21) in HLA-DQA1 (MHC)
- **Rationale:** Multiple genetic associations including strong MHC signal, EBV is ubiquitous
- **Expected Power Gain:** High - excellent bimodal distribution

#### ebv_ebna1
- **Threshold:** 250
- **Butler-Laporte Evidence:**
  - rs67886110 (chr3, P=2.1e-08) near MED12L
  - rs6927022 (chr6, P=1.9e-76) in HLA-DQA1 (MHC)
- **Rationale:** Strongest genetic signal among EBV antigens, multiple associations
- **Expected Power Gain:** High - clear genetic control demonstrated

#### ebv_zebra
- **Threshold:** 100
- **Butler-Laporte Evidence:**
  - rs34034915 (chr6, P=3.3e-09) near RP1-97D16.1
  - Chr6:32597087 (P=3.0e-75) in HLA-DQA1 (MHC)
- **Rationale:** Very strong MHC association, important for EBV lytic cycle
- **Expected Power Gain:** Moderate - good separation but complex MHC genetics

#### ebv_ead
- **Threshold:** 100
- **Butler-Laporte Evidence:**
  - rs2316515 (chr6, P=4.2e-08) near IRF4
  - rs2395192 (chr6, P=3.1e-16) in HLA-DRB9 (MHC)
- **Rationale:** IRF4 association is biologically plausible for EBV, additional MHC signal
- **Expected Power Gain:** Moderate - good component separation

### H. pylori

#### hp_omp
- **Threshold:** 170
- **Butler-Laporte Evidence:** rs3104361 (chr6, P=6.5e-10) in HLA-DQB1 (MHC)
- **Rationale:** Strong MHC association, H. pylori is important bacterial pathogen
- **Expected Power Gain:** Moderate - bacterial antigens often show good separation

#### hp_urea
- **Threshold:** 130
- **Butler-Laporte Evidence:** rs71569678 (chr6, P=3.0e-08) near RP11-439H9.1
- **Rationale:** Additional H. pylori antigen, different genetic architecture than OMP
- **Expected Power Gain:** Moderate - good for bacterial pathogen representation

### T. gondii

#### toxo_sag1
- **Threshold:** 160
- **Butler-Laporte Evidence:**
  - rs541989586 (chr3, P=3.9e-08) near SOX2-OT
  - rs148929820 (chr5, P=8.5e-10) near RP11-510I6.1
  - rs11881343 (chr19, P=3.9e-08) near NOTCH3
- **Rationale:** Multiple genetic associations outside MHC, important protozoal pathogen
- **Expected Power Gain:** High - clear genetic signal with good distribution

### Polyomaviruses

#### bkv_vp1
- **Threshold:** 250
- **Butler-Laporte Evidence:** rs492602 (chr19, P=4.3e-09) in FUT2
- **Rationale:** Strongest FUT2 association, important for polyomavirus family
- **Expected Power Gain:** High - clear genetic signal with FUT2

#### jcv_vp1
- **Threshold:** 250
- **Butler-Laporte Evidence:** rs2432132 (chr19, P=8.8e-15) in FUT2
- **Rationale:** Additional FUT2 association, different from BKV
- **Expected Power Gain:** High - strong genetic signal

#### mcv_vp1
- **Threshold:** 250
- **Butler-Laporte Evidence:** rs55792153 (chr5, P=3.6e-10) near TMEM173
- **Rationale:** TMEM173 association, important for polyomavirus family
- **Expected Power Gain:** Moderate - good genetic signal

### Additional Herpesviruses

#### vzv
- **Threshold:** 100
- **Butler-Laporte Evidence:** Multiple MHC hits including rs1766 (chr6, P=1.1e-11) in HLA-DQB1
- **Rationale:** Strong MHC associations, important herpesvirus
- **Expected Power Gain:** Moderate - good distribution but mainly MHC signal

#### hhv6_ie1a
- **Threshold:** 100
- **Butler-Laporte Evidence:** rs13079586 (chr3, P=4.6e-08) near ITGA9
- **Rationale:** Clear genetic signal outside MHC, good distribution
- **Expected Power Gain:** High - non-MHC genetic control

#### hhv6_ie1b
- **Threshold:** 100
- **Butler-Laporte Evidence:** rs28752523 (chr6, P=7.9e-09) in HLA-DQA1
- **Rationale:** MHC association, additional HHV-6 antigen
- **Expected Power Gain:** Moderate - MHC genetics

#### hhv7_u14
- **Threshold:** 100
- **Butler-Laporte Evidence:** rs139299944 (chr6, P=4.1e-12) in HLA-DQA1
- **Rationale:** Strong MHC association, high prevalence (~95%)
- **Expected Power Gain:** Moderate - strong MHC signal

### Chlamydia trachomatis

#### ct_pgp3
- **Threshold:** 200
- **Butler-Laporte Evidence:** Chr7:66874490 (P=4.2e-08) near AC006480.1
- **Rationale:** Bacterial pathogen representation, clear genetic signal
- **Expected Power Gain:** Moderate - good for bacterial pathogens

### Additional H. pylori Antigens

#### hp_caga
- **Threshold:** 400
- **Butler-Laporte Evidence:** No specific hits in Table 2
- **Rationale:** Important virulence factor, high threshold may reduce misclassification
- **Expected Power Gain:** Low - high threshold reduces power gain potential

#### hp_vaca
- **Threshold:** 100
- **Butler-Laporte Evidence:** No specific hits in Table 2
- **Rationale:** Additional H. pylori antigen, lower threshold
- **Expected Power Gain:** Moderate - good for bacterial pathogen representation

#### toxo_p22
- **Threshold:** 100
- **Butler-Laporte Evidence:** No specific hits in Table 2
- **Rationale:** Additional T. gondii antigen, lower threshold
- **Expected Power Gain:** Moderate - good for protozoal pathogen representation

## Excluded Antigens

### Chlamydia trachomatis (Additional antigens)
- **Reason:** Limited genetic signal for most antigens beyond pGP3
- **Butler-Laporte Evidence:** Only pGP3 shows clear genetic signal
- **Note:** momp A, momp D, tarp-D F1, tarp-D F2, PorB included for completeness but low priority

### H. pylori (Additional antigens)
- **Reason:** Limited genetic signal for some antigens
- **Butler-Laporte Evidence:** Only OMP and UreA show clear genetic signals
- **Note:** CagA, VacA, GroEL, Catalase included for completeness but lower priority

### T. gondii p22
- **Reason:** Limited genetic signal compared to sag1
- **Butler-Laporte Evidence:** No specific hits in Table 2
- **Note:** Included for completeness but sag1 is the primary T. gondii antigen

## Summary

The selected antigens provide:
1. **Diverse Pathogen Coverage:** Herpesviruses, bacteria, and protozoa
2. **Strong Genetic Evidence:** Multiple genome-wide significant associations
3. **Good Distribution Properties:** Clear bimodality with overlap requiring mixture models
4. **High Prevalence:** Most antigens have >20% seropositivity ensuring adequate sample sizes

This selection maximizes the potential for discovering novel genetic associations while ensuring robust statistical power through the mixture model approach. 