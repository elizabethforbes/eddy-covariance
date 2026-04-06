---
editor_options: 
  markdown: 
    wrap: 72
---

# Crop Staging Logic Reference Tables

## Cover Crop (Winter Rye) Phenology Thresholds

### Fall/Winter Planted Rye (September - November)

| Stage | Entry Conditions | Temperature Thresholds | NDVI Thresholds | Duration Logic | Exit Conditions |
|------------|------------|------------|------------|------------|------------|
| **Pre-planting** | `days_since_planting < -1` | N/A | N/A | Until planting | New planting date |
| **Germination/Emergence** | `days_since_planting ≤ 10` | N/A | N/A | Fixed: 10 days | Day counter |
| **Fall Tillering** | Post-emergence AND<br>- Soil ≥15°C OR Air_10 ≥15°C<br>- NDVI \>0.30 | **Soil:** ≥15°C (7-day smooth)<br>**Air:** ≥15°C (10-day smooth) | **NDVI:** \>0.30<br>(smoothed) | Active growth | Soil \<15°C in Oct+ |
| **Fall Hardening** | Oct+ AND<br>- Soil 5-15°C OR Air_10 \<10°C<br>- NDVI \>0.20 | **Soil:** 5-15°C (7-day smooth)<br>**Air:** \<10°C (10-day smooth) | **NDVI:** \>0.20 | **Sustained:** 3 consecutive days | Soil \<5°C sustained |
| **Dormant** | **Immediate:** Soil \<2°C OR Air_7 \<0°C<br>**Sustained:** Soil \<5°C for 10 days<br>**Forced:** Dec 1 if still active | **Immediate:** Soil \<2°C OR Air_7 \<0°C<br>**Sustained:** Soil \<5°C (7-day smooth)<br>**Force:** Dec 1 regardless | **NDVI:** N/A | **Sustained:** 10 consecutive days \<5°C<br>**Warning if:** Soil \>5°C on Dec 1 | Spring conditions OR May+ |
| **Spring Growth** | (Month \>3 OR Mar 20+) AND<br>- Soil \>5°C OR Air_10 \>10°C<br>- NDVI slope \>0.001<br>- NDVI \>0.25 | **Soil:** \>5°C (7-day smooth)<br>**Air:** \>10°C (10-day smooth) | **NDVI:** \>0.25<br>**NDVI slope:** \>0.001 (7-day) | **Sustained:** 10 consecutive days<br>**Calendar gate:** March 20+ | Soil \>10°C + rapid NDVI |
| **Rapid Growth** | From dormant/spring AND<br>- Soil \>10°C OR Air_10 \>15°C<br>- NDVI slope \>0.003<br>- NDVI \>0.35<br>- Month ≥3 | **Soil:** \>10°C (7-day smooth)<br>**Air:** \>15°C (10-day smooth) | **NDVI:** \>0.35<br>**NDVI slope:** \>0.003 (7-day) | Rapid green-up | Termination |
| **Terminated** | First tillage ≥60 days after planting | N/A | N/A | Permanent | End of record |

------------------------------------------------------------------------

### Warm-Season Planted Rye (March - August, e.g., July planting)

| Stage | Entry Conditions | Temperature Thresholds | NDVI Thresholds | Duration Logic | Exit Conditions |
|------------|------------|------------|------------|------------|------------|
| **Germination/Emergence** | `days_since_planting ≤ 10` | N/A | N/A | Fixed: 10 days | Day counter |
| **Seedling** | `days_since_planting 11-30` | N/A | N/A | Fixed: 20 days | Day counter |
| **Vegetative** | `days_since_planting 31-90` | N/A | N/A | Fixed: 60 days | Day counter |
| **Mature** | `days_since_planting >90` AND<br>Not meeting temp thresholds | Holding phase | N/A | Until fall hardening conditions | Temp-driven transitions |
| **Fall Hardening** | Oct+ AND<br>- Soil 5-15°C OR Air_10 \<10°C<br>- NDVI \>0.20 | **Soil:** 5-15°C (7-day smooth)<br>**Air:** \<10°C (10-day smooth) | **NDVI:** \>0.20 | **Sustained:** 3 consecutive days | Soil \<5°C |
| **Dormant** | **Immediate:** Soil \<2°C OR Air_7 \<0°C<br>**Sustained:** Soil \<5°C for 10 days<br>**Forced:** Dec 1 | **Immediate:** Soil \<2°C OR Air_7 \<0°C<br>**Sustained:** Soil \<5°C (7-day smooth) | N/A | **Sustained:** 10 consecutive days | Spring conditions |
| **Spring Growth** | Same as fall/winter planted | Same as fall/winter planted | Same as fall/winter planted | Same as fall/winter planted | Same as fall/winter planted |
| **Rapid Growth** | Same as fall/winter planted | Same as fall/winter planted | Same as fall/winter planted | Same as fall/winter planted | Termination |

**Important Note:** Warm-season planted rye does NOT have flowering or
seed production stages because it lacks the vernalization requirement
(30-60 days cold) needed to trigger reproductive development.

------------------------------------------------------------------------

## Alfalfa Hay Phenology Thresholds

### Temperature-Based Dormancy Override (October - April)

| Stage | Entry Conditions | Temperature Thresholds | NDVI Thresholds | Duration Logic | Exit Conditions |
|------------|------------|------------|------------|------------|------------|
| **Fall Hardening** | Oct-Nov AND<br>- Air_7 4.4-10°C<br>- Not yet dormant | **Air:** 4.4-10°C (7-day smooth)<br>**Months:** Oct-Nov only | N/A | **Sustained:** 3 consecutive days | Air \<0°C OR soil \<4°C OR NDVI \<0.20 |
| **Dormant** | **Direct:** NDVI \<0.20 OR Air_7 \<0°C OR Soil \<4°C<br>**From hardening:** Natural progression | **Air:** \<0°C (7-day smooth)<br>**Soil:** \<4°C (7-day smooth) | **NDVI:** \<0.20 (smoothed) | **Minimum:** 30 days before spring exit allowed | Spring conditions after Feb 1 |
| **Spring Growth** | **Calendar:** After Feb 1<br>**Duration:** 30+ days dormant<br>**Temp:** Air_10 \>5°C OR Soil \>5°C<br>**NDVI:** \>0.20 | **Air:** \>5°C (10-day smooth)<br>**Soil:** \>5°C (7-day smooth) | **NDVI:** \>0.20 | Transition phase | Air_10 \>10°C + NDVI slope \>0.003 |
| **Rapid Growth** | From dormant/spring/hardening AND<br>- Air_10 \>10°C OR Soil \>5°C<br>- NDVI slope \>0.003 | **Air:** \>10°C (10-day smooth)<br>**Soil:** \>5°C (7-day smooth) | **NDVI slope:** \>0.003 (7-day rolling) | Exits dormancy window | NDVI slope \<0.003 (plateaus) |

------------------------------------------------------------------------

### Days-Based Logic (Non-Dormancy Periods)

#### Post-Harvest Regrowth (after any harvest)

| Stage | Entry Conditions | Temperature | NDVI | Duration Logic | Notes |
|------------|------------|------------|------------|------------|------------|
| **Post-Cut Recovery** | `days_since_harvest 0-5` | N/A | N/A | Fixed: 5 days | High Reco phase |
| **Seedling Regrowth** | `days_since_harvest 6-10` | N/A | N/A | Fixed: 5 days | Early regrowth |
| **Vegetative** | `days_since_harvest 11-30` | N/A | N/A | Fixed: 20 days | Active growth |
| **Bud** | `days_since_harvest 31-35` | N/A | N/A | Fixed: 5 days | Bud formation |
| **Flowering** | `days_since_harvest 36-65` | N/A | N/A | Fixed: 30 days | Bloom period |
| **Maturity** | `days_since_harvest >65` | N/A | N/A | Until next harvest | Harvest window |

#### Post-Dormancy Spring (NO harvest yet - first year)

| Stage | Entry Conditions | Temperature | NDVI | Duration Logic | Notes |
|------------|------------|------------|------------|------------|------------|
| **Spring Growth** | `days_since_dormancy_exit 0-10` | Temp-triggered exit | Increasing | Fixed: 10 days | Early green-up |
| **Vegetative** | `days_since_dormancy_exit 11-30` | N/A | N/A | Fixed: 20 days | Canopy development |
| **Bud** | `days_since_dormancy_exit 31-35` | N/A | N/A | Fixed: 5 days | First reproductive stage |
| **Flowering** | `days_since_dormancy_exit 36-65` | N/A | N/A | Fixed: 30 days | Spring bloom |
| **Maturity** | `days_since_dormancy_exit >65` | N/A | N/A | Until first harvest | First harvest window |

#### True Establishment (before first dormancy window)

| Stage | Entry Conditions | Temperature | NDVI | Duration Logic | Notes |
|------------|------------|------------|------------|------------|------------|
| **Germination** | `days_since_planting 0-10` | N/A | N/A | Fixed: 10 days | Emergence |
| **Seedling** | `days_since_planting 11-30` | N/A | N/A | Fixed: 20 days | Early establishment |
| **Vegetative** | `days_since_planting 31-72` | N/A | N/A | Fixed: 42 days | Active growth |
| **Flowering** | `days_since_planting 73-93` | N/A | N/A | Fixed: 21 days | Fall bloom |
| **Seed Production** | `days_since_planting 94-128` | N/A | N/A | Fixed: 35 days | Seed set |
| **Maturity** | `days_since_planting >128` | N/A | N/A | Until dormancy | Pre-winter maturation |

------------------------------------------------------------------------

## Variable Definitions

### Temperature Variables

-   **`soil_temp_smooth`**: 7-day rolling mean of daily average soil
    temperature (°C) from ERA5
-   **`air_temp_smooth_7`**: 7-day rolling mean of daily average air
    temperature (°C) from ERA5
-   **`air_temp_smooth_10`**: 10-day rolling mean of daily average air
    temperature (°C) from ERA5

### NDVI Variables

-   **`ndvi_smooth`**: 7-day rolling mean of NDVI (smoothed values)
-   **`ndvi_slope`**: 7-day rolling OLS slope of NDVI (units/day),
    indicates rate of greening/senescence

### Day Counters

-   **`days_since_planting`**: Days since current planting date
-   **`days_since_last_harvest_same_crop`**: Days since most recent
    harvest of the same crop (NA if no harvest yet)
-   **`days_since_dormancy_exit`**: Days since dormancy window exited
    (alfalfa only, NA if not exited)

### Date Variables

-   **`month(date)`**: Month number (1-12)
-   **`day(date)`**: Day of month (1-31)
-   **Calendar gates**: Prevent false spring triggers (e.g., March 20+)
    or force dormancy (Dec 1)

------------------------------------------------------------------------

## Key Differences: Cover Crop vs Alfalfa

| Feature | Cover Crop (Winter Rye) | Alfalfa Hay |
|------------------|-----------------------------------|-------------------|
| **Primary Driver** | Temperature + NDVI for dormancy transitions;<br>Day-count for warm-season planted | Temperature for dormancy override;<br>Day-count for everything else |
| **Dormancy Entry** | Soil \<5°C sustained 10 days OR<br>Soil \<2°C immediate OR<br>Forced Dec 1 | NDVI \<0.20 OR<br>Air \<0°C OR<br>Soil \<4°C |
| **Spring Exit** | Soil \>5°C OR Air_10 \>10°C<br>+ NDVI slope \>0.001<br>+ March 20+ gate<br>+ 10 consecutive days | Air_10 \>5°C OR Soil \>5°C<br>+ NDVI \>0.20<br>+ After Feb 1<br>+ 30+ days dormant minimum |
| **Harvest Effects** | N/A (terminated, not harvested) | Resets to post_cut_recovery<br>Uses days_since_harvest for staging |
| **Vernalization** | Required for flowering<br>(warm-season planted rye lacks it) | Not required<br>(perennial, photoperiod-driven) |
| **Fall Hardening** | Soil 5-15°C, 3 consecutive days | Air 4.4-10°C, 3 consecutive days<br>(Oct-Nov only) |

------------------------------------------------------------------------

## Sustained Logic Implementation

### Cover Crop Day Counters

``` r
dormancy_entry_days <- 0L    # Requires 10 consecutive days soil <5°C
spring_exit_days    <- 0L    # Requires 10 consecutive days spring conditions
hardening_days      <- 0L    # Requires 3 consecutive days hardening conditions
```

**How it works:** - Counter increments when conditions met - Counter
resets to 0 when conditions broken - Transition only when counter
reaches threshold - Prevents single warm/cold days from triggering false
transitions

### Alfalfa Day Counters

``` r
hardening_days <- 0L         # Requires 3 consecutive days air 4.4-10°C
```

**Additional checks:** - Minimum 30 days in dormancy before spring exit
allowed - Dormancy entry/exit dates tracked explicitly - Uses
`days_in_dormancy` counter for spring exit eligibility

------------------------------------------------------------------------

## crop_stage_simple Mapping

### Cover Crop Stages → Simplified Categories

| crop_stage | crop_stage_simple | Flux Characteristics |
|------------------|-------------------------|-----------------------------|
| germination/emergence | early | Low NEE, low GPP |
| seedling | early | Low NEE, low GPP |
| fall tillering | vegetative | Moderate NEE uptake, increasing GPP |
| fall hardening | vegetative | Moderate NEE, root allocation |
| rapid growth | vegetative | High NEE uptake, high GPP |
| mature | mature | Declining NEE uptake, senescence |
| dormant | dormant | Near-zero NEE/GPP, low Reco |
| spring growth | regrowth\* | Increasing NEE uptake, greening<br>\*Only if harvest occurred |
| spring growth | vegetative | First-year post-dormancy (no harvest) |
| terminated | fallow | Soil respiration only |

### Alfalfa Stages → Simplified Categories

| crop_stage | crop_stage_simple | Flux Characteristics |
|------------------|-------------------------|-----------------------------|
| germination | early | Low NEE, establishment |
| seedling | early | Low NEE, low GPP |
| post_cut_recovery | recovery | HIGH Reco (+6 to +8), carbon source |
| seedling_regrowth | regrowth | Moderate Reco (+2 to +3), early regrowth |
| spring growth | regrowth\* OR vegetative\*\* | \*After harvest<br>\*\*Before first harvest |
| vegetative | vegetative | Moderate NEE uptake, active growth |
| bud | reproductive | Moderate NEE, reproductive allocation |
| flowering | reproductive | Moderate NEE, peak bloom |
| maturity | mature | Declining NEE, harvest-ready |
| fall hardening | vegetative | Resource allocation to roots |
| dormant | dormant | Near-zero fluxes |

------------------------------------------------------------------------

## Hudson Valley Climate Context

**Location:** Hudson Valley, NY (Zone 5b/6a)

**Typical Timeline:** - **Fall tillering end:** Late October (soil drops
to 15°C) - **Fall hardening:** Late Oct - Mid Nov (soil 5-15°C) -
**Dormancy entry:** Late November (soil \<5°C sustained) - **Dormancy
duration:** \~120-140 days (late Nov - late March) - **Spring growth:**
Late March - Early April (soil \>5°C sustained after March 20) - **First
alfalfa harvest:** Early June (60-70 days post-dormancy) - **Second
alfalfa harvest:** Early August (60-65 days post-first-harvest)

**Temperature Ranges (ERA5):** - **Winter (Dec-Feb):** Soil -3 to 2°C,
Air -20 to 10°C - **Spring (Mar-May):** Soil 3 to 15°C, Air 0 to 25°C -
**Summer (Jun-Aug):** Soil 15 to 25°C, Air 15 to 35°C - **Fall
(Sep-Nov):** Soil 5 to 20°C, Air 5 to 25°C

------------------------------------------------------------------------

## Validation Checklist

Use this checklist to verify staging is working correctly:

### Cover Crop

-   [ ] Spring growth NOT before March 20
-   [ ] Spring growth starts when soil \>5°C (not 3-4°C)
-   [ ] Dormancy entry when soil \<5°C sustained (\~late Nov)
-   [ ] Fall hardening appears in late October
-   [ ] No flowering/seed production in July/August planted rye
-   [ ] Warm-season planted: germination → seedling → vegetative →
    mature → fall hardening → dormant

### Alfalfa

-   [ ] Spring growth (pre-first-harvest) → vegetative, NOT regrowth
-   [ ] Post-harvest: recovery (days 0-5) → regrowth (days 5-10) →
    vegetative (days 11+)
-   [ ] Recovery only when harvest occurred (check
    days_since_last_harvest_same_crop is not NA)
-   [ ] Dormancy entry Oct-Nov (air \<0°C OR soil \<4°C OR NDVI \<0.20)
-   [ ] Spring exit after Feb 1 + 30 days minimum dormancy
-   [ ] Bud appears \~35 days post-dormancy or post-harvest

------------------------------------------------------------------------

## References

**Temperature Calibration:** - Hudson Valley, NY climate data (ERA5) -
Observed spring/fall transition dates (2018-2021)

**Phenology Sources:** - Alfalfa: UC Davis Extension, Cornell CALS Small
Grains Guide - Winter Rye: USDA PLANTS Database, Cornell Extension -
Vernalization requirements: Salisbury & Ross (1992) Plant Physiology

**Implementation:** -
`merge_management_data_covercropsandperennials_v3_6_0.r` - Lines 36-393
(cover crop function) - Lines 395-648 (alfalfa function) - Lines
1020-1080 (crop_stage_simple mapping)
