---
status: >
    Accepted
description : >
    Data looks good, but there were some issues with technical replicates. See 
    notes section of README for more info. 

---

# 2021-09-12 (Run 1) Single KO Glucose-Acetate Diauxic Shift

## Purpose
This is an experiment measuring the lag-time for the single KOs of  "useless" proteins 
during the glucose-acetate diauxic shfit in minimal medium.

## Materials

### Growth Media
| **Label** | **Buffer Base** | **Carbon Source & Concentration** |
|:--:|:--:|:--:|
| ga_shift | N-C- + micronutrients | 0.6 mM glucose + 30 mM acetate|

### Strains 
| **Label** | **Parent Strain**|  **Genotype** | **Location(s)**|
|:--: | :--:| :--:| :--:|
|∆opp ∆rbs| NCM3722 | oppABCDF::attL-FRT-attR rbsDACB::attL-FRT-attR| `GC064`|
|∆mgl ∆glt| NCM3722 | mglBAC::attL-FRT-attR gltIJKL::attL-FRT-attR| `GC065`|
|∆his ∆rbs| NCM3722 | hisJQMP::attL-FRT-attR rbsDACB::attL-FRT-attR | `GC066`|
|∆opp ∆flh| NCM3722 | oppABCDF::attL-FRT-attR flhDC::attL-FRT-attR | `GC067`|
|∆his ∆flh | NCM3722 | hisJQMP::attL-FRT-attR flhDC::attL-FRT-attR | `GC068`|
|∆flh ∆glt | NCM3722 | flhDC::attL-FRT-attR gltIJKL::attL-FRT-attR | `GC069`|
|∆opp ∆glt | NCM3722 | oppABCDF::attL-FRT-attR gltIJKL::attL-FRT-attR| `GC070`|
|∆dpp ∆flh| NCM3722 | dppABCDF::attL-FRT-attR flhDC::attL-FRT-attR| `GC071`|
|∆nmp ∆rbs | NCM3722 | nmpC::attL-FRT-attR rbsDACB::attL-FRT-attR| `GC072`|
|∆rbs ∆glt | NCM3722 | rbsDACB::attL-FRT-attR gltIJKL::attL-FRT-attR | `GC073`|
|∆his ∆glt| NCM3722 | hisJQMP::attL-FRT-attR gltIJKL::attL-FRT-attR | `GC074`|
|∆dpp ∆glt| NCM3722 | dppABCDF::attL-FRT-attR gltIJKL::attL-FRT-attR| `GC075`|
|∆mgl ∆opp| NCM3722 | mglBAC::attL-FRT-attR oppABCDF::attL-FRT-attR | `GC076`|

### Instrument Settings
| Instrument | BioTek Epoch2 Microplate Reader|
|:--:| :--:|
| Temperature| 37° C|
| Shaking Speed| 1096 cpm (1mm) |
| Shaking Mode | Linear |
| Shaking Duration| 7m00s|
|Read Speed| Normal|
| Read Time | 0m32s|
| Total Interval | 7m32s |
| Number of Measurements |  96| 

### Plate Layout
| **Wells** | **Label** | **Identifier** |
|:--: | :--:  | :--: |
|C2, D2, E2 | ∆opp ∆rbs | `GC064`|
|C3, D3, E3 | ∆mgl ∆glt | `GC065` | 
|C4, D4, E4 | ∆his ∆rbs | `GC066` |
|C5, D5, E5 | ∆opp ∆flh | `GC067` |
|C6, D6, E6 | ∆his ∆flh | `GC068` |
|C7, D7, E7 | ∆flh ∆glt | `GC069` |
|C8, D8, E8 | ∆opp ∆glt | `GC070`| 
|C9, D9, E9 | ∆dpp ∆flh | `GC071` |
|C10, D10, E10 | ∆nmp ∆rbs| `GC072` |
|C11, D11, E11 | ∆rbs ∆glt| `GC073` |
|F3, F4, F5 | ∆his ∆glt | `GC074` |
|F6, F7, F8 | ∆dpp ∆glt | `GC075` |
|F9, F10, F11 | ∆mgl ∆opp | `GC076` |


## Notes & Results

Samples for ∆opp and ∆rbs aappeared to show *multiple* diauxic shifts, which doesn't 
make sense in this growth condition. For now, I've dropped these data from analysis.
Replicate 1 for ∆opp ∆flh was notably slower in the glucose-growth phase than 
the two other replicates, indicating some technical problem. THus, this replicate 
was dropped. Finally, the first measurement for ∆rbs ∆glt was very, very high, 
indicating that there was a measurement error. This point was dropped. 

### Lag Time Inference

| **Strain** | **Glucose growth rate, µ [hr]** | **Acetate growth rate, µ [per hr]** | **Lag Time, δ [hr]** | 
|:--: |:--:| :--: | :--: |
|∆opp ∆rbs | Not Determined | Not Determined | Not Determined |
|∆mgl ∆glt | 0.71 ± 0.06|  0.52 ± 0.03| 2.09 ± 0.09|
|∆his ∆rbs | Not Determined| Not Determined |  Not Determined |
|∆opp ∆flh | 0.83 ± 0.04| 0.33 ± 0.098| 2.2 ± 0.5|
|∆his ∆flh | 0.89 ± 0.04| 0.43 ± 0.03| 2.2 ± 0.1|
|∆flh ∆glt | 0.83 ± 0.03| 0.41 ± 0.02| 2.12 ± 0.03|
|∆opp ∆glt | 0.76 ± 0.03| 0.51 ± 0.01| 2.40 ± 0.02| 
|∆dpp ∆flh | 0.88 ± 0.02| 0.410 ± 0.002| 2.25 ± 0.02|    
|∆nmp ∆rbs | 0.85 ± 0.03| 0.369 ± 0.007| 2.427 ± 0.004| 
|∆rbs ∆glt | 0.79 ± 0.03| 0.37 ± 0.02| 2.48 ± 0.09|
|∆his ∆glt | 0.72 ± 0.03| 0.534 ± 0.002| 1.97 ± 0.04|
|∆dpp ∆glt | 0.74 ± 0.04| 0.468 ± 0.006| 2.5 ± 0.2|
|∆mgl ∆opp | 0.57 ± 0.04| 0.38 ± 0.02| 1.62 ± 0.09|


### Plots

![](output/2021-09-13_r1_SingleKO_glucose-acetate_shift_plot.png)

## Protocol 
1. Precultures were prepared by inoculating 3 mL of N-C- minimal medium (10 mM glucose) + 3 µL of LB with a single colony from a fresh (< 2 week old) plate.
2. Precultures were grown for 5 hours at 37° C with 250 rpm shaking until an OD 
of ≈ 0.5 was reached.
3. 1 mL of each preculture was transferred to a eppendorf tube and centrifuged 
 at 13 xg for 60 seconds. Supernatant was removed.
4. The cell pellet was resuspended with 1 mL of N-C- + 0.6 mM glucose + 30 mM acetate.
This washing step was repeated once more. 
5. THe cell pellet was resuspended in 1 mL of N-C- + 0.6 mM glucose + 30 mM acetate 
and was further diluted 1:20 into fresh shift medium prewarmed to 37° C.
4. A fresh 96 well plate was filled with water in blank wells. The remaining wells 
were filled with 200 µL of diluted and mixed cultures as appropriate and described in 
the section "Plate Layout".
5. The lid of the plate was loosely sealed to the plate by applying 4 strips of 
lab tape to the sides, preventing grinding of the plate while shaking. 
6. Plate was placed in the BioTek Epoch2 Plate reader and a kinetic cycle was begun 
as described in "Instrument Settings".
7. Data was saved, backed-up, exported, and analyzed using the `processing.py` and 
`analysis.py` Python scripts.
