![](kcl.png)
![](perron.png)
![](murdoch.png)

![](projectMine.png)
![](mnda.png)


# Guide to using the App

## Contacts

- Tom Spargo <thomas.spargo@kcl.ac.uk>
- Guy Hunt <guy.hunt@kcl.ac.uk>
- Alfredo Iacoangeli <alfredo.iacoangeli@kcl.ac.uk>

## Phenotypic Clustering

Phenotypic file upload, clustering settings and figures of the clustering results.

The phenotypic file must contain the following columns:

CTYDEL                               - Diagnostic delay (in years)

Country_of_Diagnosis                 - The country where the patient was diagnosed in [ISO ALPHA-2 format](https://www.iso.org/obp/ui/#search)

Sex_at_birth                         - Sex (1 = Male, 2 = Female)

Site_of_Onset                        - Site of onset (0 = Non-bulbar, 1 =  Bulbar)

Age_at_onset_years                   - Age at symptom onset (in years)

Time_to_death_or_last_followup_years - Disease duration until death or censoring (in years)

Phenotype_1                          - Clinical phenotype (0= Other, 1=ALS)

Phenotype_2                          - Clinical phenotype (0= Other, 1=PLS)

Phenotype_3                          - Clinical phenotype (0= Other, 1=PMA)

As diagnostic delay varies with the country of diagnosis, the webserver will standardize the variable to have a mean of 0 and a standard deviation of 1 according to country of origin. If the country is one of the 13 countries originally analysed, then standardisation is relative to the mean and standard deviation from the original sample. If a given country is new for the analysis, then standardisation is performed based on the data provided when 30 or more individuals from that country are included in the data uploaded. If fewer than 30 people are sampled from that country, the standardisation is performed against the average of all countries; this is to avoid a circumstance where standardisation is performed based on data insufficient to represent the distribution of diagnostic delay in a given country.

## Phenotypic Comparison

Phenotypic comparison settings and figures of the results.

## Affiliations, Funding and Acknowledgements

1. Department of Basic and Clinical Neuroscience, Maurice Wohl Clinical Neuroscience Institute, Institute of Psychiatry, Psychology and Neuroscience, King's College London, London, SE5 9NU, UK
2. Department of Biostatistics and Health Informatics, Institute of Psychiatry, Psychology and Neuroscience, King's College London, London, UK
3. Perron Institute for Neurological and Translational Science, Nedlands, WA 6009, Australia
4. Centre for Molecular Medicine and Innovative Therapeutics, Murdoch University, Murdoch, WA 6150, Australia
5. National Institute for Health Research Biomedical Research Centre and Dementia Unit at South London and Maudsley NHS Foundation Trust and King's College London, London, UK

**ROSALIND**: We acknowledge use of the research computing facility at King’s College London, Rosalind (https://rosalind.kcl.ac.uk), which is delivered in partnership with the National Institute for Health Research (NIHR) Biomedical Research Centres at South London & Maudsley and Guy’s & St. Thomas’ NHS Foundation Trusts, and part-funded by capital equipment grants from the Maudsley Charity (award 980) and Guy’s & St. Thomas’ Charity (TR130505). 

**CREATE**: We acknowledge use of the King's Computational Research, Engineering and Technology Environment ([CREATE](https://doi.org/10.18742/rnvf-m076)).

A.I is funded by the Motor Neurone Disease Association. This study represents independent research partly funded by the National Institute for Health Research (NIHR) Biomedical Research Centre at South London and Maudsley NHS Foundation Trust and King's College London. 

The views expressed are those of the authors and not necessarily those of the NHS, the NIHR, King’s College London, or the Department of Health and Social Care.

We acknowledge the use of and are thankful to the people who have generously contributed to Project MINE and Strength ALS consortia.
