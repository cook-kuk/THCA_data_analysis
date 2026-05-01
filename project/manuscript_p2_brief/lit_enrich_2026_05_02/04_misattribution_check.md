# 04 Misattribution check

DOI fingerprint cross-check (year, first-author surname, journal).
Surfaces problems like the Krishnamoorthy / Landa misattribution found earlier.


## ⚠️ 2 entries with mismatches


### `Haugen2016`  (DOI: `10.1089/thy.2015.0020`)

- OpenAlex year=2015 vs bib=2016

  - bib       : {'year': '2016', 'first_author': 'Haugen', 'journal': 'thyroid'}
  - CrossRef  : {'year': '2016', 'first_author': 'Haugen', 'journal': 'thyroid'}
  - OpenAlex  : {'year': '2015', 'first_author': 'Haugen', 'journal': 'thyroid'}

### `TCGA2014`  (DOI: `10.1016/j.cell.2014.09.050`)

- first author mismatch: bib={Cancer Genome Atlas Research Network}, CrossRef=Agrawal

  - bib       : {'year': '2014', 'first_author': '{Cancer Genome Atlas Research Network}', 'journal': 'cell'}
  - CrossRef  : {'year': '2014', 'first_author': 'Agrawal', 'journal': 'cell'}
  - OpenAlex  : {'year': '2014', 'first_author': 'Agrawal', 'journal': 'cell'}