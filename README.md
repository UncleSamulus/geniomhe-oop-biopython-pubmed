# PubMed

PubMed keyword search between two dates, export in csv file, and make a representation of that using holoviz.

## Goals

- [ ] Search PubMed between two dates for a specific keyword
  - Function for user input
  - Query constructor
  - API Query
  - Save records as list
  - Save as CSV
    Columns:
    - PMCID
    - Titre;
    - Keywords;
    - Authors;
    - Date;
    - Journal;
- [ ] Representations
  - [ ] Number of publication per year;
  - [ ] Network of keywords
  - [ ] Network of authors

## Install

```bash
pip install -e .
```

## Unit testing

```bash
pytest --pyargs pubmed
```
