# PubMed

> **General Guideline**: PubMed keyword search between two dates, export in csv file, and make a representation of that using holoviz.

A project in Object Oriented Programming teaching unit at GENIOMHE Master 1 at Université Paris-Saclay, Université d'Évry val d'Essonne.

_Due date_: 15 min oral presentation on 2023-11-15.

## Milestones

- [ ] Search PubMed between two dates for a specific keyword
  - Function for user input
  - Query constructor
  - API Query
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

## Authors and acknowledgment

- Abdou Jeylani
- Benamad Kader Houssein
- Naïa Périnelle
- Shanthosh Muruganantham
- Stéphanie Vincent
- Samuel Ortion

## References

- <https://github.com/guillaumelobet/citation-graph-pubmed>