import pandas as pd

import pubmed.explorer

from Bio import Entrez

Entrez.email = "samuel.ortion@etud.univ-evry.fr"


def main():
    start_date = 2010
    end_date = 2022
    query = "cancer"
    explorer = pubmed.explorer.Explorer()
    data = list(
        map(
            lambda item: explorer.extract_info(item[1], item[0]),
            explorer.query(query, start_date, end_date),
        )
    )
    if None in data:
        data.remove(None)
    data_series = dict(
        pmid=[],
        title=[],
        doi=[],
        main_author=[],
        date=[],
        keywords=[],
    )

    for row in data:
        for key in row:
            data_series[key].append(row[key])

    df = pd.DataFrame(data_series)
    df.to_csv("tmp/test.csv")


if __name__ == "__main__":
    main()
