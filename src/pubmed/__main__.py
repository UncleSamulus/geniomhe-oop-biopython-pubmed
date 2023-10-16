import pandas as pd

import pubmed.explorer


from Bio import Entrez

Entrez.email = "samuel.ortion@etud.univ-evry.fr"


def main():
    start_date = 2020
    end_date = 2022
    query = "single cell"
    explorer = pubmed.explorer.Explorer()
    df_parts = list(
        map(
            lambda item: explorer.extract_info(item[1], item[0]),
            explorer.query(query, start_date, end_date),
        )
    )
    df = pd.concat(df_parts)
    df.to_csv("test.csv")


if __name__ == "__main__":
    main()
