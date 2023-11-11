import pubmed.explorer


# def test_explorer_query_term():
#     start_date = 2020
#     end_date = 2022
#     query = "single cell"

#     explorer = pubmed.explorer.Explorer()
#     print(explorer.query_term(start_date, end_date, query))


# def test_explorer_query_records():
#     start_date = 2020
#     end_date = 2022
#     query = "single cell"

#     explorer = pubmed.explorer.Explorer()
#     print(len(list(explorer.query(start_date, end_date, query, maxret=5))))


def test_to_pandas():
    start_date = 2020
    end_date = 2022
    query = "single cell"

    explorer = pubmed.explorer.Explorer()

    for pmid, record in explorer.query(start_date, end_date, query, maxret=2):
        print(explorer.extract_info(record, pmid))
