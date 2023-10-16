import pandas as pd
from Bio import Entrez


class Explorer:
    def __init__(self):
        pass

    def query_term(self, keyword, start_date, end_date):
        """
        Construct PubMed Query
        """
        term = f"({keyword}[Title/Abstract]) AND {start_date}:{end_date}[PDATE]"
        return term

    def query(self, keyword, start_date, end_date, maxret=None):
        """
        Fetch publication records for keyword
        """
        # Retrieve PMCIDs
        handle = Entrez.esearch(
            db="pubmed",
            term=self.query_term(keyword, start_date, end_date),
            maxret=maxret,
        )
        record = Entrez.read(handle)
        handle.close()

        pmids = record["IdList"]

        # Fetch records
        for pmid in pmids:
            handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
            record = Entrez.read(handle)
            handle.close()
            yield pmid, record

    def extract_info(self, record, pmid):
        medline_citation = record["PubmedArticle"][0]["MedlineCitation"]
        keywords = list(
            map(lambda keyword: str(keyword), medline_citation["KeywordList"])
        )
        title = medline_citation["Article"]["ArticleTitle"]
        date = medline_citation["Article"]["ArticleDate"]
        if len(date) == 0:
            print("No date")
            date = None
        else:
            date = date[0]
            date = f"{date['Year']}-{date['Month']}-{date['Day']}"
        author_list = medline_citation["Article"]["AuthorList"]
        author = f"{author_list[0]['LastName']}, {author_list[0]['ForeName']}"
        return pd.DataFrame(
            dict(
                pmid=pmid, title=title, main_author=author, date=date, keywords=keywords
            )
        )
