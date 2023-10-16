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
            yield record
