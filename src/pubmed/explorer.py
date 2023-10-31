import queue

from Bio import Entrez

from .node import Node

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
            
            yield record

    def efetch(self, pmid):
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        return record

    def extract_info(self, record):
        if "PubmedArticle" in record:
            medline_citation = record["PubmedArticle"][0]["MedlineCitation"]
            pmid = medline_citation["PMID"]
            if len(medline_citation["KeywordList"]) > 0:
                keywords = list(
                    map(lambda keyword: str(keyword), medline_citation["KeywordList"][0])
                )
            else:
                keywords = None
            title = medline_citation["Article"]["ArticleTitle"]
            date = medline_citation["Article"]["ArticleDate"]
            doi = str(medline_citation["Article"]["ELocationID"][0])
            if len(date) == 0:
                date = None
            else:
                date = date[0]
                date = f"{date['Year']}-{date['Month']}-{date['Day']}"
            author_list = medline_citation["Article"]["AuthorList"]
            author = f"{author_list[0]['LastName']}, {author_list[0]['ForeName']}"
            return dict(
                pmid=pmid,
                title=title,
                doi=doi,
                main_author=author,
                date=date,
                keywords=keywords,
            )
        else:
            return None
    
    def get_cited_pmids(self, efetch_record: dict) -> list[str]:
        reference_pmids = []
        try:
            reference_list = efetch_record["PubmedArticle"][0]["PubmedData"]["ReferenceList"][0]
            for reference in reference_list["Reference"]:
                    article_id_list = reference["ArticleIdList"]
                    for pmid in article_id_list:
                        if pmid.attributes["IdType"] == "pubmed":
                                reference_pmids.append(str(pmid))
        except KeyError as e:
            pass
        except IndexError as e:
            pass
        return reference_pmids

    def breadth_first_exploration(self, pmid: str, max_depth: int = 5, max_iterations=50):
        """
        Perform a breadth first search in the citation graph, 
        based on successive Entrez PubMed queries
        """
        frontier: queue.Queue = queue.Queue()
        node = Node(dict(pmid=pmid, depth=0))
        root = node
        frontier.put(node)
        explored: set[str] = set()
        iteration: int = 0
        while not frontier.empty() and iteration < max_iterations:
            node = frontier.get()
            pmid = node.data["pmid"]
            explored.add(pmid)
            depth = node.data["depth"]+1
            if depth > max_depth:
                continue
            # Expand
            record = self.efetch(pmid)
            successors = self.get_cited_pmids(record)
            for successor in successors:
                if successor not in explored:
                    successor_node = Node(dict(pmid=successor, depth=depth), [], node)
                    node.children.append(successor_node)
                    frontier.put(successor_node)
            iteration += 1
        return root