import queue
import multiprocessing as mp

from Bio import Entrez
import networkx as nx
import pandas as pd

from .node import Node, BiNode

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

        return pmids


    def efetch(self, pmid):
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        return record


    def elink(self, pmid):
        handle = Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs", from_uid=pmid)
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
        
    def search_info_csv(self, keyword, start_date, end_date):
        data = [
            self.extract_info(self.efetch(pmid))
            for pmid in self.query(keyword, start_date, end_date)
        ]
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
        return df

    def get_cited_pmids(self, pmid) -> list[str]:
        reference_pmids = []
        try:
            efetch_record = self.efetch(pmid)
            reference_list = efetch_record["PubmedArticle"][0]["PubmedData"]["ReferenceList"][0]
            for reference in reference_list["Reference"]:
                    article_id_list = reference["ArticleIdList"]
                    for pmid in article_id_list:
                        if pmid.attributes["IdType"] == "pubmed":
                                reference_pmids.append(str(pmid))
        except Exception:
            pass
        return reference_pmids
    
    def get_citing_pmids(self, pmid: str) -> list[str]:
        citing_pmids = []
        try:
            elink_record = self.elink(pmid)
            if len(elink_record[0]["LinkSetDb"]) > 0:
                citing_list = elink_record[0]["LinkSetDb"][0]["Link"]
                for citing in citing_list:
                    citing_pmids.append(citing["Id"])
        except Exception:
            pass
        return citing_pmids

    def breadth_first_exploration(self, pmid: str, max_depth: int = 5, max_iterations=50, sense="cited"):
        """
        Perform a breadth first search in the citation graph, 
        based on successive Entrez PubMed queries
        """
        if sense == "cited":
            get_successors = self.get_cited_pmids
        elif sense == "citing":
            get_successors = self.get_citing_pmids
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
            successors = get_successors(pmid)
            for successor in successors:
                if successor not in explored:
                    successor_node = Node(dict(pmid=successor, depth=depth), [], node)
                    node.children.append(successor_node)
                    frontier.put(successor_node)
            iteration += 1
        return root
    
    def bidirectional_breadth_first_exploration(self, pmid: str, max_depth: int = 10, log=True):
        q: queue.Queue = queue.Queue()
        root: BiNode = BiNode(dict(pmid=pmid, depth=0))
        q.put(root)
        explored = set()
        current_depth = 0
        while not q.empty():
            node = q.get()
            pmid = node.data["pmid"]
            explored.add(pmid)
            depth = node.data["depth"]+1
            if depth > current_depth:
                current_depth = depth
                if log:
                    print(f"Depth {current_depth}")
            if depth > max_depth:
                continue
            # Expand
            successors = self.get_cited_pmids(pmid)
            for successor in successors:
                successor_node = BiNode(dict(pmid=successor, depth=depth), [], [node])
                node.children.append(successor_node)
                if successor not in explored:
                    q.put(successor_node)
            predecessors = self.get_citing_pmids(pmid)
            for predecessor in predecessors:
                predecessor_node = BiNode(dict(pmid=predecessor, depth=depth), [node], [])
                node.parents.append(predecessor_node)
                if predecessor not in explored:
                    q.put(predecessor_node)
        return root
    
 
    def root_to_digraph(self, root: Node) -> nx.DiGraph:
        """
        Reconstruct a directed graph from the exploration node tree
        """
        G = nx.DiGraph()
        q = queue.Queue()
        current = root
        q.put(current)
        while not q.empty():
            current = q.get()
            for child in current.children:
                q.put(child)
                G.add_edge(current.data["pmid"], child.data["pmid"])
        return G
    
    def biroot_to_digraph(self, root: BiNode) -> nx.DiGraph:
        """
        """
        G = nx.DiGraph()
        q = queue.Queue()
        current = root
        q.put(current)
        while not q.empty():
            current = q.get()
            for child in current.children:
                q.put(child)
                G.add_edge(current.data["pmid"], child.data["pmid"])
            for parent in current.parents:
                q.put(parent)
                G.add_edge(parent.data["pmid"], current.data["pmid"])
        return G