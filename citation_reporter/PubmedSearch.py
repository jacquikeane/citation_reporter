import logging

from Bio import Entrez, Medline
from datetime import datetime
from StringIO import StringIO

class Searcher(object):
  @classmethod
  def get_publications(cls, author, start_year=None, end_year=None):
    logger = logging.getLogger(__name__)
    if start_year == None:
      start_year = 2010
    if end_year == None:
      end_year = datetime.now().year+1
    search_query="""\
    (
      ({author_query}) AND
      ({affiliation}[Affiliation])
    ) AND (
      \"{start_year}/1/1\"[Date - Publication] : \"{end_year}\"[Date - Publication]
    )""".format(author_query=author.format_pubmed_query(),
                affiliation=author['affiliation'],
                start_year=start_year, end_year=end_year)
    logger.debug("Search query: %s" % search_query)
    Entrez.email = "Your.Name.Here@example.org"
    handle=Entrez.esearch(db="pubmed", term=search_query, retmax=1000)
    records=Entrez.read(handle)
    handle.close()
    logger.info("Found {count} records".format(count=len(records["IdList"])))
    return {pubmed_id: Publication(pubmed_id) for pubmed_id in records["IdList"]}

class Publication(dict):
  @classmethod
  def get_details(cls, publications):
    pubmed_ids = publications.keys()
    handle = Entrez.efetch(db="pubmed", id=pubmed_ids, retmode="text",
                           rettype="medline", retmax=10000)
    records = Medline.parse(handle)
    for record in records:
      pubmed_id = record['PMID']
      publications[pubmed_id].update(record)
    return publications

  def __init__(self, pubmed_id, data=None):
    self.logger = logging.getLogger(__name__)
    data = {} if data == None else data
    self.pubmed_id = pubmed_id
    for key, value in data.items():
      self[key] = value
    self.affiliated_authors = {}

  @classmethod
  def format_header_row(cls):
    return ["Pubmed ID", "Location Identifier","Title","Authors",
            "E-publication Date", "Publication Date", "Publication Type",
            "Journal", "Journal Abbreviation", "Volumne", "Issue",
            "Pages", "Publication Year", "Affiliated Authors"]

  def format_row(self):
    outlist=[]
    for key in ["PMID", "LID", "TI", "AU", "DEP", "DP", "PT", "JT", "TA", "VI", "IP", "PG"]:
      value = self.get(key, "")
      if isinstance(value, list):
        value="; ".join(map(str,value))
      value = str(value).replace(",", ";")
      outlist.append(value)
    publication_date = str(self.get("DP", "").split()[0])
    outlist.append(publication_date)
    authors_text = "; ".join([author["full_name"] for author in
                              self.most_likely_affiliated_authors()])
    outlist.append(authors_text)
    return outlist

  def update_authors(self, authors):
    for author_string in self["AU"]:
      author_string=author_string.strip()
      self.affiliated_authors.setdefault(author_string, [])
      for author in authors.values():
        if author.is_pseudonym(author_string):
          self.affiliated_authors[author_string].append((author, 1.0))

    likely_authors = self.most_likely_affiliated_authors()
    if len(likely_authors) > 0:
      matching_authors_string = ', '.join(author["full_name"] for author in
                                         likely_authors)
      self.logger.info("%s matches %s authors: %s" % (self["TI"],
                                                      len(likely_authors),
                                                      matching_authors_string))
    else:
      self.logger.info("%s matches no authors in file" % self["TI"])

  def has_affiliated_authors(self):
    return len(self.most_likely_affiliated_authors()) > 0

  def most_likely_affiliated_authors(self):
    def by_probability(author_probability):
      author, probability = author_probability
      return probability
    likely_authors = []
    for author_string, author_probabilities in self.affiliated_authors.items():
      if len(author_probabilities) == 0:
        continue
      most_likely = max(author_probabilities, key=by_probability)
      author, probability = most_likely
      if probability > 0:
        likely_authors.append(author)
    return likely_authors
