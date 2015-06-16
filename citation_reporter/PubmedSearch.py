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
    )""".format(author_query=author.query(), affiliation=author['affiliation'],
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
    data = {} if data == None else data
    self.pubmed_id = pubmed_id
    for key, value in data.items():
      self[key] = value

  def format(self):
    outlist=[]
    for key in ["PMID", "LID", "TI", "AU", "DEP", "DP", "PT", "JT", "TA", "VI", "IP", "PG"]:
      value = self.get(key, "")
      if isinstance(value, list):
        value="; ".join(map(str,value))
      value = str(value).replace(",", ";")
      outlist.append(value)
    publication_date = str(self.get("DP", "").split()[0])
    outlist.append(publication_date)
    return ','.join(outlist)
