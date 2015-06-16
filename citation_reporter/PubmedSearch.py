import logging

from Bio import Entrez
from datetime import datetime

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

class Publication(object):
  def __init__(self, pubmed_id):
    self.pubmed_id = pubmed_id
