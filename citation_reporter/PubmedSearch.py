import logging
import os

from Bio import Entrez, Medline
from datetime import datetime
from StringIO import StringIO

from citation_reporter.Author import Author

class Searcher(object):
  @classmethod
  def get_publications(cls, user, start_year=None, end_year=None):
    logger = logging.getLogger(__name__)
    if start_year == None:
      start_year = 2010
    if end_year == None:
      end_year = datetime.now().year+1
    search_query="""\
    (
      ({user_query}) AND
      ({affiliation}[Affiliation])
    ) AND (
      \"{start_year}/1/1\"[Date - Publication] : \"{end_year}\"[Date - Publication]
    )""".format(user_query=user.format_pubmed_query(),
                affiliation=user.affiliation,
                start_year=start_year, end_year=end_year)
    logger.debug("Search query: %s" % search_query)
    Entrez.email = os.environ.get("EntrezEmail", "Your.Name.Here@example.org")
    handle=Entrez.esearch(db="pubmed", term=search_query, retmax=1000)
    records=Entrez.read(handle)
    handle.close()
    logger.info("Found {count} records".format(count=len(records["IdList"])))
    return {pubmed_id: Publication(pubmed_id) for pubmed_id in records["IdList"]}

class Publication(dict):
  @classmethod
  def get_details(cls, publications):
    pubmed_ids = publications.keys()
    Entrez.email = os.environ.get("EntrezEmail", "Your.Name.Here@example.org")
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
    authors_text = "; ".join([author.user.full_name() for author in
                              self.most_likely_affiliated_authors()])
    outlist.append(authors_text)
    return outlist

  def update_authors(self, users):
    for author_string in self["AU"]:
      author_string=author_string.strip()
      if self._get_confirmed_author(author_string):
        # We already know who this author is, move on
        continue
      self.affiliated_authors.setdefault(author_string, [])
      for user in users.values():
        if self._denies_author(author_string, user.ID):
          # Someone has already said this isn't the right User, move on
          continue
        if user.is_pseudonym(author_string):
          author = Author(author_string, user.ID, user, Author.POSSIBLE)
          self.affiliated_authors[author_string].append(author)

    likely_authors = self.most_likely_affiliated_authors()
    if len(likely_authors) > 0:
      matching_authors_string = ', '.join(author.user.full_name() for author in
                                         likely_authors)
      self.logger.info("%s matches %s authors: %s" % (self["TI"],
                                                      len(likely_authors),
                                                      matching_authors_string))
    else:
      self.logger.info("%s matches no authors in file" % self["TI"])

  def _get_confirmed_author(self, author_string):
    for author in self.affiliated_authors.get(author_string, []):
      if author.confirmation_status == Author.CONFIRMED:
        return author
    return None

  def _get_possible_authors(self, author_string):
    def is_possible(author):
      return author.confirmation_status == Author.POSSIBLE
    return filter(is_possible, self.affiliated_authors.get(author_string, []))

  def _denies_author(self, author_string, user_id):
    for author in self.affiliated_authors.get(author_string, []):
      if author.ID == user_id and author.confirmtion_status == Author.DENIED:
        return True
    return False

  def has_affiliated_authors(self):
    return len(self.most_likely_affiliated_authors()) > 0

  def most_likely_affiliated_authors(self):
    likely_authors = []
    for author_string in self.affiliated_authors.keys():
      author = self._get_confirmed_author(author_string)
      if author != None:
        likely_authors.append(author)
      else:
        try:
          author = self._get_possible_authors(author_string)[0]
          likely_authors.append(author)
        except IndexError:
          # No possible authors
          continue

    return likely_authors
