import logging
import os
import yaml

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
    """Updates self.affiliated_authors

    self.affiliated_authors is a dictionary mapping the author strings returned
    by pubmed to possible Users via a list of Author objects.  Author objects
    include the pubmed provided string, a user_id, a User object and whether or
    not a human has confirmed that this is a correct match"""
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
    """Returns a list of Author objects where a human has confirmed that one of
    our Users is actually an Author of this Publication"""
    for author in self.affiliated_authors.get(author_string, []):
      if author.confirmation_status == Author.CONFIRMED:
        return author
    return None

  def _get_possible_authors(self, author_string):
    """Returns a list of all of the possible Authors of a Publication.

    i.e. it excludes all mappings of Publication to User which a human has said
    are either definitly accurate or definitly not accurate"""
    def is_possible(author):
      return author.confirmation_status == Author.POSSIBLE
    return filter(is_possible, self.affiliated_authors.get(author_string, []))

  def _denies_author(self, author_string, user_id):
    """Checks whether a human has denied that the given user_id is not the
    reported author_string for this Publication"""
    for author in self.affiliated_authors.get(author_string, []):
      if author.user_id == user_id and author.confirmtion_status == Author.DENIED:
        return True
    return False

  def has_affiliated_authors(self):
    """Checks if any of our Users is a confirmed or possible Author of the
    Publication"""
    return len(self.most_likely_affiliated_authors()) > 0

  def most_likely_affiliated_authors(self):
    """Returns the most likely Users who may have authored the Publication.

    If a human has confirmed a User is an Author, this is returned.  Otherwise
    each of the authors returned by pubmed is matched to at most one User (or
    none if no sensible matches are possible"""
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

  def to_dict(self):
    keys = self.format_header_row()
    values = self.format_row()
    data = dict(zip(keys, values))
    affiliated_authors_data = {}
    for author_string, authors in self.affiliated_authors.items():
      affiliated_authors_data[author_string] = [author.to_dict() for author in
                                                authors]
    data['affiliated_authors_data'] = affiliated_authors_data
    return data

  @classmethod
  def to_yaml(cls, publications):
    data = [publication.to_dict() for publication in publications.values()]
    return yaml.dump(data, default_flow_style=False)
