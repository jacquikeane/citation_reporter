import csv
import logging
import os
import yaml

from Bio import Entrez, Medline
from datetime import datetime
from StringIO import StringIO

from citation_reporter.Author import Author

class Searcher(object):
  @classmethod
  def get_pubmed_ids_for_user(cls, user, start_year=None, end_year=None):
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
    logging.debug("Search query: %s" % search_query)
    Entrez.email = os.environ.get("EntrezEmail", "Your.Name.Here@example.org")
    handle=Entrez.esearch(db="pubmed", term=search_query, retmax=1000)
    records=Entrez.read(handle)
    handle.close()
    logging.info("Found {count} records for {name}".format(count=len(records["IdList"]),
                                                           name=user.full_name()))
    pubmed_ids = records["IdList"]
    return pubmed_ids

class Publications(dict):
  def __init__(self, publications=None):
    publications = [] if publications == None else publications
    for publication in publications:
      pubmed_id = publication['PMID']
      self[pubmed_id] = publication

  @classmethod
  def from_pubmed_ids(cls, pubmed_ids):
    Entrez.email = os.environ.get("EntrezEmail", "Your.Name.Here@example.org")
    handle = Entrez.efetch(db="pubmed", id=pubmed_ids, retmode="text",
                           rettype="medline", retmax=10000)
    records = Medline.parse(handle)
    publications = Publications()
    for record in records:
      pubmed_id = record['PMID']
      publications[pubmed_id] = Publication(pubmed_id, record)
    return publications

  def to_csv(self, output_file):
    output_csv = csv.writer(output_file, lineterminator='\n') # remove line terminator to be windows friendly
    output_csv.writerow(Publication.format_header_row())
    for publication in self.values():
      output_csv.writerow(publication.format_row())

  def to_yaml(self):
    data = [publication.to_dict() for publication in self.values()]
    return yaml.dump(data, default_flow_style=False)

  @classmethod
  def from_yaml(cls, publications_yaml):
    data = yaml.load(publications_yaml)
    publications = Publications()
    data = [] if data == None else data # publications file might be empty
    for publication_data in data:
      pubmed_id = publication_data["Pubmed ID"]
      publication = Publication(pubmed_id)
      for internal_key, external_key in zip(publication._internal_keys(),
                                            publication.format_header_row()):
        value = publication_data.get(external_key)
        if not value is None:
          publication[internal_key] = value

      affiliated_authors_data = publication_data.get('affiliated_authors_data', {})
      publication.affiliated_authors = publication._affiliated_authors_from_dict(affiliated_authors_data)
      publication.confirmation_status = publication_data.get('confirmation_status', Publication.POSSIBLE)

      publications[pubmed_id] = publication
    logging.debug("Loaded %s publications" % len(publications))
    return publications

  @classmethod
  def merge(cls, old_publications, new_publications):
    """If a publication already exists in old_publications, do not update it
    with data from new_publications"""
    publications = new_publications
    publications.update(old_publications)
    return publications

  def denied(self):
    return Publications([publication for publication in self.values() if
            publication.confirmation_status == Publication.DENIED])

  def not_denied(self):
    return Publications([publication for publication in self.values() if
            publication.confirmation_status != Publication.DENIED])

class Publication(dict):

  CONFIRMED='confirmed'
  POSSIBLE='possible'
  DENIED='denied'

  def __init__(self, pubmed_id, data=None):
    data = {} if data == None else data
    self.pubmed_id = pubmed_id
    for key, value in data.items():
      self[key] = value
    affiliated_authors_data = data.get('affiliated_authors_data', {})
    self.affiliated_authors = self._affiliated_authors_from_dict(affiliated_authors_data)
    self.confirmation_status = Publication.POSSIBLE

  def deny(self):
    """Deny that the publication is relevant

    Set the status of all authors to denied and update the status of the
    publication itself"""
    self.confirmation_status = Publication.DENIED
    for author_string, authors in self.affiliated_authors.items():
      for author in authors:
        author.confirmation_status = Author.DENIED

  @classmethod
  def _affiliated_authors_from_dict(cls, affiliated_authors_data):
    affiliated_authors = {}
    for author_string, authors_list in affiliated_authors_data.items():
      authors = [Author.from_dict(author_data) for author_data in authors_list]
      affiliated_authors[author_string] = authors
    return affiliated_authors

  @classmethod
  def format_header_row(cls):
    return ["Pubmed ID", "Location Identifier","Title","Authors",
            "E-publication Date", "Publication Date", "Publication Type",
            "Journal", "Journal Abbreviation", "Volume", "Issue",
            "Pages", "Publication Year", "Affiliated Authors"]

  @classmethod
  def _internal_keys(cls):
    return ["PMID", "LID", "TI", "AU", "DEP", "DP", "PT", "JT", "TA", "VI",
            "IP", "PG"]

  def format_row(self):
    outlist=[]
    for key in self._internal_keys():
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
      for user in users.values():
        if self._denies_author(author_string, user.ID):
          # Someone has already said this isn't the right User, move on
          continue
        if self._user_already_author(author_string, user.ID):
          # We already know this is a user, move on
          continue
        if user.is_pseudonym(author_string):
          author = Author(author_string, user.ID, user, Author.POSSIBLE)
          self.affiliated_authors.setdefault(author_string, []).append(author)

    if not self.has_affiliated_authors():
      self.deny()

  def update_author_status(self, author_string, user_id, status):
    if not status in [Author.DENIED, Author.CONFIRMED]:
      raise KeyError("Status can only be set to %s or %s" % (Author.DENIED,
                                                             Author.CONFIRMED))
    if user_id == 'all' and status == Author.DENIED:
      for author in self.affiliated_authors[author_string]:
        author.confirmation_status = status
    elif user_id != "all" and status == Author.CONFIRMED:
      for author in self.affiliated_authors[author_string]:
        if author.user_id == user_id:
          author.confirmation_status = Author.CONFIRMED
        else:
          author.confirmation_status = Author.DENIED
    elif status == Author.DENIED:
      for author in self.affiliated_authors[author_string]:
        if author.user_id == user_id:
          author.confirmation_status = Author.DENIED

    if not self.has_affiliated_authors():
      self.confirmation_status = Publication.DENIED
    elif self.has_confirmed_author():
      self.confirmation_status = Publication.CONFIRMED
    else:
      self.confirmation_status = Publication.POSSIBLE

  def _user_already_author(self, author_string, user_id):
    for author in self.affiliated_authors.get(author_string, []):
      if author.user_id == user_id:
        return True
    return False

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
      if author.user_id == user_id and author.confirmation_status == Author.DENIED:
        return True
    return False

  def has_confirmed_author(self):
    for authors in self.affiliated_authors.values():
      for author in authors:
        if author.confirmation_status == Author.CONFIRMED:
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
    data = {}
    # Most data is stored in a dictionary which mirrors the format used by the
    # pubmed API.  This is not very human readable so when we store the data on
    # disk we use nicer human readable keys, as used in the CSV output.
    internal_keys = self._internal_keys()
    external_keys = self.format_header_row()
    for internal_key, external_key in zip(internal_keys,
                                          external_keys):
      value = self.get(internal_key)
      if not value is None:
        data[external_key] = value
    affiliated_authors_data = {}
    for author_string, authors in self.affiliated_authors.items():
      affiliated_authors_data[author_string] = [author.to_dict() for author in
                                                authors]
    data['affiliated_authors_data'] = affiliated_authors_data
    data['confirmation_status'] = self.confirmation_status
    return data
