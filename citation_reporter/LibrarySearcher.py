import os, sys, re, requests
from datetime import datetime
import logging

START_YEAR=2010
END_YEAR=datetime.now().year

class LibrarySearcher(object):
  @classmethod
  def get_pubmed_ids(cls, start_year=START_YEAR, end_year=END_YEAR):
    pubmed_ids = set()
    logging.info("Getting pubmed_ids from the Library")
    try:
      for year in range(start_year, end_year + 1):
        r = requests.get("http://www.sanger.ac.uk/research/publications/%s.html" %
                         year)
        new_pubmed_ids = re.findall("http://ukpmc.ac.uk/abstract/MED/(\d{8})",
                                    r.text)
        pubmed_ids.update(new_pubmed_ids)
    except Exception as e:
      logging.error("Could not fetch pubmed_ids from the Library: %s" % e)
    logging.info("Fetched %s ids from the library" % len(pubmed_ids))
    return pubmed_ids
