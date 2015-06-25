#!/usr/bin/env python

import os, sys
import xml.etree.ElementTree as ET
from optparse import OptionParser
from datetime import datetime
import logging
import tempfile

logging.basicConfig(level=logging.INFO)

from citation_reporter.Author import User
from citation_reporter.PubmedSearch import Searcher, Publication, Publications

START_YEAR=2010
END_YEAR=datetime.now().year+1
DEFAULT_AFFILIATION="Sanger"


def main():
  usage = "usage: %prog [options]"
  parser = OptionParser(usage=usage)


  parser.add_option("-o", "--output", action="store", dest="outputfile",
                    help="output file name (in csv format)", type="string", metavar="FILE", default="")
  parser.add_option("-u", "--users", action="store", dest="usersfile",
                    help="Yaml file of the following data for each user: surname (required), first name (required), middle initial (optional), affiliation(required), ORCID ID (optional), ResearchGate ID (optional)", type="string", metavar="FILE", default="")
  parser.add_option("-p", "--publications", action="store", dest="publicationsfile",
                    help="publications file to be updated (in yaml format as output by this script)", type="string", metavar="FILE", default="")
  parser.add_option("-s", "--start_year", action="store", dest="start", help="Year to start search from [default = %default]", type="int", metavar="YEAR", default=START_YEAR)
  parser.add_option("-e", "--end_year", action="store", dest="end", help="Year to end search [default = present (<%default)]", type="int", metavar="YEAR", default=END_YEAR)
  parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                    help="Verbose mode", metavar="VERBOSE", default=False)


  return parser.parse_args()

def check_command_line():
  do_exit=False

  if options.usersfile=="":
    print "No users file specified"
    do_exit=True
  if not os.path.isfile(options.usersfile):
    print "Cannot find users file:", options.usersfile
    do_exit=True
  if options.outputfile=="":
    print "No output file specified"
    do_exit=True
  if options.start>datetime.now().year:
    print "Start year must be <", datetime.now().year
    do_exit=True
  if options.start<1900:
    print "Start year must be >1900"
    do_exit=True
  if options.end<1900:
    print "End year must be >1900"
    do_exit=True
  if options.start>options.end:
    print "Start year must be <= end year"
    do_exit=True
  
  if do_exit:
    print"Exiting"
    print
    sys.exit(1)


if __name__=="__main__":

  (options, args)=main()
  check_command_line()

  try:
    with open(options.publicationsfile, 'r') as publications_input:
      publications = Publications.from_yaml(publications_input.read())
  except IOError:
    # The file is missing
    logging.info("Could not load existing publications from %s, skipping" %
                options.publicationsfile)
    publications={}
  users={}

  with open(options.usersfile, "r") as usersfile:
    users = User.from_yaml(usersfile.read())
  users_count = len(users)
  logging.info("Loaded %s users from %s" % (users_count, options.usersfile))

  pubmed_ids = set()
  logging.info("Searching pubmed for publications")
  for user in users.values():
    new_pubmed_ids = Searcher.get_pubmed_ids_for_user(user, options.start, options.end)
    pubmed_ids.update(new_pubmed_ids)
  pubmed_id_count = len(pubmed_ids)
  logging.info("Found %s citations by searching pubmed" % pubmed_id_count)

  new_publications = Publications.from_pubmed_ids(list(pubmed_ids))
  publications = Publications.merge(publications, new_publications)

  for publication in publications.values():
    publication.update_authors(users)

  output=open(options.outputfile, "wb")
  affiliated_publications = publications.not_denied()
  affiliated_publications.to_csv(output)
  output.close()

  logging.info("%s citations with at least one user matching the input queries have been printed to %s" % (
    len(affiliated_publications), options.outputfile))
  try:
    with open(options.publicationsfile, 'w') as publications_output:
      publications_output.write(publications.to_yaml())
  except IOError:
    # The file is missing
    logging.info("Could not write existing publications to %s, skipping" %
                options.publicationsfile)
  logging.info("Finished")
