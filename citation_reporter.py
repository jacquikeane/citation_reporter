#!/usr/bin/env python

import csv
import os, sys
import xml.etree.ElementTree as ET
from optparse import OptionParser
from datetime import datetime
import logging
import tempfile

logging.basicConfig(level=logging.INFO)

from citation_reporter.Author import User
from citation_reporter.PubmedSearch import Searcher, Publication

START_YEAR=2010
END_YEAR=datetime.now().year+1
DEFAULT_AFFILIATION="Sanger"


def main():
  usage = "usage: %prog [options]"
  parser = OptionParser(usage=usage)
        
        
  parser.add_option("-o", "--output", action="store", dest="outputfile", help="output file name (in csv format)", type="string", metavar="FILE", default="")
  parser.add_option("-u", "--users", action="store", dest="usersfile", help="csv file containing two to six columns for each user: surname (required), first name (required), middle initial (optional), affiliation(optional), ORCID ID (optional), ResearchGate ID (optional)", type="string", metavar="FILE", default="")
  parser.add_option("-p", "--publications", action="store", dest="publicationsfile",
                    help="publications file to be updated (in yaml format as output by this script)", type="string", metavar="FILE", default="")
  parser.add_option("-s", "--start_year", action="store", dest="start", help="Year to start search from [default = %default]", type="int", metavar="YEAR", default=START_YEAR)
  parser.add_option("-e", "--end_year", action="store", dest="end", help="Year to end search [default = present (<%default)]", type="int", metavar="YEAR", default=END_YEAR)
  parser.add_option("-A", "--affiliation", action="store", dest="affiliation", help="Default affiliation (applies to any user without an affiliation specified in the users file. [default = %default]", metavar="AFFILIATION", default=DEFAULT_AFFILIATION)
  parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="Verbose mode", metavar="AFFILIATION", default=False)
  
  
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
      publications = Publication.from_yaml(publications_input.read())
  except IOError:
    # The file is missing
    logging.info("Could not load existing publications from %s, skipping" %
                options.publicationsfile)
    publications={}
  users={}

  with open(options.usersfile, "rU") as usersfile:
    users = User.load_csv(usersfile, options)
  users_count = len(users)
  pubmed_ids = set()
  for user in users.values():
    new_pubmed_ids = Searcher.get_pubmed_ids_for_user(user, options.start, options.end)
    pubmed_ids.update(new_pubmed_ids)

  if options.verbose:
    print "\nFound", len(publications), "citations with at least one user matching the input queries"
  
  new_publications = Publication.from_pubmed_ids(list(pubmed_ids))
  publications = Publication.merge_publications(publications, new_publications)

  output=open(options.outputfile, "wb")
  csv_writer = csv.writer(output, lineterminator='\n') # remove line terminator to be windows friendly
  csv_writer.writerow(Publication.format_header_row())
  
  out_count=0
  for publication in publications.values():
    publication.update_authors(users)

    if not publication.has_affiliated_authors():
      continue
    
    publication_row = publication.format_row()
    
    out_count+=1
    csv_writer.writerow(publication_row)
    
  output.close()
  print "\n", out_count, "citations with at least one user matching the input queries have been printed to", options.outputfile
  try:
    with open(options.publicationsfile, 'w') as publications_output:
      publications_output.write(Publication.to_yaml(publications))
  except IOError:
    # The file is missing
    logging.info("Could not write existing publications to %s, skipping" %
                options.publicationsfile)
  print "\nFinished\n"
