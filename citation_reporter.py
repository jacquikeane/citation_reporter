#!/usr/bin/env python

import os, sys
import xml.etree.ElementTree as ET
from optparse import OptionParser
from datetime import datetime
import logging
import tempfile

from citation_reporter.Author import Author
from citation_reporter.PubmedSearch import Searcher, Publication

START_YEAR=2010
END_YEAR=datetime.now().year+1
DEFAULT_AFFILIATION="Sanger"


def main():
  usage = "usage: %prog [options]"
  parser = OptionParser(usage=usage)
        
        
  parser.add_option("-o", "--output", action="store", dest="outputfile", help="output file name (in csv format)", type="string", metavar="FILE", default="")
  parser.add_option("-a", "--authors", action="store", dest="authorsfile", help="csv file containing two to six columns for each author: surname (required), first name (required), middle initial (optional), affiliation(optional), ORCID ID (optional), ResearchGate ID (optional)", type="string", metavar="FILE", default="")
  parser.add_option("-i", "--include", action="store", dest="includefile", help="File containing PMIDs to add to matches (optional)", type="string", metavar="FILE", default="")
  parser.add_option("-x", "--exclude", action="store", dest="excludefile", help="File containing list of PMIDs to exclude from matches (optional)", type="string", metavar="FILE", default="")
  parser.add_option("-s", "--start_year", action="store", dest="start", help="Year to start search from [default = %default]", type="int", metavar="YEAR", default=START_YEAR)
  parser.add_option("-e", "--end_year", action="store", dest="end", help="Year to end search [default = present (<%default)]", type="int", metavar="YEAR", default=END_YEAR)
  parser.add_option("-A", "--affiliation", action="store", dest="affiliation", help="Default affiliation (applies to any author without an affiliation specified in the authors file. [default = %default]", metavar="AFFILIATION", default=DEFAULT_AFFILIATION)
  parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="Verbose mode", metavar="AFFILIATION", default=False)
  
  
  return parser.parse_args()
  
def check_command_line():
  do_exit=False
  
  if options.authorsfile=="":
    print "No authors file specified"
    do_exit=True
  if not os.path.isfile(options.authorsfile):
    print "Cannot find authors file:", options.authorsfile
    do_exit=True
  if options.includefile!="" and not os.path.isfile(options.includefile):
    print "Cannot find file:", options.includefile
    do_exit=True
  if options.excludefile!="" and not os.path.isfile(options.excludefile):
    print "Cannot find file:", options.excludefile
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
  logger = logging.getLogger(__name__)
  logger.setLevel(logging.DEBUG)

  (options, args)=main()
  check_command_line()
  
  include_list=set([])
  
  if options.includefile!="":
    for line in open(options.includefile):
      include_list.add(line.strip())
  
  exclude_list=set([])
  if options.excludefile!="":
    for line in open(options.excludefile):
      exclude_list.add(line.strip())
  
  publications={pubmed_id: Publication(pubmed_id, {'whitelist': True}) for pubmed_id in include_list}
  authors={}

  with open(options.authorsfile, "rU") as authorsfile:
    authors = Author.load_csv(authorsfile, options)
  authors_count = len(authors)
  for author in authors.values():
    publications.update(Searcher.get_publications(author, options.start,
                                                  options.end))

  if options.verbose:
    print "\nFound", len(publications), "citations with at least one author matching the input queries"
  
  publications={publication.pubmed_id: publication for publication in
                publications.values() if publication.pubmed_id not in
                exclude_list}
  
  publications = Publication.get_details(publications)
  
  output=open(options.outputfile, "w")
  print >> output, ','.join(["Pubmed ID", "Location Identifier","Title","Authors","E-publication Date", "Publication Date", "Publication Type", "Journal", "Journal Abbreviation", "Volumne", "Issue", "Pages", "Publication Year", "Affiliated Authors"])
  
  out_count=0
  
  for publication in publications.values():
    author_matches=[]
    for au in publication["AU"]:
      au=au.strip()
      for authorid in authors:
        author=authors[authorid]
        if au in author["all_names"]:
          author_matches.append(author["ID"])
          
    author_name_matches=[]
    for match in author_matches:
      author_name_matches.append(authors[match]["full_name"])
    
    if options.verbose:
      if len(author_matches)==0:
        print publication["TI"], "matches no authors in file. Removing..."
        continue
      else:
        if len(author_matches)==1:
          print publication["TI"], "matches", len(author_matches), "author:", ', '.join(author_name_matches)
        else:
          print publication["TI"], "matches", len(author_matches), "authors:", ', '.join(author_name_matches)
    else:
      if len(author_matches)==0:
        continue
    
    publication_text = publication.format()
    author_text = "; ".join(author_name_matches)
    all_text = "{publication},{author_names}".format(
                 publication=publication_text,
                 author_names=author_text)
    
    out_count+=1
    print >> output, all_text
    
  output.close()
  print "\n", out_count, "citations with at least one author matching the input queries have been printed to", options.outputfile
  print "\nFinished\n"
