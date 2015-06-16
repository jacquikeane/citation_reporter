#!/usr/bin/env python

import os, sys
import xml.etree.ElementTree as ET
from optparse import OptionParser
from datetime import datetime
from Bio import Entrez, Medline
import tempfile

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


def run_esearch(author, start_year, end_year, affiliation):
	if options.verbose:
		print "Searching for publications by", author, "affiliated with", affiliation, "between", start_year, "and", end_year
	
	
	esearch_query='(('+author+'[Author]) AND '+affiliation+'[Affiliation]) AND ("'+str(start_year)+'/1/1"[Date - Publication] : "'+str(end_year)+'"[Date - Publication])'
	#esearch_query='(('+author+'[Author]) AND ("'+str(start_year)+'/1/1"[Date - Publication] : "'+str(end_year)+'"[Date - Publication])'
	#esearch_args = shlex.split(esearch_command)
	#esearch_result=subprocess.check_output(esearch_args)
	
	Entrez.email = "Your.Name.Here@example.org"
	handle=Entrez.esearch(db="pubmed", term=esearch_query, retmax=1000)
	records=Entrez.read(handle)
	handle.close()
	if options.verbose:
		print "Found", records["Count"], "records", len(records["IdList"])
	return records["IdList"]
	

def read_author_info(line, author_id):
	words=line.strip().split(",")
	if len(words)<2 or len(words)>6:
		print "Names file must be a csv containing two to six columns for each quthor: surname (required), first name (required), middle initial (optional), affiliation(optional), ORCID ID (optional), ResearchGate ID (optional)"
		return 1
	author={}
	author["ID"]=author_id
	author["surname"]=words[0].strip()
	author["first_name"]=words[1].strip()
	if len(words)>2:
		author["middle_initials"]=words[2].strip()
	else:
		author["middle_initials"]=""
	if len(words)>3 and words[3].strip()!="":
		author["affiliation"]=words[3].strip()
	else:
		author["affiliation"]=options.affiliation
	if len(words)>4:
		author["ORCID"]=words[4].strip().replace("-","")
	else:
		author["ORCID"]=""
	if len(words)>5:
		author["Researchgate"]=words[4].strip().replace("-","")
	else:
		author["Researchgate"]=""
	author["first_initial"]=author["first_name"][0]
	author["primary_name"]=author["surname"]+" "+author["first_initial"]
	author["all_names"]=[author["primary_name"], author["surname"]+" "+author["first_name"]]
	if author["middle_initials"]!="":
		author["all_names"].append(author["surname"]+" "+author["first_initial"]+author["middle_initials"])
		author["all_names"].append(author["surname"]+" "+author["first_name"]+" "+author["middle_initials"])
		author["full_name"]=author["first_name"]+" "+author["middle_initials"]+" "+author["surname"]
	else:
		author["full_name"]=author["first_name"]+" "+author["surname"]
	return author

if __name__=="__main__":

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
	
	
	pmids={}
	authors={}
	print
	author_count=0
	for line in open(options.authorsfile, "rU"):
		author=read_author_info(line, author_count)
		if author==1:
			continue
		author_count+=1
		authors[author["ID"]]=author
		author_pmids=run_esearch(author["primary_name"], options.start, options.end, author["affiliation"])
		for pmid in author_pmids:
			if not pmid in pmids:
				pmids[pmid]=[author["ID"]]
			else:
				pmids[pmid].append(author["ID"])
	
	if options.verbose:
		print "\nFound", len(pmids), "citations with at least one author matching the input queries"
	
	for key in pmids:
		include_list.add(key)
	
	pmidlist=list(include_list.difference(exclude_list))
	
	handle = Entrez.esummary(db="pubmed", id=pmidlist[0])
	record = Entrez.read(handle)
	record[0].keys()
	
	handle = Entrez.efetch(db="pubmed", id=pmidlist, retmode="text", rettype="medline", retmax=10000)
	fp = tempfile.TemporaryFile()
	#print records
	#tmpfile=open("tmp", "w")
	fp.write(handle.read())
	#tmpfile.close()
	fp.seek(0)
	records=Medline.parse(fp)
	
	output=open(options.outputfile, "w")
	print >> output, ','.join(["Pubmed ID", "Location Identifier","Title","Authors","E-publication Date", "Publication Date", "Publication Type", "Journal", "Journal Abbreviation", "Volumne", "Issue", "Pages", "Publication Year", "Affiliated Authors"])
	
	out_count=0
	
	for record in records:
		author_matches=[]
		for au in record["AU"]:
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
				print record["TI"], "matches no authors in file. Removing..."
				continue
			else:
				if len(author_matches)==1:
					print record["TI"], "matches", len(author_matches), "author:", ', '.join(author_name_matches)
				else:
					print record["TI"], "matches", len(author_matches), "authors:", ', '.join(author_name_matches)
		else:
			if len(author_matches)==0:
				continue
		
		outlist=[]
		for key in ["PMID", "LID", "TI", "AU", "DEP", "DP", "PT", "JT", "TA", "VI", "IP", "PG"]:
			if key in record:
				if type(record[key]) is list:
					keylist="; ".join(map(str,record[key]))
					outlist.append(keylist.replace(",", ";"))
				else:
					outlist.append(str(record[key]).replace(",", ";"))
			else:
				outlist.append("")
		
		try:
			if int(record["DP"].split()[0])>=options.start and int(record["DP"].split()[0])<=options.end:
				outlist.append(str(record["DP"].split()[0]))
			else:
				outlist.append("")
		except ValueError:
			outlist.append("")
		
		outlist.append("; ".join(author_name_matches))
		
		out_count+=1
		print >> output, ",".join(outlist)
		

			
		
	fp.close()
	handle.close()
	output.close()
	print "\n", out_count, "citations with at least one author matching the input queries have been printed to", options.outputfile
	print "\nFinished\n"
