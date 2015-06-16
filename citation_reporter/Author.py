import csv
import logging


class Author(dict):
  @classmethod
  def load_csv(cls, authors_file, options):
    logger = logging.getLogger(__name__)
    authors_file.seek(0)
    authors_csv = csv.reader(authors_file)
    authors = {}
    for row_number,row in enumerate(authors_csv):
      author = Author()
      if len(row) < 2 or len(row) > 6:
        logging.warn("Row %s: Names file must be a csv containing two to six columns for each quthor: surname (required), first name (required), middle initial (optional), affiliation(optional), ORCID ID (optional), ResearchGate ID (optional)" % row_number)
        continue

      # Pad the row to 6 fields
      row.append([""*6])
      row = row[:6]
  
      author["ID"]=row_number
      author["surname"]=row[0].strip()
      author["first_name"]=row[1].strip()
      author["middle_initials"]=row[2].strip()
      if row[3].strip()!="":
        author["affiliation"]=row[3].strip()
      else:
        author["affiliation"]=options.affiliation
      author["ORCID"]=row[4].strip().replace("-","")
      author["Researchgate"]=row[4].strip().replace("-","")
  
      author["first_initial"]=author["first_name"][0]
      author["primary_name"]=author["surname"]+" "+author["first_initial"]
      author["all_names"]=[author["primary_name"], author["surname"]+" "+author["first_name"]]
      if author["middle_initials"]!="":
        author["all_names"].append(author["surname"]+" "+author["first_initial"]+author["middle_initials"])
        author["all_names"].append(author["surname"]+" "+author["first_name"]+" "+author["middle_initials"])
        author["full_name"]=author["first_name"]+" "+author["middle_initials"]+" "+author["surname"]
      else:
        author["full_name"]=author["first_name"]+" "+author["surname"]
  
      authors[row_number] = author

    return authors

  def query(self):
    return "%s[Author]" % self["primary_name"]
