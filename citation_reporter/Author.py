import csv
import logging


class Author(object):
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

      author.ID=row_number
      author.surname=row[0].strip()
      author.first_name=row[1].strip()
      author.middle_initials=row[2].strip()
      author.middle_initials.replace(" ", "")
      if row[3].strip()!="":
        author.affiliation=row[3].strip()
      else:
        author.affiliation=options.affiliation
      author.ORCID=row[4].strip().replace("-","")
      author.Researchgate=row[4].strip().replace("-","")

      authors[row_number] = author

    return authors

  def full_name(self):
    if self.middle_initials:
      return "%s %s %s" % (self.first_name, self.middle_initials, self.surname)
    else:
      return "%s %s" % (self.first_name, self.surname)

  def primary_name(self):
    return "%s %s" % (self.surname, self.first_name[0])

  def all_names(self):
    names = [self.primary_name()]
    names.append("%s %s" % (self.surname, self.first_name))
    names.append("%s %s" % (self.first_name, self.surname))
    names.append("%s %s" % (self.first_name[0], self.surname))
    if self.middle_initials:
      names.append("%s %s%s" % (self.surname, self.first_name[0],
                                self.middle_initials))
      names.append("%s %s %s" % (self.first_name, self.middle_initials,
                                 self.surname))
      names.append("%s%s %s" % (self.first_name, self.middle_initials,
                                 self.surname))
    return names

  def is_pseudonym(self, name):
    return name in self.all_names()

  def format_pubmed_query(self):
    return "%s[Author]" % self.primary_name()
