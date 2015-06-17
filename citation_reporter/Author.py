import csv
import logging

class Author(object):
  """An Author is used to link a publication to a potential User who contributed
  to it

  The object includes the string returned by PubMed to describe the user, a
  user_id for a User object (if available) and whether this link has been
  confirmed or denied by a human"""

  CONFIRMED='confirmed'
  POSSIBLE='possible'
  DENIED='denied'

  def __init__(self, pubmed_string, user_id=None, user=None,
               confirmation_status=None):
    confirmation_status = Author.POSSIBLE if confirmation_status == None else confirmation_status
    self.pubmed_string = pubmed_string
    self.user_id = user_id
    self.user = user
    self.confirmation_status = confirmation_status

class User(object):
  """A User is an object representing a person we're interested in.  This object
  has methods which can be used to evaluate whether the user is an Author of a
  Publication."""
  @classmethod
  def load_csv(cls, users_file, options):
    logger = logging.getLogger(__name__)
    users_file.seek(0)
    users_csv = csv.reader(users_file)
    users = {}
    for row_number,row in enumerate(users_csv):
      user = User()
      if len(row) < 2 or len(row) > 6:
        logging.warn("Row %s: Names file must be a csv containing two to six columns for each quthor: surname (required), first name (required), middle initial (optional), affiliation(optional), ORCID ID (optional), ResearchGate ID (optional)" % row_number)
        continue

      # Pad the row to 6 fields
      row.append([""*6])
      row = row[:6]

      user.ID=row_number
      user.surname=row[0].strip()
      user.first_name=row[1].strip()
      user.middle_initials=row[2].strip()
      user.middle_initials.replace(" ", "")
      if row[3].strip()!="":
        user.affiliation=row[3].strip()
      else:
        user.affiliation=options.affiliation
      user.ORCID=row[4].strip().replace("-","")
      user.Researchgate=row[4].strip().replace("-","")

      users[row_number] = user

    return users

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
