import csv
import logging
import yaml

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

  def to_dict(self):
    data = {
      "pubmed_string": self.pubmed_string,
      "user_id": self.user_id,
      "user": self.user.to_dict(),
      "confirmation_status": self.confirmation_status
    }
    return data

  @classmethod
  def from_dict(cls, author_data):
    return Author(
                  author_data["pubmed_string"],
                  author_data["user_id"],
                  User.from_dict(author_data["user"]),
                  author_data["confirmation_status"],
                 )

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
      if len(row) < 3 or len(row) > 7:
        logging.warn("Row %s: Names file must be a csv containing three to seven columns for each quthor: user_id (required), surname (required), first name (required), middle initial (optional), affiliation(optional), ORCID ID (optional), ResearchGate ID (optional)" % row_number)
        continue

      # Pad the row to 6 fields
      row += [""]*7
      row = row[:7]

      user.ID=row[0].strip()
      if user.ID in users:
        logging.warn("Row %s: Could not add user, user_id (first column) must be unique" % row_number)
        continue
      user.surname=row[1].strip()
      user.first_name=row[2].strip()
      user.middle_initials=row[3].strip()
      user.middle_initials.replace(" ", "")
      if row[4].strip()!="":
        user.affiliation=row[4].strip()
      else:
        user.affiliation=options.affiliation
      user.ORCID=row[5].strip().replace("-","")
      user.Researchgate=row[5].strip().replace("-","")

      users[row_number] = user

    return users

  def full_name(self):
    if self.middle_initials:
      return "%s %s %s" % (self.first_name, self.middle_initials, self.surname)
    else:
      return "%s %s" % (self.first_name, self.surname)

  def primary_name(self):
    return "%s %s" % (self.surname, self.first_name[0])

  def ordered_name(self):
    # Name used for sorting users
    if self.middle_initials:
      return "%s, %s %s" % (self.surname, self.first_name, self.middle_initials)
    else:
      return "%s, %s" % (self.surname, self.first_name)

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
      names.append("%s%s %s" % (self.first_name[0], self.middle_initials,
                                 self.surname))
    return names

  def is_pseudonym(self, name):
    return name in self.all_names()

  def format_pubmed_query(self):
    terms = ["%s[Author]" % name for name in self.all_names()]
    return " OR ".join(terms)

  def to_dict(self):
    user_data = {
             "ID": self.ID,
             "surname": self.surname,
             "first_name": self.first_name,
             "middle_initials": self.middle_initials,
             "affiliation": self.affiliation,
             "ORCID": self.ORCID,
             "Researchgate": self.Researchgate
    }
    return user_data

  @classmethod
  def to_yaml(self, users):
    data = []
    for user in users.values():
      user_data = user.to_dict()
      data.append(user_data)
    return yaml.dump(data, default_flow_style=False)

  @classmethod
  def from_dict(cls, user_data):
    user = User()
    user.ID = str(user_data.get("ID", ""))
    user.surname = user_data.get("surname", "")
    user.first_name = user_data.get("first_name", "")
    user.middle_initials = user_data.get("middle_initials", "")
    user.affiliation = user_data.get("affiliation", "")
    user.ORCID = user_data.get("ORCID", "")
    user.Researchgate = user_data.get("Researchgate", "")
    return user

  @classmethod
  def from_yaml(self, yaml_data):
    data = yaml.load(yaml_data)
    users = {}
    for user_data in data:
      users[user_data['ID']] = User.from_dict(user_data)
    return users
