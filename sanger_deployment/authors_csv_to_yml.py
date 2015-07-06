import csv
import yaml

authors_file = open("authors.csv", 'rb')
authors_file.readline()
authors_csv = csv.reader(authors_file)
users = []

for row_n,row in enumerate(authors_csv):
  row = [el.strip() for el in row]
  row += ["" for i in range(6)]
  surname, firstname, initials, team, orchid, google = row[:6]
  user_data = {
    'ID': str(row_n),
    'surname': surname,
    'first_name': firstname,
    'middle_initials': initials,
    'affiliation': 'Sanger',
    'ORCID': orchid,
    'google_scholar': google,
    'team': team
  }
  users.append(user_data)


print yaml.dump(users, default_flow_style=False)
