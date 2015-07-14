# Citation Reporter

Citation reporter is a tool to keep track of a team's publications.  It currently only searched pubmed and the Sanger Library for publications and it matches these with a list of users in your group.  Publications by authors from the team are refered to as 'affiliated' in the UI and the code.

The tool is based on scripts developed by @simonrharris.

## Usage

Publications can be loaded through the web UI or a commandline script.  New users can only be added via the commandline script.

### Installation

```
pip install git+https://github.com/sanger-pathogens/citation_reporter.git@master
```

### scripts/citation_reporter_cli.py

```
citation_reporter_cli.py -u authors.yml -s 2014 -e 2015 -p publications.yml -o publications.csv -v
```

* `-u authors.yml` - details of authors in yml format (see examples/authors.yml)
* `-s 2014` - start year for search: only include publications from 2014 or later
* `-e 2015` - end year for search: only include publications from 2015 or earlier (defaults to next year)
* `-p publications.yml` - Yaml file to persist details of publications: will be created if it doesn't exist
* `-o publications.csv` - CSV output of publications with confirmed or potential affiliated authors
* `-v` - more verbose logs

### scripts/citation_reporter_web.py

```
CITATION_REPORTER_PORT=8080 citation_reporter_web.py
```

Environment Variables:

* `CITATION_REPORTER_PUBLICATIONS` - file to persist details of publications (defaults to 'publications.yml')
* `CITATION_REPORTER_USERS` - file containing details of users (defaults to 'authors.yml')
* `BOOTSTRAP_URL` - URL for a version of bootstrap (defaults to "//maxcdn.bootstrapcdn.com/bootstrap/3.3.5")
* `JQUERY_URL` - URL for jquery (defaults to "//code.jquery.com/jquery-2.1.4.min.js")
* `PERSIST_CITATION_REPORTER_CHANGES` - Should changes through the web UI be persisted to disk (defaults to "True", other values considered False)
* `CITATION_REPORTER_PORT` - port to run on (defaults to $PORT if present or 8080 otherwise)

## Further Info

There are further instructions on installation and Sanger specific deployment considerations in [sanger_deployment](sanger_deployment/README.md).

## Known Issues

* Can't update users from through the web UI
* Doesn't support pseudonyms (e.g. change of surname; Bill vs William)
* You could use loads of memory if:
  * someone adds loads of publications
  * someone makes updates via the web UI quicker than they can be written to disk
* There is no authentication, you have to trust that users won't mess things up horribly (best to backup `publications.yml` periodically)
* Users don't have start and end dates which makes removing them a bit un-intuative
