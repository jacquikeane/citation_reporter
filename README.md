# Citation Reporter
Tool to keep track of a team's publications. 

[![Build Status](https://travis-ci.org/sanger-pathogens/citation_reporter.svg?branch=master)](https://travis-ci.org/sanger-pathogens/citation_reporter)   
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/citation_reporter/blob/master/LICENSE)   

## Contents
  * [Introduction](#introduction)
  * [Installation](#installation)
    * [Required dependencies](#required-dependencies)
    * [Using pip](#using-pip)
    * [Running the tests](#running-the-tests)
  * [Usage](#usage)
    * [scripts/citation\_reporter\_cli\.py](#scriptscitation_reporter_clipy)
    * [scripts/citation\_reporter\_web\.py](#scriptscitation_reporter_webpy)
  * [License](#license)
  * [Feedback/Issues](#feedbackissues)
  * [Further Information](#further-information)
  * [Known Issues](#known-issues)

## Introduction
Citation reporter is a tool to keep track of a team's publications. It currently only searches pubmed and the Sanger Library for publications and matches these with a list of users in your group. Publications by authors from the team are refered to as 'affiliated' in the UI and the code.

The tool is based on scripts developed by @simonrharris.

## Installation
citation_reporter has the following dependencies:

### Required dependencies
* biopython==1.70
* boltons==0.6.4
* Flask==0.12.3
* PyYAML==4.2b1
* requests==2.20.0

Details for installing citation_reporter are provided below. If you encounter an issue when installing citation_reporter please contact your local system administrator. If you encounter a bug please log it [here](https://github.com/sanger-pathogens/citation_reporter/issues) or email us at path-help@sanger.ac.uk.

### Using pip

```
pip install git+https://github.com/sanger-pathogens/citation_reporter.git@master
```

### Running the tests
The test can be run from the top level directory:  
```
python setup.py test
```

## Usage
Publications can be loaded through the web UI or a commandline script.  New users can only be added via the commandline script.

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

## License
citation_reporter is free software, licensed under [GPLv3](https://github.com/sanger-pathogens/citation_reporter/blob/master/LICENSE).

## Feedback/Issues
Please report any issues to the [issues page](https://github.com/sanger-pathogens/citation_reporter/issues) or email path-help@sanger.ac.uk.

## Further Information
There are further instructions on installation and Sanger specific deployment considerations in [sanger_deployment](sanger_deployment/README.md).

## Known Issues
* Can't update users from through the web UI
* Doesn't support pseudonyms (e.g. change of surname; Bill vs William)
* You could use loads of memory if:
  * someone adds loads of publications
  * someone makes updates via the web UI quicker than they can be written to disk
* There is no authentication, you have to trust that users won't mess things up horribly (best to backup `publications.yml` periodically)
* Users don't have start and end dates which makes removing them a bit un-intuative
