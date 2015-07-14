# Mostly Sanger Specific Installation Notes

The following notes and scripts are here for the convenience of colleagues who may need to install or maintain this at the Wellcome Trust Sanger Institute.  I've included this here because some of these scripts and notes may be a handy starting point for other teams.

## Installation

### Dependencies

You will need the following software:
```
gcc
git
python
python-dev
logrotate
```

You probably also want to have `python-virtualenv` installed for sanity

### System requirements

The server will need to have some non-priviledged ports (ports > 1024) open internally but these should not be accessibly from the internet.  Pick one or use the default which is 8080.

The server needs between 500 and 800 MBs to hold all of Sanger's publications since 2000.  This is probably an overestimate of the resource requirements.

### Installing citation_reporter

Create a folder and copy in an `authors.yml` file (and `publications.yml` file if available).  Create a new virtualenv into which you will install citation_reporter, activate it and install.  For example:

```
virtualenv venv -p $(which python2)
. venv/bin/activate
pip install git+https://github.com/sanger-pathogens/citation_reporter.git@master
```

In this case I have installed directly from the master branch on Github.  Releases may also be available from PyPi or you might want to use another branch (e.g. `master`).  NB, if you've already installed citation_reporter, you will need to append `--upgrade` to the `pip` command else no changes will be made.

### Setup Backups

* Copy [backup.logrotate.conf](backup.logrotate.conf) onto the server;
* Amend backup.logrotate.conf to point to the correct directories;
* Ensure that a [backups](backups) folder exists and create associated files if they don't exist (empty files are fine);
* Amend [backup.logrotate.crontab](backup.logrotate.crontab) to point to the correct config; and
* Add / edit the cronjob `crontab -e`

### Install monit

* Download monit from the [monit website](https://mmonit.com/monit/#download)
* Configure it so that it is installed into a non-priviledged location (something like `./configure --prefix=my_folder`)
* Make and install it according to the latest instructions in the downloaded README
* Edit citation_reporter_web 'init' script (e.g. paths) and copy it onto the server
* Edit [monitrc](monitrc) and copy it onto the server
  * Correct paths
  * Enable / disable the monit web UI or change the port
  * Update the alert email details so that you're warned rather than me :)
* Copy the [check scripts](checks) onto the server and update relevant paths in the scripts
* Edit the paths in [monit.crontab](monit.crontab) and update cron (`crontab -e`)
* [Optional] Add the monit/bin directory to your `PATH` and / or `PATH` in your .bash_profile`
* Start monit - `monit -c path_to_monitrc -l monit.log`
* Check the status of services - `monit -c path_to_monitrc status`

The citation_reporter_web daemon should now be running and should restart pretty quickly if you were to `kill` the relevant process.

## Management

### Starting / Stopping citation_reporter_web

* Start citation_reporter - `monit -c monitrc_path_here start citation_reporter_web`
* Stop citation_reporter - `monit -c monitrc_path_here stop citation_reporter_web`

### Updating citation_reporter

* Check that no-one is using citaton_reporter (e.g. `tail -f web.log`)
* Stop citation_reporter - `monit -c monitrc_path_here stop citation_reporter_web`
* [Optional] Backup the virtualenv, `publications.yml` and `authors.yml`
* Source the virtualenv - `. path_to_venv/bin/activate`
* Update citation_reporter - `pip install git+https://github.com/sanger-pathogens/citation_reporter.git@YOUR_BRANCH_HERE --upgrade`
* Deactivate the virtualenv - `deactivate`
* Restart the web application - `monit -c monitrc_path_here start citation_reporter_web`

### Update publications

* If you know the publications' pubmed ids, the easiest solution is to copy and paste a list of them into the web UI.  These can be space or comma deliminated.
* If you want to search pubmed, make sure that the `authors.yml` list is up to date, stop the web server, run `./citation_reporter_web update-publications` and restart the web server.  NB this script would stop and start the server itself but monit would get confused so it is important to stop and start the server with monit.

### Editing authors

* To add authors, edit `authors.yml` and then follow the instructions from 'Update publications' using the cli
* Removing authors is a bit harder.  If you remove them from the `authors.yml` then the cli script will not automatically add new publications using them.  However if you add a publication through the web UI it tries to match the authors against the authors in the `authors.yml` file and which have previously been matched against a publication.  I personally think this is the least bad behaviour but it would be solved by giving authors start and end dates.
* If you want to completely remove an author you will need to stop the web server, delete the relevant YAML from `publications.yml` and `authors.yml` and restart the server.

### Viewing the monit checks in a web browser

Monit has a handy web view to show you how the web server is behaving.  This can only be accessed from the web server itself so you will need to setup an SSH tunnel.  The default credentials are "admin/monit" but these can be changed in [monitrc](monitrc).

You can setup the SSH tunnel with the following command:

```
ssh ssh -L 28120:localhost:2812 <THE_SERVER_YOU'RE_HOSTED_ON>
```

and then visit `localhost:28120` in your browser

### Adding additional checks

The [monit documentation](https://mmonit.com/monit/documentation/monit.html) is pretty good.  My advice is to write a script which exits with a non-zero status and prints a short failure message and put it into the [checks](checks) directory.  The edit [monitrc](monitrc) to execute that check (by copying one of the existing examples).

### Changing alert recipients

Edit [monitrc](monitrc)

### Testing in a Dev environment

* Create a small list of authors
* Install citation_reporter into a new directory with a new virtualenv (`deactivate` the production one first)
* Run `citation_reporter_cli.py -u authors.yml -s 2015 -e 2015 -v` to find some publications
* It would be a good idea to check that the server you are running on has 300MB+ of RAM before starting and that you don't have too many publications.
* Run `CITATION_REPORTER_PORT=8081 PERSIST_CITATION_REPORTER_CHANGES=False citation_reporter_web.py`
* Visit the dev version on port '8081'
* Turn off the dev environment

### Fixing errors in publications

* Stop the server
* Backup `publications.yml`
* Diff against on of the 'latest' [backups](backups)
* Make some changes
* Maybe check the changes in a dev environment
* Start the server

### Removing publications

In most cases it is fine just to click on the 'trash' icon in the web browser.  

It there is a publication that you never want anyone to see through the web UI you can also manually edit the `publications.yml` file to remove all references to authors and set the confirmation status to 'denied'.  This will mean that the publication will never be shown even if new authors are added in the future or if someone tries to manually add the publication through the web UI.
