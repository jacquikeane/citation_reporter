#!/usr/bin/env python

import datetime
import logging
import os
import pkg_resources
import re
import time

from citation_reporter.LibrarySearcher import LibrarySearcher
from citation_reporter.PubmedSearch import Publication, Publications
from citation_reporter.Author import Author, User

import flask
from boltons.strutils import slugify
from collections import OrderedDict
from flask import Flask, render_template, request, jsonify, \
                         make_response, url_for, redirect, \
                         flash, abort
from functools import wraps
from StringIO import StringIO
from threading import Thread
from Queue import Queue

logging.basicConfig(level=logging.DEBUG)

parent_folder = os.path.abspath(os.path.dirname(__file__))
template_folder = os.path.join(parent_folder, '..', 'templates')
if not os.path.isdir(template_folder):
  template_folder = pkg_resources.resource_filename('citation_reporter',
                                                    'templates')
logging.debug("Loading templates from %s" % template_folder)

app = Flask(__name__, template_folder=template_folder)

def parse_pubmed_ids(publications_string):
  publication_ids = publications_string.replace(',', ' ').replace(';', ' ').split(' ')
  publication_ids = [pid for pid in publication_ids if re.match('^\d{8}$', pid)]
  return publication_ids

def new_publications_stats(new_publications):
  new_trash_publications_set = set(new_publications.denied().keys())
  new_publications_set = set(new_publications.keys())
  existing_publication_set = set(publications.keys())
  new_publications_count = len(new_publications_set.difference(existing_publication_set))
  new_trash_publications_count = len(new_trash_publications_set.difference(existing_publication_set))
  duplicate_publications_count = len(new_publications_set.intersection(existing_publication_set))
  return (new_publications_count,
          new_trash_publications_count,
          duplicate_publications_count)

def message_about_new_publications(new_publications):
  new_count, trash_count, duplicate_count = new_publications_stats(new_publications)
  if new_count == 1:
    message = "Added 1 new publication"
  else:
    message = "Added %s new publications" % new_count
  if trash_count > 0:
    message += "; %s had no potential authors found" % trash_count
  if duplicate_count > 0:
    message += "; %s duplicates of an existing publication" % duplicate_count
  logging.debug(message)
  return message

def save_publications(data_to_write, publication_data_filename):
  while(True):
    for i in range(20):
      if data_to_write.qsize() > 1:
        data_to_write.get()
      else:
        break
    timestamp, latest_publications = data_to_write.get()
    data_to_write.publications = None
    # We want an "atomic" write so write to a temporary file and then move it to
    # the file we actually want to persist to.
    with open(publication_data_filename + ".part", 'w') as publications_file:
      publications_file.write(latest_publications.sorted_by_date().to_yaml())
      publications_file.flush()
      os.fsync(publications_file.fileno())
    os.rename(publication_data_filename + ".part", publication_data_filename)
    logging.debug("Have saved publications from '%s' to disk.  Queue is %s long" %
                  (timestamp.isoformat(), data_to_write.qsize()))

def save_changes(func):
  @wraps(func)
  def decorator(*args, **kwargs):
    result = func(*args, **kwargs)
    if app.config['PERSIST_CHANGES'] and not request.method in ["GET", "HEAD"]:
      global data_to_write
      data_to_write.put((datetime.datetime.now(), publications))
      logging.debug("Queued publications to be saved to disk. Queue is %s long" % data_to_write.qsize())
    elif request.method in ["GET", "HEAD"]:
      logging.debug("%s request, no changes to publications" % request.method)
    else:
      logging.warning("Persistance disabled, have not saved to disk")
    return result
  return decorator

@app.route('/', methods=["GET", "POST"])
@save_changes
def publications_page():
  global publications
  if request.method == 'POST':
    publication_ids = parse_pubmed_ids(request.form['pubmed_ids'])
    logging.debug("Request to add publications: %s" % "; ".join(map(str, publication_ids)))
    new_publications = Publications.from_pubmed_ids(publication_ids)
    for publication in new_publications.values():
      publication.update_authors(users)
      logging.debug("Found %s potential authors for %s" %
                    (len(publication.most_likely_affiliated_authors()),
                     publication.pubmed_id))
    flash(message_about_new_publications(new_publications))
    publications = Publications.merge(publications, new_publications)
    return redirect(url_for('publications_page'))
  return render_template('affiliated.html',
                        publications=publications.not_denied(),
                        users=users,
                        user_title="All affiliated publications",
                        download_link=url_for('download'))

@app.route('/user/<user_id>/')
def user_page(user_id):
  if not user_id in users:
    return render_template('affiliated.html',
                          publications=Publications(),
                          users=users,
                          user_title="Unknown user",
                          download_link='#'), 404
  return render_template('affiliated.html',
                        publications=publications.filter_by_user_id(user_id),
                        users=users,
                        user_title="Publications for %s" % users[user_id].full_name(),
                        download_link=url_for('download', user_id=user_id))

@app.route('/publications.csv')
def download():
  output_file = StringIO()
  user_id = request.args.get("user_id")
  if user_id:
    if user_id not in users:
      abort(404)
    output_publications = publications.filter_by_user_id(user_id,
                                                         include_denied=False)
  else:
    output_publications = publications.not_denied()
  output_csv = output_publications.to_csv(output_file)
  output = make_response(output_file.getvalue())
  if user_id:
    user_name = slugify(users[user_id].full_name())
    output.headers["Content-Disposition"] = "attachment; filename=%s_publications.csv" % user_name
  else:
    output.headers["Content-Disposition"] = "attachment; filename=publications.csv"
  output.headers["Content-type"] = "text/csv"
  return output

@app.route('/trash/')
def trash_page():
  return render_template('trash.html',
                        publications=publications.denied().has_potential_author(),
                        user_title="Trash",
                        users=users,
                        download_link=url_for('download_trash'))

@app.route('/trash.csv')
def download_trash():
  output_file = StringIO()
  output_publications = publications.denied().has_potential_author()
  output_csv = output_publications.to_csv(output_file)
  output = make_response(output_file.getvalue())
  output.headers["Content-Disposition"] = "attachment; filename=trash_publications.csv"
  output.headers["Content-type"] = "text/csv"
  return output

@app.route('/unaffiliated/')
def unaffiliated():
  return render_template('trash.html',
                        publications=publications.denied().has_no_potential_author(),
                        users=users,
                        user_title="Unaffiliated publications",
                        download_link=url_for('download_unaffiliated'))

@app.route('/unaffiliated.csv')
def download_unaffiliated():
  output_file = StringIO()
  output_publications = publications.denied().has_no_potential_author()
  output_csv = output_publications.to_csv(output_file)
  output = make_response(output_file.getvalue())
  output.headers["Content-Disposition"] = "attachment; filename=unaffiliated_publications.csv"
  output.headers["Content-type"] = "text/csv"
  return output

def refresh_publication_page(publication):
  if publication.has_affiliated_authors():
    isTrash=False
  else:
    isTrash=True
  return render_template('publication.html',
                         publication=publication,
                         isTrash=isTrash,
                         pubmed_id=publication.pubmed_id)

@app.route('/publication/<pubmed_id>/<author_string>/', methods=['GET', 'PUT'])
@save_changes
def update_publication_authors(pubmed_id, author_string):
  global publications
  publication = publications.get(pubmed_id)
  if publication == None:
    error = {"error": "Couldn't find publication with pubmed_id: %s" %
            pubmed_id}
    return jsonify(**error), 404
  elif not author_string in publication.affiliated_authors.keys():
    error = {"error": "Couldn't find potential author '%s' for pubmed_id: %s" %
            (author_string, pubmed_id)}
    return jsonify(**error), 404
  else:
    for user_id, status in request.get_json().items():
      try:
        publication.update_author_status(author_string, user_id, status)
      except KeyError:
        logging.warning("Could not set status '%s' for user '%s' on publication '%s'; skipping" %
                        (status, user_id, pubmed_id))
        continue
    return refresh_publication_page(publication)

@app.route('/publication/<pubmed_id>/', methods=['GET', 'DELETE'])
@save_changes
def publication(pubmed_id):
  global publications
  if request.method == 'GET':
    publication = publications[pubmed_id]
    return refresh_publication_page(publication)
  try:
    publications[pubmed_id].deny()
  except KeyError:
    error = {"error": "Could not find publication with id '%s'" % pubmed_id}
    return jsonify(**error), 404
  except:
    raise
  else:
    message = {"success": "Moved publication '%s' to trash" % pubmed_id}
    return jsonify(**message)

@app.route('/missing')
def missing():
  """Hidden URL for missing publications

  I've intentionally not publicised this feature because it makes lots of HTTP
  requests to the Sanger website and I probably only want me / administrators to
  use it on an infrequent basis.

  This endpoint lists all of the publicatons which are currently not in the
  Sanger Library's dataset"""
  library_pubmed_ids = set(LibrarySearcher.get_pubmed_ids())
  possible_pubmed_ids = set(publications.not_denied().keys())
  missing_pubmed_ids = possible_pubmed_ids.difference(library_pubmed_ids)
  missing_publications = Publications([publications[pubmed_id] for pubmed_id in
                                       missing_pubmed_ids])
  return render_template('affiliated.html',
                        publications=missing_publications,
                        users=users,
                        user_title="Publications not in the Sanger Library")

@app.route('/stats')
def stats():
  """Hidden URL for statistics

  I've intentionally not publicised this feature because it makes lots of HTTP
  requests to the Sanger website and I probably only want me / administrators to
  use it on an infrequent basis"""
  all_pubmed_ids = set(publications.keys())
  confirmed_pubmed_ids = set(publications.confirmed().keys())
  denied_pubmed_ids = set(publications.denied().keys())
  to_be_confirmed_ids = all_pubmed_ids.difference(confirmed_pubmed_ids, denied_pubmed_ids)

  library_pubmed_ids = set(LibrarySearcher.get_pubmed_ids())
  confirmed_in_library_ids = library_pubmed_ids.intersection(confirmed_pubmed_ids)
  to_be_confirmed_in_library_ids = library_pubmed_ids.intersection(to_be_confirmed_ids)

  not_in_libary_ids = all_pubmed_ids.difference(library_pubmed_ids)
  confirmed_not_in_library = not_in_libary_ids.intersection(confirmed_pubmed_ids)
  to_be_confirmed_not_in_library = not_in_libary_ids.intersection(to_be_confirmed_ids)

  statistics = {
    'tracked: total': len(all_pubmed_ids),
    'tracked: confirmed': len(confirmed_pubmed_ids),
    'tracked: possible': len(to_be_confirmed_ids),
    'tracked: denied': len(denied_pubmed_ids),
    'library: total': len(library_pubmed_ids),
    'library: confirmed': len(confirmed_in_library_ids),
    'library: possible': len(to_be_confirmed_in_library_ids),
    'other: total': len(confirmed_not_in_library) + len(to_be_confirmed_not_in_library),
    'other: confirmed': len(confirmed_not_in_library),
    'other: possible': len(to_be_confirmed_not_in_library)
  }

  return jsonify(**statistics)

def update_config():
  global app
  default_publication_data_filename = os.path.join(os.getcwd(), 'publications.yml')
  app.config['PUBLICATIONS_FILE'] = os.environ.get('CITATION_REPORTER_PUBLICATIONS',
                                                    default_publication_data_filename)

  default_user_data_filename = os.path.join(os.getcwd(), 'authors.yml')
  app.config['USERS_FILE'] = os.environ.get('CITATION_REPORTER_USERS',
                                                    default_user_data_filename)
  app.config['BOOTSTRAP_URL'] = os.environ.get('BOOTSTRAP_URL',
                                                    "//maxcdn.bootstrapcdn.com/bootstrap/3.3.5")
  app.config['JQUERY_URL'] = os.environ.get('JQUERY_URL',
                                             "//code.jquery.com/jquery-2.1.4.min.js")
  app.config['PERSIST_CHANGES'] = os.environ.get('PERSIST_CITATION_REPORTER_CHANGES',
                                                    "True") == "True"
  app.config['PORT'] = int(os.environ.get('CITATION_REPORTER_PORT',
                                           os.environ.get('PORT',
                                                           8080)))

def load_publications():
  try:
    with open(app.config['PUBLICATIONS_FILE'], 'r') as publication_data_file:
      publications = Publications.from_yaml(publication_data_file.read())
    logging.info("Loaded publications from %s" %
                 app.config['PUBLICATIONS_FILE'])
  except IOError:
    logging.warning("Could not open %s to read publications" %
                    app.config['PUBLICATIONS_FILE'])
    publications = Publications([])
  except Exception as e:
    logging.error("There was a problem parsing publications from %s" %
                  app.config['PUBLICATIONS_FILE'], exc_info=e)
    publications = Publications([])
  return publications

def load_users():
  try:
    with open(app.config['USERS_FILE'], 'r') as user_data_file:
      users = User.from_yaml(user_data_file.read())
    logging.info("Loaded %s users from %s" % (len(users),
                                              app.config['USERS_FILE']))
  except IOError:
    logging.warning("Could not open %s to read users" %
                    app.config['USERS_FILE'])
    users = {}
  except Exception as e:
    logging.error("There was a problem parsing users from %s" %
                      app.config['USERS_FILE'], exc_info=e)
    users = {}

  publication_users = publications.get_users()
  publication_users.update(users)
  users = publication_users

  def by_name(user_id__user):
    user_id, user = user_id__user
    return user.ordered_name()
  users = OrderedDict(sorted(users.items(), key=by_name))

  logging.info("Merge users from userfile with users from publications, %s found" % len(users))

  return users

if __name__ == '__main__':
  update_config()
  publications = load_publications()
  users = load_users()

  for publication in publications.values():
    publication.update_authors(users)
  logging.info("Updated publication authors")

  data_to_write = Queue()
  save_publications_process = Thread(target=save_publications,
                                      args=(data_to_write,
                                            app.config['PUBLICATIONS_FILE']))
  save_publications_process.daemon = True
  save_publications_process.start()

  app.secret_key = os.environ.get("FLASK_SECRET_KEY", os.urandom(24))
  app.run(host="0.0.0.0", port=app.config['PORT'])
