#!/usr/bin/env python

import logging
import os
import re

from citation_reporter.PubmedSearch import Publication, Publications
from citation_reporter.Author import Author

import flask
from boltons.strutils import slugify
from flask import Flask, render_template, request, jsonify, \
                         make_response, url_for, redirect, \
                         flash, abort
from StringIO import StringIO

parent_folder = os.path.abspath(os.path.dirname(__file__))
template_folder = os.path.join(parent_folder, '..', 'templates')
static_folder = os.path.join(parent_folder, '..', 'static')
publication_data_filename = os.path.join(parent_folder, '..', 'publications.yml')

app = Flask(__name__, template_folder=template_folder,
            static_folder=static_folder)
logging.basicConfig(level=logging.DEBUG)

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

@app.route('/', methods=["GET", "POST"])
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
def update_publication_authors(pubmed_id, author_string):
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
def publication(pubmed_id):
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

if __name__ == '__main__':
  try:
    with open(publication_data_filename, 'r') as publication_data_file:
      publications = Publications.from_yaml(publication_data_file.read())
  except IOError:
    logging.warning("Could not open %s to read publications" %
                    publication_data_filename)
    publications = Publications([])
  except Exception as e:
    logging.error("There was a problem parsing publications from %s: %s" %
                  (publication_data_filename, e))
    publications = Publications([])

  users = publications.get_users()
  for publication in publications.values():
    publication.update_authors(users)

  app.secret_key = os.environ.get("FLASK_SECRET_KEY", os.urandom(24))
  app.run(host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))
