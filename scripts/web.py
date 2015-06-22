#!/usr/bin/env python

import logging
import os

from citation_reporter.PubmedSearch import Publication, Publications
from citation_reporter.Author import Author

import flask
from flask import Flask, render_template, request, jsonify, make_response
from StringIO import StringIO

parent_folder = os.path.abspath(os.path.dirname(__file__))
template_folder = os.path.join(parent_folder, '..', 'templates')
static_folder = os.path.join(parent_folder, '..', 'static')
publication_data_filename = os.path.join(parent_folder, '..', 'publications.yml')

app = Flask(__name__, template_folder=template_folder,
            static_folder=static_folder)
logging.basicConfig(level=logging.DEBUG)

@app.route('/')
def publications():
  return render_template('index.html',
                        publications=publications.not_denied(),
                        download_link='/publications.csv')

@app.route('/publications.csv')
def download():
  output_file = StringIO()
  output_publications = publications.not_denied()
  output_csv = output_publications.to_csv(output_file)
  output = make_response(output_file.getvalue())
  output.headers["Content-Disposition"] = "attachment; filename=publications.csv"
  output.headers["Content-type"] = "text/csv"
  return output

@app.route('/trash')
def trash():
  return render_template('index.html',
                        publications=publications.denied(),
                        download_link='/trash.csv')

@app.route('/trash.csv')
def download_trash():
  output_file = StringIO()
  output_publications = publications.denied()
  output_csv = output_publications.to_csv(output_file)
  output = make_response(output_file.getvalue())
  output.headers["Content-Disposition"] = "attachment; filename=trash_publications.csv"
  output.headers["Content-type"] = "text/csv"
  return output

@app.route('/publication/<pubmed_id>/<author_string>/', methods=['GET', 'PUT'])
def update_publication_authors(pubmed_id, author_string):
  def refresh_publication_page():
    return render_template('publication.html',
                           publication=publication,
                           pubmed_id=pubmed_id)
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
    return refresh_publication_page()

@app.route('/publication/<pubmed_id>/', methods=['DELETE'])
def delete_publication(pubmed_id):
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
    publications = {}
  app.run(host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))
