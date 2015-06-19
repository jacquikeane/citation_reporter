#!/usr/bin/env python

import logging
import os

from citation_reporter.PubmedSearch import Publication

import flask
from flask import Flask, render_template, request, jsonify

parent_folder = os.path.abspath(os.path.dirname(__file__))
template_folder = os.path.join(parent_folder, '..', 'static', 'templates')
publication_data_filename = os.path.join(parent_folder, '..', 'publications.yml') 

app = Flask(__name__, template_folder=template_folder)
logging.basicConfig(level=logging.DEBUG)

@app.route('/')
def index():
  return render_template('index.html',
                        publications=publications)

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

if __name__ == '__main__':
  try:
    with open(publication_data_filename, 'r') as publication_data_file:
      publications = Publication.from_yaml(publication_data_file.read())
  except IOError:
    logging.warning("Could not open %s to read publications" %
                    publication_data_filename)
    publications = {}
  except:
    logging.error("There was a problem parsing publications from %s" %
                  publication_data_filename)
    publications = {}
  app.run(host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))
