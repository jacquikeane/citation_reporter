#!/usr/bin/env python

import os

from citation_reporter.PubmedSearch import Publication

import flask
from flask import Flask, render_template

parent_folder = os.path.abspath(os.path.dirname(__file__))
template_folder = os.path.join(parent_folder, '..', 'static', 'templates')
publication_data_filename = os.path.join(parent_folder, '..', 'publications.yml') 

app = Flask(__name__, template_folder=template_folder)

@app.route('/')
def index():
  return render_template('index.html',
                        publications=publications)

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
