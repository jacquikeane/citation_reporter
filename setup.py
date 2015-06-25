import os
from setuptools import setup, find_packages
import multiprocessing

template_files = os.listdir(os.path.join(os.path.dirname(__file__), 'templates'))
templates = [os.path.join('templates', template) for template in template_files]

example_files = os.listdir(os.path.join(os.path.dirname(__file__), 'examples'))
examples = [os.path.join('examples', example) for example in example_files]

setup(name='citation_reporter',
      version='0.0.0',
      scripts=[
        'scripts/citation_reporter_cli.py', 
        'scripts/citation_reporter_web.py'
      ],
      test_suite='nose.collector',
      tests_require=[
        'nose',
        'mock'
      ],
      install_requires=[
        'biopython',
        'boltons',
        'Flask',
        'PyYAML',
        'requests'
      ],
      include_package_data=True,
      package_data={
        'templates': 'templates/*',
        'examples': 'examples/*'
      },
      packages=find_packages(),
      zip_safe=False
)
