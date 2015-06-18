from setuptools import setup
import multiprocessing

setup(name='citation_reporter',
      version='0.0.0',
      scripts=['scripts/web.py'],
      test_suite='nose.collector',
      tests_require=[
        'nose',
        'mock'
      ]
)
