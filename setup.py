from setuptools import setup
import multiprocessing

setup(name='citation_reporter',
      version='0.0.0',
      scripts=['scripts/citation_reporter.py', 'scripts/citation_reporter_web.py'],
      test_suite='nose.collector',
      tests_require=[
        'nose',
        'mock'
      ]
)
