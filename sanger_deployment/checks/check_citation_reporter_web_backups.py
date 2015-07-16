#!/usr/bin/env python

import os
import re
import sys
import datetime

from collections import OrderedDict

now = datetime.datetime.utcnow()

def delta_seconds(filepath):
  atime = datetime.datetime.utcfromtimestamp(os.path.getatime(filepath))
  delta = now - atime
  return float(delta.total_seconds())

def delta_hours(filepath):
  return delta_seconds(filepath) / (60*60)

def delta_days(filepath):
  return delta_hours(filepath) / 24

def delta_months(filepath):
  return delta_days(filepath) / (365 / 12)

def check_latest(filepath):
  return delta_hours(filepath) < 2

def check_daily(filepath):
  return delta_days(filepath) < 2

def check_weekly(filepath):
  return delta_days(filepath) < 10

def check_monthly(filepath):
  return delta_months(filepath) < 2

def newest(filepath):
  return delta_seconds(filepath)

backup_groups = OrderedDict([
  ('latest', []),
  ('daily', []),
  ('weekly', []),
  ('monthly', []),
  ('web', []),
  ('monit', [])
])

def group_by_filename(filepath):
  if os.path.getsize(filepath) == 0:
    return # file is empty, not a backup
  filename = os.path.basename(filepath)
  if re.findall('latest\.\d+\.yml.tgz$', filename):
    backup_groups['latest'].append(filepath)
  elif re.findall('daily-\d+\.yml.tgz$', filename):
    backup_groups['daily'].append(filepath)
  elif re.findall('weekly-\d+\.yml.tgz$', filename):
    backup_groups['weekly'].append(filepath)
  elif re.findall('monthly-\d+\.yml.tgz$', filename):
    backup_groups['monthly'].append(filepath)
  elif re.findall('web-\d+\.log$', filename):
    backup_groups['web'].append(filepath)
  elif re.findall('monit-\d+\.log$', filename):
    backup_groups['monit'].append(filepath)

parent_dir = os.path.abspath(sys.argv[1].strip())
for base_dir, dirs, files in os.walk(parent_dir):
  for filename in files:
    filepath = os.path.join(base_dir, filename)
    group_by_filename(filepath)

ok_backup_periods = []
failed_backup_periods = []
check_functions = {
  'latest': check_latest,
  'daily': check_daily,
  'weekly': check_weekly,
  'monthly': check_monthly,
  'web': check_daily,
  'monit': check_daily
}

for backup_period, backup_files in backup_groups.items():
  if len(backup_files) == 0:
    failed_backup_periods.append(backup_period)
    continue
  most_recent = sorted(backup_files, key=newest)[0]
  check_function = check_functions[backup_period]
  if check_function(most_recent):
    ok_backup_periods.append(backup_period)
  else:
    failed_backup_periods.append(backup_period)

ok = ", ".join(ok_backup_periods) if ok_backup_periods else "None"
fail = ", ".join(failed_backup_periods) if failed_backup_periods else "None"
print "OK => '%s', FAILED => '%s'" % (ok, fail)
if len(failed_backup_periods) == 0:
  exit(0)
else:
  exit(1)
