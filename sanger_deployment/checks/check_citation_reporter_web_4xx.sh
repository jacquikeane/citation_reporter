#!/bin/bash

WORKING_FOLDER=/var/opt/citation_reporter_prod
LOGFILE=${WORKING_FOLDER}/web.log

# Cap response at the first 100 failures so that we don't run out of memory if 
# it all goes horribly wrong
RELEVANT_LINES=$(/bin/egrep '^INFO:werkzeug.+4[0-9]{2} -$' $LOGFILE | head -100)
COUNT=$(echo "$RELEVANT_LINES" | wc -l)

SUMMARY=$(/bin/echo "$RELEVANT_LINES" | /usr/bin/awk -F'"' '{print $2,$3}' | /usr/bin/sort | /usr/bin/uniq -c | sed 's/^\s\+//' | sort -rn)

if [ ! -z "$RELEVANT_LINES" ]; then
  echo "$COUNT matching lines"
  echo "$SUMMARY" && /usr/bin/test $COUNT -lt 10
else
  echo "No matching lines"
fi
