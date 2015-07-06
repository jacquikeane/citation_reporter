#!/bin/bash

URL="localhost:8080?test=True"
RESPONSE=$(curl --head --silent $URL)
EXPECTED_RESPONSE='HTTP/1.0 200 OK'

if [ -z "$RESPONSE" ]; then
  echo "Empty response from server, server probably down"
else
  echo "$RESPONSE"
fi

echo "$RESPONSE" | grep "$EXPECTED_RESPONSE" > /dev/null
