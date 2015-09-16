set -ex

WORKING_DIR="/var/opt/citation_reporter_prod"
VENV="${WORKING_DIR}/venv"
MONIT_CONFIG="${WORKING_DIR}/monitrc"

. ${VENV}/bin/activate
monit -c $MONIT_CONFIG stop all
${WORKING_DIR}/citation_reporter_web update-publications
monit -c $MONIT_CONFIG start citation_reporter_web
sleep 60
monit -c $MONIT_CONFIG start all
