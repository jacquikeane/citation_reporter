set -ex

WORKING_DIR="/var/opt/citation_reporter_prod"
VENV="${WORKING_DIR}/venv"
MONIT_CONFIG="${WORKING_DIR}/monitrc"
MONIT_BINARY="${WORKING_DIR}/monit/bin/monit"

. ${VENV}/bin/activate
$MONIT_BINARY -c $MONIT_CONFIG stop all
${WORKING_DIR}/citation_reporter_web update-publications
$MONIT_BINARY -c $MONIT_CONFIG start citation_reporter_web
sleep 60
$MONIT_BINARY -c $MONIT_CONFIG start all
