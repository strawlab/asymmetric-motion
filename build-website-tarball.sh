#!/bin/bash
set -o errexit

docker build -t strawlab/asymmetric-motion .

docker run \
        --name asymmetric-motion-build-report \
        strawlab/asymmetric-motion

TEMPDIR=`mktemp -d`
FNAME="asymmetric-motion-website-files.tar.gz"
docker cp asymmetric-motion-build-report:/asymmetric-motion/tutorial-src ${TEMPDIR}
mv ${TEMPDIR}/tutorial-src/${FNAME} .
docker rm asymmetric-motion-build-report
rm -rf ${TEMPDIR}
echo "your file awaits you at ${FNAME}"
