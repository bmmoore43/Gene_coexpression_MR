#!/bin/bash
# Simple user epilogue file that echoes everything it learns
# about the job it was attached to.
# 'Job Account' is almost always going to be blank.

echo "\$1 - Job ID: $1"
echo "\$2 - User Name: $2"
echo "\$3 - Group Name: $3"
echo "\$4 - Job Name: $4"
echo "\$5 - Session ID: $5"
echo "\$6 - Resources Requested: $6"
echo "\$7 - Resources Used: $7"
echo "\$8 - Queue: $8"
echo "\$9 - Job Account: $9"
# Note - we can't use '$10' - we need to put positions > 9 in braces
echo "\$10 - Job exit status: ${10}"

exit 0

