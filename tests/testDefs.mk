# force sort to be consistent
export LC_ALL=C
.SECONDARY:

# remove {make out exists ... } from job file so it can be executed as commands
jobsToCmds = sed -Ee 's/\{check out exists (.+)\}/\1/'

# basic command to dump a database created by a test
sqldumpcmd = sqlite3 -header -batch output/$@.db

# call function with table as argument to dump table and diff with expected, ignoring the id column
sqldumpdiff = ${sqldumpcmd} 'select * from $(1)' | cut -f 2- > output/$@.$(1).tsv && ${diff} expected/$@.$(1).tsv output/$@.$(1).tsv

hostname = $(shell hostname)
