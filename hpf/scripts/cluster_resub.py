#!/usr/bin/env python

# Quick script to slowly re-submit bad tasks

from time import sleep
from subprocess import call

handle = open('redo01.tasks')
tasks = list()

for line in handle:
	tasks.append(line.rstrip())
handle.close()

for task in tasks:
	print "Submitting job {0}".format(task)
	call(["qsub", "-t",  task, "jobscript.pbs"])
	sleep(1)

print "Complete"

