#!/share/apps/EPD_64bit/bin/python2.7
__author__ = 'Ryan Jacobs'
__version__ = '2.0'
__date__ = 'last updated April 14th, 2017'

import os
import sys
import logging
from VASP_Analyzer import JobMonitor, DirectoryUtilities

def main(resubmit_incomplete_jobs=False, submit_queued_jobs=False, resubmit_crashed_jobs=False):
    # Create log file
    logging.basicConfig(filename='job_monitor.log', level=logging.INFO)

    # Get current working directory
    cwd = os.getcwd()

    # Get all directories in the current working directory
    directory_list = DirectoryUtilities().get_downmost_directory_list()

    # Create an instance of JobMonitor
    jobmonitor = JobMonitor()

    # Parse the full directory list and only include directories containing VASP input files
    parsed_directory_list = jobmonitor._parse_directory_list(directory_list_to_parse=directory_list)

    running_job_dirs, queued_job_dirs, pruned_directory_list = jobmonitor._get_running_and_queued_jobs(directory_list=parsed_directory_list)
    completed_job_dirs, incomplete_job_dirs, pruned_directory_list2 = jobmonitor._get_complete_and_incomplete_jobs(directory_list=pruned_directory_list)
    crashed_job_dirs = jobmonitor._get_crashed_jobs(directory_list=parsed_directory_list)
    nonstarted_job_dirs, old_job_dirs, current_job_dirs = jobmonitor._get_nonstarted_and_old_jobs(directory_list=parsed_directory_list)

    if resubmit_incomplete_jobs == bool(True):
        resubmitted_jobs = jobmonitor._resubmit_incomplete_jobs(incomplete_job_dirs=incomplete_job_dirs, old_job_dirs=old_job_dirs)
        logging.info('Finished resubmitting incomplete jobs')
    elif resubmit_incomplete_jobs == bool(False):
        logging.info('You have chosen not to resubmit incomplete jobs')
        resubmitted_jobs = []
    if submit_queued_jobs == bool(True):
        submitted_new_jobs = jobmonitor._submit_nonstarted_jobs(nonstarted_job_dirs=nonstarted_job_dirs, old_job_dirs=old_job_dirs)
        logging.info('Finished submitting new jobs')
    elif submit_queued_jobs == bool(False):
        logging.info('You have chosen not to submit jobs that have not been started yet')
        submitted_new_jobs = []
    if resubmit_crashed_jobs == bool(True):
        resubmitted_crashed_jobs = jobmonitor._resubmit_crashed_jobs(crashed_job_dirs=crashed_job_dirs, old_job_dirs=old_job_dirs)
    elif resubmit_crashed_jobs == bool(False):
        logging.info('You have chosen not resubmit jobs that previously crashed')
        resubmitted_crashed_jobs = []

    jobmonitor._write_job_status_report(parent_directory=cwd, crashed_job_dirs=crashed_job_dirs,
                                        running_job_dirs=running_job_dirs, incomplete_job_dirs=incomplete_job_dirs,
                                        completed_job_dirs=completed_job_dirs, resubmitted_job_dirs=resubmitted_jobs,
                                        submitted_new_job_dirs=submitted_new_jobs, nonstarted_job_dirs=nonstarted_job_dirs,
                                        old_job_dirs=old_job_dirs)
    logging.info('Your job status report has been completed!')

if __name__ == "__main__":
    input = []
    try:
        if len(sys.argv[0]) > 0:
            input.append(sys.argv[0])
        elif len(sys.argv[0]) == 0:
            input.append(bool(False))
        if len(sys.argv[1]) > 0:
            input.append(sys.argv[1])
        elif len(sys.argv[1]) == 0:
            input.append(bool(False))
        if len(sys.argv[2]) > 0:
            input.append(sys.argv[2])
        elif len(sys.argv[2]) == 0:
            input.append(bool(False))
        main(resubmit_incomplete_jobs=bool(input[0]), submit_queued_jobs=bool(input[1]), resubmit_crashed_jobs=bool(input[2]))
    except(IndexError):
        main(resubmit_incomplete_jobs=bool(False), submit_queued_jobs=bool(False), resubmit_crashed_jobs=bool(False))