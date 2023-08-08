import logging
import multiprocessing
import os
import subprocess
from subprocess import Popen
import time
import datetime

# region run external programs
def execute_command(input):
    """Execute command and measure execution time

    Args:
        cmd (list of str): Execution arguments
    """
    cmd=input[0]
    log_file=input[1]
    
    command = " ".join(cmd)
    logging.debug(f"Running command: {command}")
    start = time.perf_counter()

    log_file = open(log_file, "w") if log_file else subprocess.DEVNULL
    subprocess.run(cmd, check=True, stderr=log_file, stdout=log_file)

    runtime = time.perf_counter() - start
    runtime_str = datetime.timedelta(seconds=runtime)

    logging.debug(f"Finished in: {runtime_str}.")
    
def execute_command_parallel(inputs):
    """Execute command and measure execution time                                                                                      
    Args:
        cmd (list of str): Execution arguments                                                                                         
    """
    
    for i in inputs:
        cmd=i[0]
        log_file=i[1]
        try:
            
            command = " ".join(cmd)
            logging.debug(f"Running command: {command}")
            start = time.perf_counter()

            log_file = open(log_file, "w") if log_file else subprocess.DEVNULL
            procs = [ Popen(cmd, stderr=log_file, stdout=log_file) ]
            runtime = time.perf_counter() - start
            runtime_str = datetime.timedelta(seconds=runtime)

            logging.debug(f"Finished in: {runtime_str}.")


        except subprocess.CalledProcessError:
            raise Exception(
                f"Error occured while running GREMLIN for sequence {unp_id}."
            )
        
    for p in procs:
        p.wait()
    
    return procs
