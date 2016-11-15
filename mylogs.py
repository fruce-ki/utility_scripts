#!~/bin/python3

"""mylogs.py

Author: Kimon Froussios
Compatibility tested: python 2.7.3, 3.5.2
Last revised: 29/08/2016

Library for custom logging.

XXXstring() functions return formatted strings for custom use by the caller.
log_XXX() functions format and record messages to log-files.

NOTE: XXXstring() return values are not intended to be passed as arguments to 
log_XXX(), as this will cause duplication and will mess up the formatting. Use 
XXXstring() to print info to the screen or within your output files. Use 
log_XXX() to record to log-files.
"""

import sys, datetime, os, subprocess


# Do everything in ONE write operation, to prevent race conditions when accessing the logfile. 


def tstamp():
    """Readable current date-time.

    Returns:
        (str)
    """
    return str(datetime.datetime.now())


def paramstring(message=""):
    """Execution parameters log-string.
    
    I normally use this at the beginning of my output.
    
    Args:
        message(str): Optional message to add.
    Returns:
        (str)
    """
    message.rstrip().replace("\n"," ")
    return "###INFO### "+ tstamp() +"\t"+ " ".join(sys.argv) +"\n###INFO### \t\t\t\t### CWD: " + os.getcwd() + "  ### PYTHON: " + sys.executable + ("\n###INFO### " if message else "") + message + "\n"


def donestring(message=""):
    """Done log-string.
    
    I normally use this at the end of tasks.
    
    Args:
        message(str): A custom message to add.
    Returns:
        (str)
    """
    message.rstrip().replace("\n"," ")
    return "###INFO### " + tstamp() + "\t" + os.path.basename(sys.argv[0]) + " - " + "Done "+ message +".\n"


def infostring(message=""):
    """Info log-string.
    
    I normally use this at the end of tasks.
    
    Args:
        message(str): A custom message to add.
    Returns:
        (str)
    """
    message.rstrip().replace("\n"," ")
    return "###INFO### " + tstamp() + "\t" + os.path.basename(sys.argv[0]) + " - " + message +"\n"


def warnstring(message=""):
    """Warning log-string.
    
    Args:
        message(str): The message to add.
    Returns:
        (str)
    """
    message.rstrip().replace("\n"," ")
    return "#WARNING!# " + tstamp() + "\t" + os.path.basename(sys.argv[0]) + " - " + message +"\n"


def errstring(message=""):
    """Error log-string.
    
    Args:
        message(str): The message to add.
    Returns:
        (str)
    """
    message.rstrip().replace("\n"," ")
    return "##!ERROR!# " + tstamp() + "\t" + os.path.basename(sys.argv[0]) + " - " + message +"\n"


def log_command(message="", logfile = "./commands.log" ):
    """Record timestamp, command-line call and optional message.
    
    This function obtains the command from sys.argv[].
    
    Args:
        message(str): Message to record in the log. (Default empty)
        logfile(str): File to write to (./commands.log).
    """
    with open(logfile,'a') as comlog:
        if message == "":
            comlog.write(tstamp() + "\t" + str(sys.executable) + " " + " ".join(sys.argv) + "\n")
        else:
            comlog.write(tstamp() + "\t" + str(sys.executable) + " " + " ".join(sys.argv) + "\n" + "          " + message.rstrip().replace("\n","\n          ") + "\n")
        

def log_message(message="", logfile="./messages.log"):
    """Record timestamp and message.
    
    Records timestamp and message to specified log-file.
    
    Args:
        message(str): Message to record in the log. (Default empty)
        logfile(str): File to write to. (./messages.log)
    """
    with open(logfile,'a') as comlog:
        comlog.write(tstamp() + "\t" + message.rstrip().replace("\n","\n          ") + "\n")
    


######################
##### EXECUTABLE #####
######################


if __name__ == "__main__":
    
    if sys.argv[1] == "-e":
        # Log command and run it.
        log_message(message = " ".join(sys.argv[2:]), logfile = "commands.log")
        subprocess.call(sys.argv[2:], stdout = sys.stdout, shell = False)
    else:
        # Log message.
        try:
            log_message(message = sys.argv[1], logfile = sys.argv[2])
        except IndexError:
            # No logfile provided, use default.
            log_message(message = sys.argv[1])
    
    sys.exit(0)
    
#EOF