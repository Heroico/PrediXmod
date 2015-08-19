__author__ = 'heroico'

import os
import errno

def ensure_file_path(path):
    if not "/" in path:
        return

    dirname = os.path.dirname(path)
    try:
        print "building file path "+dirname
        os.makedirs(dirname)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            print "error building file path "+dirname
            raise
        else:
            print "Already ensured file path "+path

def ensure_folder_path(path):
    if not "/" in path:
        return

    dirname = path
    try:
        print "building folder path "+dirname
        os.makedirs(dirname)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            print "error building folder path "+dirname
            raise
        else:
            print "Already ensured folder path "+path