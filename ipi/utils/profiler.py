#!/usr/bin/env python2

import threading
import time
import re

TIMING = time.time
PROF_LOCK = threading.Lock()

profile = {}
function_dict = {}

class Timer(object):
    """ Class that create a timer and stores timing.
    """
    def __init__(self, name=None):
        self.ncalls = 0
        self.ntime = 0
        if name:
            profile[name] = self

    def start(self):
        PROF_LOCK.acquire()
        self.ncalls += 1
        self.ntime -= TIMING()
        PROF_LOCK.release()

    def stop(self):
        PROF_LOCK.acquire()
        self.ntime += TIMING()
        PROF_LOCK.release()

    def __str__(self):
        msg = ' \t=>\t '+ str(self.ncalls)
        if self.ntime > 0:
            msg += ' \t  '+str(self.ntime)
        else:
            msg += ' \t  NA'
        return msg

FILENAME_REGEX = re.compile('.*/(ipi/.*).py')

def timethis(func):
    """ Decorator to timing functions.
    """
    def inner(*args, **kwargs):
        if __debug__:
            try:
                name = function_dict[func.func_code]
            except KeyError:
                code_obj = func.func_code
                filepath = FILENAME_REGEX.match(str(code_obj.co_filename)).group(1)
                fileline = str(code_obj.co_firstlineno)
                funcname = str(code_obj.co_name)
                name = funcname+'@'+filepath+':'+fileline
                function_dict[func.func_code] = name
                profile[name] = Timer()
            profile[name].start()
            result = func(*args, **kwargs)
            profile[name].stop()
            return result
        else:
            return func(*args, **kwargs)
    return inner
