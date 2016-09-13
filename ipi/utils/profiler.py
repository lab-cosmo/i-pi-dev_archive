#!/usr/bin/env python2

import threading
import time

timing = time.clock
prof_lock = threading.Lock()

profile = {}

def start(name):
    prof_lock.acquire()
    if name in profile.keys():
        profile[name][0] += 1
        profile[name][1] -= timing()
    elif name not in profile.keys():
        profile[name] = [1, -timing(), True]
    prof_lock.release()

def stop(name):
    prof_lock.acquire()
    if name in profile.keys():
        profile[name][1] += timing()
    else:
        raise TimerError('@PROFILING: Timer for '+name+' not started!\n Try to stop a non-started clock\n')
    prof_lock.release()

def timethis(func):
    """ Decorator to timing functions.
    """
    def inner(*args, **kwargs):
        if __debug__:
            code_obj = func.func_code
            filepath = str(code_obj.co_filename)
            fileline = str(code_obj.co_firstlineno)
            funcname = str(code_obj.co_name)
            name = funcname+'@'+filepath+':'+fileline
            start(name)
            result = func(*args, **kwargs)
            stop(name)
            return result
        else:
            pass
    return inner


class TimerError(Exception):
    def __init__(self, msg):
        self.value = msg
    def __str__(self):
        return self.value


