#!/usr/bin/env python2

import threading
import time

timing = time.clock
prof_lock = threading.Lock()

profile = {}
function_dict = {}

class Timer(object):    
    def __init__(self):
        self.ncalls = 0
        self.ntime =0
        
    def start(self):
        prof_lock.acquire()
        self.ncalls += 1
        self.ntime -= timing()
        prof_lock.release()

    def stop(self):
        prof_lock.acquire()
        self.ntime += timing()
        prof_lock.release()

def timethis(func):
    """ Decorator to timing functions.
    """
    def inner(*args, **kwargs):
        if __debug__:
            try:
                name = function_dict[func.func_code]
            except KeyError:
                code_obj = func.func_code
                filepath = str(code_obj.co_filename)
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


class TimerError(Exception):
    def __init__(self, msg):
        self.value = msg
    def __str__(self):
        return self.value


