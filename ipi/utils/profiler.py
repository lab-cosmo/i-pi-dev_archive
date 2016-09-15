#!/usr/bin/env python2

import threading
import time

timing = time.clock
prof_lock = threading.Lock()

profile = {}
function_dict = {}

class Timer(object):    
    def __init__(self, name=None):
        self.ncalls = 0
        self.ntime =0
        if name:
            profile[name] = self
        
    def start(self):
        prof_lock.acquire()
        self.ncalls += 1
        self.ntime -= timing()
        prof_lock.release()

    def stop(self):
        prof_lock.acquire()
        self.ntime += timing()
        prof_lock.release()

    def __str__(self):
        msg = ' \t=>\t '+ str(self.ncalls)
        if self.ntime > 0:
            msg += ' \t  '+str(self.ntime)
        else:
            msg += ' \t  NA'
        return msg

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


