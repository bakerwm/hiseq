

import sys
import logging
import functools


# Test for Decorator

class Logger(object):
    def __init__(self, level='INFO'):
        logging.basicConfig(
            format='[%(asctime)s %(levelname)s] %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            stream=sys.stdout)
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(level)

    def __call__(self, fn):
        @functools.wraps(fn)
        def decorated(*args, **kwargs):
            try:
                self.logger.info('{0} - {1} - {2}'.format(fn.__name__, args, kwargs))
                result = fn(*args, **kwargs)
                # self.logger.info(result)
                return result
            except Exception as ex:
                self.logger.info('Exception {0}'.format(ex))
                raise ex
            return result
        return decorated



@Logger('INFO')
def funcA(a, b, c):
    return a + b + c


@Logger('INFO')
def funcB(a, b, c):
    return a + b + c

funcA(1, 2, 3)

funcB(3, 4, 5)
