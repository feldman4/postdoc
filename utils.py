import time
import re

def timestamp(filename='', fmt='%Y%m%d_%H%M%S', sep='.'):
    stamp = time.strftime(fmt)
    pat= r'(.*)\.(.*)'
    match = re.findall(pat, filename)
    if match:
        return sep.join([match[0][0], stamp, match[0][1]])
    elif filename:
        return sep.join([filename, stamp])
    else:
        return stamp