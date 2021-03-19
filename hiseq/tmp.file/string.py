#!/usr/bin/env python3

"""Functions for string, list, dict, ...
transform to another version 
attribution of the vectors


"""



import os



def unique_list(seq, sorted=True, idfun=None):
    """Get the unique items of list
    seq: a list with items
    sorted: whether sort the output(unique)

    get the unique of inlist

    see1: Markus
    remove duplicates from a list while perserving order
    https://stackoverflow.com/a/480227/2530783

    def f7(seq):
        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]

    see2: ctcherry
    https://stackoverflow.com/a/89202/2530783

    def f5(seq, idfun=None):
        # order preserving
        if idfun is None:
            def idfun(x): return x
        seen = {}
        result = []
        for item in seq:
            marker = idfun(item)
            # in old Python versions:
            # if seen.has_key(marker)
            # but in new ones:
            if marker in seen: continue
            seen[marker] = 1
            result.append(item)
        return result
    """
    if idfun is None:
        def idfun(x): return x # for None
    if not isinstance(seq, list):
        log.error('list required, but get {}'.format(type(seq)))
        out = [] # blank
    elif sorted is True:
        out = list(set(seq))
    else:
        seen = set()
        out = [x for x in seq if x not in seen and not seen.add(x)]
    return out



##########################
## check files          ##
##########################
def update_obj(obj, d, force=True, remove=False):
    """Update the object, by dict
    d: dict
    force: bool, update exists attributes
    remove: bool, remove exists attributes
    Update attributes from dict
    force exists attr
    """
    if remove is True:
        for k in obj.__dict__:
            delattr(obj, k)
    # add attributes
    if isinstance(d, dict):
        for k, v in d.items():
            if not hasattr(obj, k) or force:
                setattr(obj, k, v)
    return obj
