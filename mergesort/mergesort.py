'''
Created on May 6, 2016

@author: covingto
'''

import json, tempfile
import shutil, os.path

class PeekWrapper(object):
    def __init__(self, myiter):
        self.iter = myiter
        self._next = None
        self.take()
    
    def peek(self):
        return self._next
    
    def take(self):
        old = self._next
        try:
            self._next = self.iter.next()
        except StopIteration:
            self._next = None
        return old
    
    def __iter__(self):
        return self 
    
    def next(self):
        if self._next is not None:
            return self.take()
        else:
            raise StopIteration
    
def chunk(myiter, size):
    batch = []
    for i in myiter:
        batch.append(i)
        if len(batch) > size:
            yield batch
            batch = []
    yield batch

def mergepeek(initers, keyfun):
    inpeeks = [PeekWrapper((r for r in i if r != '')) for i in initers]
    inpeeks = sorted(inpeeks, key = lambda x: keyfun(x.peek()))
    nsorts = 1
    while len(inpeeks) > 1:
        if inpeeks[0].peek() is None:
            inpeeks = sorted(inpeeks[1:], key = lambda x: keyfun(x.peek()))
            continue # iterate again
        if keyfun(inpeeks[0].peek()) > keyfun(inpeeks[1].peek()):
            print 'mergepeek sorting number %d' % nsorts
            nsorts += 1
            inpeeks = sorted(inpeeks, key = lambda x: keyfun(x.peek()))
        yield inpeeks[0].take()
    lastpeek = inpeeks[0]
    while lastpeek.peek() is not None:
        yield lastpeek.take()
    
        
def mergesort(myiter, keyfun, maxrecords = 100000, maxpaths = 10):
    tmpdir = tempfile.mkdtemp('.mergesort', dir=os.path.abspath('.'))
    fnum = 0
    fpaths = []
    for c in chunk(myiter, maxrecords):
        #print "chunk returned: %s" % str(c)
        _opath = os.path.join(tmpdir, 'tmp-%d.jsonl' % fnum)
        fnum += 1
        csort = sorted([{'key': keyfun(cc), 'value': cc} for cc in c], key = lambda x: x['key'])
        with open(_opath, 'w') as fo:
            fo.write('\n'.join( [json.dumps(cc) for cc in csort] ))
        fpaths.append(_opath)
        if len(fpaths) > maxpaths:
            # now read back in and make one file ...
            _opath = os.path.join(tmpdir, 'tmp-%d.jsonl' % fnum)
            fnum += 1
            with open(_opath, 'w') as fo:
                fhandles = [open(p, 'r') for p in fpaths]
                for r in mergepeek(fhandles, lambda x: json.loads(x)['key']):
                    fo.write(r.strip() + '\n')
                # clean up the open file handles
                for h in fhandles:
                    h.close()
                for p in fpaths:
                    os.remove(p) # we won't use these again
            fpaths = [_opath]
    fhandles = [open(p, 'r') for p in fpaths]
    for r in mergepeek(fhandles, lambda x: json.loads(x)['key']):
        yield json.loads(r)['value']
    for h in fhandles:
        h.close()
    shutil.rmtree(tmpdir, False)

def test():
    print "starting test"
    r = [v for v in mergesort(range(100), lambda x: x, maxrecords=10, maxpaths = 5)]
    assert r == [v for v in range(100)]
    print "Success"

