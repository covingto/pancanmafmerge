
import os, os.path, sys
import socket
import threading
import time
import json, glob, csv
import SocketServer
import logging, subprocess

logger = logging.getLogger('dispatch')
logger.addHandler(logging.NullHandler())
logger.setLevel(logging.DEBUG)

PACKAGE = os.path.dirname(os.path.abspath(__file__))

def recv_timeout(the_socket,timeout=5):
    #make socket non blocking
    the_socket.setblocking(0)
    
    #total data partwise in an array
    total_data=[];
    data='';
    
    #beginning time
    begin=time.time()
    while 1:
        #if you got some data, then break after timeout
        if total_data and (time.time()-begin > timeout or total_data[-1][-1] == '}'):
            break
        
        #if you got no data at all, wait a little longer, twice the timeout
        elif time.time()-begin > timeout*2:
            break
        
        #recv something
        try:
            data = the_socket.recv(8192)
            if data:
                total_data.append(data)
                #change the beginning time for measurement
                begin=time.time()
            else:
                #sleep for sometime to indicate a gap
                time.sleep(0.1)
        except:
            pass
    
    #join all parts to make final string
    return ''.join(total_data)

class DispatchTCPClientServer(SocketServer.ThreadingMixIn, SocketServer.TCPServer):
    allow_reuse_address = True # we let the address to be resused when we are done with it
    def __init__(self, server_address, bind_and_activate=True):
        class ThreadedTCPRequestHandler(SocketServer.BaseRequestHandler):
            parent = self
            def handle(self):
                data = json.loads(recv_timeout(self.request, timeout = 5))
                response = json.dumps(self.parent.process(data))
                self.request.sendall(response)
        if server_address is None:
            server_address = socket.gethostbyname(socket.gethostname()), 0
        SocketServer.TCPServer.__init__(self, server_address, ThreadedTCPRequestHandler, bind_and_activate=bind_and_activate)
        self.ip, self.port = self.server_address
        self.server_thread = threading.Thread(target = self.serve_forever)
        self.server_thread.daemon = True
        self.server_thread.start()
        
    ##
    # The client portion of the class.  This function sends dicts (encoded in json)
    # and returns the response.
    @staticmethod
    def client(ip, port, message):
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.connect((ip, port))
        try:
            sock.sendall(json.dumps(message))
            response = recv_timeout(sock)
            return json.loads(response)
        finally:
            sock.close()
    
    def process(self, data):
        raise NotImplementedError("the process function has not been implemented")

class Worker(DispatchTCPClientServer):
    logger = logging.getLogger('dispatch.Worker')
    logger.addHandler(logging.NullHandler())

    def __init__(self, dip, dport, server_address = None):
        DispatchTCPClientServer.__init__(self, server_address)
        
        self._shutdown = False
        self._job_key = None
        self.dip = dip
        self.dport = dport

        self.job = None
        self.job_list = None
        self.loghandle = None
        self.job_lock = threading.RLock()

        self.ping_timeout = 10
        self.ping_lock = threading.RLock()

        def _ping():
            while True:
                try:
                    with self.job_lock:
                        if self.job is None:
                            logger.info("running ping")
                            self.ping()
                except:
                    logger.exception("Hit an exception during ping")
                    self.shutdown()
                    self.server_close()
                    self._shutdown = True

                time.sleep(self.ping_timeout)

        _t = threading.Thread(target = _ping)
        _t.daemon = True
        _t.start()
        
        self.check_job_timeout = 10
        self.check_job_lock = threading.RLock()
        def _check_job():
            while True:
                try:
                    with self.job_lock:
                        if self.job is not None:
                            self.check_job()
                except:
                    logger.exception("Hit an exception during check job")
                time.sleep(self.check_job_timeout)
        _t = threading.Thread(target = _check_job)
        _t.daemon = True
        _t.start()

    ## 
    # process a request
    def process(self, request):
        try:
            if 'request' not in request:
                return {'action': 'reject', 'reason': 'no request'}
            elif request['request'] == 'check':
                return self.process_check(request)
            elif request['request'] == 'kill':
                return self.process_kill(request)
            else:
                return {'action': 'reject', 'reason': 'non-standard request'}
        except Exception as inst:
            return {'action': 'reject', 'reason': str(inst)}
    
    def process_check(self, request):
        return {'action': 'accept', 'jobkey': self._job_key}

    def process_kill(self, request):
        with self.job_lock:
            if self.job is not None:
                self.job.kill()
                self.loghandle.write("Job killed!\n")
                self.loghandle.close()
                self.job = None # clear the job, if it errored we should know about that by some other means
                return {'action': 'accept'}
            else:
                return {'action': 'accept'}

    ##
    # ping for a new job
    def ping(self):
        try:
            response = self.client(self.dip, self.dport, {'request': 'ping', 'host': self.ip, 'port':self.port})
            logger.info(response)
            # response will be the command to execute and the log file location
            if 'action' in response and response['action'] == 'reject':
                self.shutdown()
                self.server_close()
                self._shutdown = True
                return response
            elif not ('cmd' in response and 'log' in response and 'jobkey' in response):
                logger.error("response is not formatted correctly: %s", response)
                return response
            else:
                logger.info(response['cmd'])
                logger.info("Starting job, log will be in %s", response['log'])
                self.loghandle = open(response['log'], 'w')
                self.loghandle.write('Executing %s\n' % response['cmd'])
                self._job_key = response['jobkey']
                self.job = subprocess.Popen(response['cmd'], shell = True, stdout = self.loghandle, stderr = self.loghandle)
                return response
        except:
            logger.exception("Hit an exception during ping")
            raise

    def check_job(self):
        self.job.poll()
        if self.job.returncode is not None:
            logger.info("Job complete, return status was %s", self.job.returncode)
            self.loghandle.write("Job complete\n")
            self.loghandle.close()
            retcode = self.job.returncode
            jkey = self._job_key
            self._job_key = None
            self.job = None # clear the job, if it errored we should know about that by some other means
            self.client(self.dip, self.dport, {'request': 'done', 'returncode': retcode, 'jobkey': jkey})

    def is_shutdown(self):
        return self._shutdown

class Dispatcher(DispatchTCPClientServer):
    logger = logging.getLogger('dispatch.Worker')
    logger.addHandler(logging.NullHandler())
    def __init__(self, server_address = ('', 46906)):
        DispatchTCPClientServer.__init__(self, server_address)
        
        self.job_list = {}
        self.problems = {}
        self.job_list_lock = threading.RLock()

        self._shutdown = False

        def _job_monitor():
            time.sleep(300)
            while True:
                if len(self.job_list) < 1:
                    self.shutdown()
                    self.server_close()
                    self._shutdown = True
                time.sleep(300)

        _t = threading.Thread(target = _job_monitor)
        _t.daemon = True
        _t.start()

    def is_shutdown(self):
        return self._shutdown

    def process(self, request):
        try:
            if 'request' not in request:
                return {'action': 'reject', 'reason': 'no request'}
            elif request['request'] == 'ping':
                return self.process_ping(request)
            elif request['request'] == 'queue':
                return self.process_queue(request)
            elif request['request'] == 'done':
                return self.process_done(request)
            elif request['request'] == 'problems':
                return self.process_problems(request)
            elif request['request'] == 'status':
                return self.process_status(request)
            else:
                return {'action': 'reject', 'reason': 'non-standard request'}
        except Exception as inst:
            return {'action': 'reject', 'reason': str(inst)}
    
    def process_status(self, request):
        return {'action': 'accepted',
                'data': [{
                    'key': k, 
                    'resultpath': job.resultpath, 
                    'cmd': job.cmd, 'log': job.log} for k, job in self.job_list.items()],
                'problems': [{'key': k, 'resultpath': job.resultpath, 'cmd': job.cmd, 'log': job.log} for k, job in self.problems.items()]
                }

    def process_queue(self, request):
        for k in ('resultpath', 'cmd', 'log', 'jobkey'):
            if k not in request:
                logger.info("queue request missing %s", k)
                return {'action': 'reject', 'reason': 'missing argument %s' % k}
        with self.job_list_lock:
            if request['jobkey'] in self.job_list:
                if request.get('overwrite', False):
                    # if the job is running we kill it
                    job = self.job_list[request['jobkey']]
                    if job.status() == 'running':
                        logger.info("Killing running job")
                        self.client(job.host(), job.port(), {'request': 'kill'})
                else:
                    logger.info("Ignoring job queue request: %s", request['jobkey'])
                    return {'action': 'accepted'}
            self.problems.pop(request['jobkey'], None) # remove the new queue from the problems list since it isn't a problem any more...

            self.job_list[request['jobkey']] = Job(request.get('dpath'), request.get('resultpath'), request.get('cmd'), request.get('log'))
            logger.info("Appended job")
            return {"action": 'accepted'}

    def process_problems(self, request):
        return {'action': 'accepted',
                'data': [{'key': k, 'resultpath': job.resultpath, 'cmd': job.cmd, 'log': job.log} for k, job in self.problems.items()]}

    def process_ping(self, request):
        with  self.job_list_lock:
            for k, job in self.job_list.items():
                if job.status() == 'running':
                    continue
                logger.info("Sending job command")
                job.set_running(**request)
                # return and terminate the iteration, we only set one job
                return {'action': 'accepted', 'cmd': job.get_cmd(), 'log': job.get_log(), 'jobkey': k}
            return {'action': 'reject', 'reason': 'no jobs'}

    def process_done(self, request):
        with self.job_list_lock:
            if 'jobkey' not in request:
                return {'action': 'reject', 'reason': 'no jobkey'}
            else:
                job = self.job_list.pop(request['jobkey'], None)
                if job is not None and request.get('returncode', 0) is not 0:
                    self.problems[request['jobkey']] = job
                return {'action': 'accepted'}

class PathFinder(object):
    ##
    # @param pathmap_fpath - file path to the mapping
    # @param pathmap_pkey - column representing the primary key for the pathmap, records are indexed by this
    # @param caller_map - maps column names in the mapping file with caller keys used by the merge {'broad_indelocator_vcf': 'INDELOCATOR', ...}
    def __init__(self, pathmap_fpath, pathmap_pkey, basedir, caller_map):
        assert isinstance(caller_map, dict), "caller_map was not a dict"
        assert os.path.isdir(basedir), "basedir does not exist"
        self._caller_map = caller_map
        self.basedir = basedir
        with open(pathmap_fpath, 'r') as fi:
            reader = csv.DictReader(fi, delimiter = '\t')
            self._map = {r[pathmap_pkey].strip('/'):r for r in reader}
        self.keys = self._map.keys()
    
    ##
    # @param key - the key used to extract records from the _map, if no match returns None
    # @returns a tuple of paths and caller keys (lists) or None, None
    def get_paths(self, key):
        if key not in self._map:
            return None, None
        else:
            r = self._map.get(key)
            paths = []
            callers = []
            for k, v in self._caller_map.items():
                fpath = os.path.join(self.basedir, r[k])
                if not os.path.isfile(fpath):
                    raise ValueError("%s is not a path" % fpath)
                paths.append(fpath)
                callers.append(v)
            return paths, callers

class Job(object):
    def __init__(self, jobkey, resultpath, cmd, log):
        self.jobkey = jobkey
        self.resultpath = resultpath
        self.cmd = cmd
        self.log = log
        self._status = 'idle'
        self._host = None
        self._port = None
   
    def status(self):
        return self._status

    def host(self):
        return self._host

    def port(self):
        return self._port

    def set_running(self, host, port, **kwargs):
        self._host = host
        self._port = port
        self._status = 'running'

    def get_cmd(self):
        return self.cmd

    def get_log(self):
        return self.log

    @staticmethod
    def dispatch(jobkey, resultpath, fmaps):
        logger.info("Processing %s", jobkey)
        outdir = os.path.join(os.path.abspath(resultpath), jobkey)
        output = os.path.join(outdir, 'merged.maf')
        tmpdir = os.path.join(outdir, 'tmp')
        if os.path.isfile(output):
            raise ValueError("%s has already been generated, remove if you want to run again", output)
        else:
            # make the directories that we need to write to
            for d in (outdir, tmpdir):
                if not os.path.isdir(d):
                    os.makedirs(d)
        # get the callers
        vcfs = []
        callers = []
        for fmap in fmaps:
            # fmap is a PathFinder instance that will allow us to map a jobkey to a set of file paths
            pl, cl = fmap.get_paths(jobkey)
            if pl is None:
                raise ValueError("%s is not a key in a file map" % jobkey)
            vcfs += pl
            callers += cl
        
        CMD = 'python %(PACKAGE)s/merge.py --vcfs %(vcfs)s --callers %(callers)s --tmpdir %(tmpdir)s %(output)s' % {
            'output': output,
            'tmpdir': tmpdir,
            'callers': ' '.join(callers),
            'vcfs': ' '.join(vcfs),
            'PACKAGE': PACKAGE}
        
        return CMD, output + '.log'

##
# start a worker instance, the worker will work until it is shutdown
def start_worker(args):
    if not args.dip or not args.dport:
        logger.error("Must specify --dip and --dport args")
        sys.exit(1)

    worker = Worker(dip = args.dip, dport = args.dport)
    
    while not worker.is_shutdown():
        time.sleep(10)
##
# build a job, the Job is just a container of variables
def build_jobs(jobkeys, resultdir, fmaps):
    with open("problem-keys.txt", 'w') as fo:
        writer = csv.DictWriter(fo, fieldnames = ['jobkey', 'reason'], delimiter = '\t')
        writer.writeheader()
        for jobkey in jobkeys: 
            try:
                t = Job.dispatch(jobkey, resultdir, fmaps)
                if t is not None:
                    cmd, log = t
                    yield Job(jobkey, resultdir, cmd, log)
            except Exception as inst:
                writer.writerow({'jobkey': jobkey, 'reason': str(inst)})


def build_fmaps(cpath):
    fmaps = []
    with open(cpath, 'r') as fi:
        config = json.load(fi)
        for mapping in config['fmaps']:
            # looks like: 
            #   {'basepath': path,
            #    'fmap': path,
            #    'pkey': key,
            #    'mapping': [[column_key, caller_key], [column_key, caller_key], ...]
            # def __init__(self, pathmap_fpath, pathmap_pkey, basedir, caller_map):
            fmaps.append(PathFinder(mapping['fmap'], mapping['pkey'], mapping['basepath'], mapping['mapping']))
    return fmaps
            
##
# starts a dispatcher, this will run until we run out of jobs
# future calls to queue may add more jobs to our list
def start_dispatcher(args):
    dispatcher = Dispatcher()
    args.dip = dispatcher.ip
    args.dport = dispatcher.port
    if args.config:
        queue(args)
    while not dispatcher.is_shutdown():
        time.sleep(10)
    logger.info("Started dispatcher: host: %s, port: %s", dispatcher.ip, dispatcher.port)

def queue(args):
    if not (args.dip and args.dport and args.resultdir and args.config):
        logger.error("Must specify all options, see help")
        sys.exit(1)
    fmaps = build_fmaps(args.config)
    for job in build_jobs(fmaps[0].keys, args.resultdir, fmaps):
        response = DispatchTCPClientServer.client(args.dip, args.dport, {
            'request': 'queue',
            'jobkey': job.jobkey,
            'resultpath': job.resultpath,
            'cmd': job.cmd,
            'log': job.log})
        if not response.get('action', None) == 'accepted':
            logger.error("action not accepted %s", response)
            sys.exit(1)

if __name__ == '__main__':
    import argparse
    ch = logging.StreamHandler(stream=sys.stdout)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)  
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers()

    parser_worker = subparser.add_parser('worker', help = 'start a worker node')
    parser_dispatcher = subparser.add_parser('dispatcher', help = 'start a dispatcher node')
    parser_queue = subparser.add_parser('queue', help = 'add records to a queue')

    parser_worker.add_argument('--dip', type = str, help = 'dispatcher ip')
    parser_worker.add_argument('--dport', type = int, help = 'dispatcher port number')
    parser_worker.set_defaults(func = start_worker)

    parser_dispatcher.add_argument('--resultdir', type = str, help = 'result dir')
    parser_dispatcher.add_argument('--config', type = str, help = 'config file path')
    parser_dispatcher.set_defaults(func = start_dispatcher)

    parser_queue.add_argument('--dip', type = str, help = 'dispatcher ip')
    parser_queue.add_argument('--dport', type = int, help = 'dispatcher port number')
    parser_queue.add_argument('--resultdir', type = str, help = 'result dir')
    parser_queue.add_argument('--config', type = str, help = 'config file path')
    parser_queue.set_defaults(func = queue)
    
    args = parser.parse_args()

    args.func(args)
