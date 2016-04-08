
import dispatch_server
import unittest

class TestDispatchServer(unittest.TestCase):
    def setUp(self):
        self.dispatcher = dispatch_server.Dispatcher()
        

    def tearDown(self):
        # close the server
        self.dispatcher.shutdown()
        self.dispatcher.server_close()
        self.dispatcher._shutdown = True

    def test_queue(self):
        # build a job
        job = dispatch_server.Job('test', '.', 'echo "test"', 'test.log')
        # add some things to the queue
        with self.dispatcher.job_list_lock:
            self.dispatcher.job_list['test'] = job
        self.assertTrue('test' in self.dispatcher.job_list, "didn't add the test job to the queue")
        # add via the client
        for i in xrange(5):
            i = str(i)
            job = dispatch_server.Job('test%s' % i, '.', 'echo "test%s"' % i, 'test%s.log' % i)
            result = dispatch_server.DispatchTCPClientServer.client(self.dispatcher.ip, self.dispatcher.port, {
                'jobkey': job.jobkey,
                'request': 'queue',
                'resultpath': job.resultpath,
                'cmd': job.cmd,
                'log': job.log})
            self.assertTrue(result['action'] == 'accepted', "result not accepted: %s" % str(result))
            self.assertTrue('test%s' % i in self.dispatcher.job_list, "failed to add a job, job list is %s" % self.dispatcher.job_list.keys())
    
        # build path finder
        import csv
        with open('test.fmap.txt', 'w') as fo:
            writer = csv.DictWriter(fo, fieldnames = ['pkey', 'indel', 'snv'], delimiter = '\t')
            writer.writeheader()
            for i in xrange(5):
                writer.writerow({
                    'pkey': 'fmap%s' % i,
                    'indel': 'indel%s.txt' % i,
                    'snv': 'snv%s.txt' % i
                    })
        fmap = dispatch_server.PathFinder('test.fmap.txt', 'pkey', '.', dict([('indel', 'INDEL'), ('snv', 'SNV')]))
        for key in fmap._map.keys():
            p, c = fmap.get_paths(key)
            self.assertTrue(p is not None, "paths were none")
            self.assertTrue(c is not None, "callers were none")
            self.assertTrue(len(p) == len(c), "paths and callers are not the same length")
            self.assertTrue(len(p) == 2, "2 paths not returned (as expected by the test")
            # test that this can be converted to a Job
            j = dispatch_server.Job.dispatch(key, '.', [fmap])
            self.assertTrue(j is not None, "job.dispatch returned None")
            cmd, log = j
            self.assertTrue(isinstance(cmd, basestring), "cmd was not a string")
            self.assertTrue(isinstance(log, basestring), "log was not a string")
    
    def test_ping(self):
        job = dispatch_server.Job('test', '.', 'echo "test"', 'test.log')
        # add some things to the queue
        with self.dispatcher.job_list_lock:
            self.dispatcher.job_list['test'] = job
        self.assertTrue('test' in self.dispatcher.job_list, "didn't add the test job to the queue")
        response = dispatch_server.DispatchTCPClientServer.client(self.dispatcher.ip, self.dispatcher.port, {'request': 'ping', 'host': '0.0.0.0', 'port': '64941'})
        self.assertTrue(response is not None, "no response")
        self.assertTrue('action' in response, "response missing action: %s" % str(response))
        self.assertTrue(response.get('action') != 'reject', "action was rejected: %s" % str(response))
        for k in ('cmd', 'log', 'jobkey'):
            self.assertTrue(k in response, "response missing %s" % k)



    def test_worker(self):
        # start a worker
        worker = dispatch_server.Worker(self.dispatcher.ip, self.dispatcher.port)
        with worker.job_lock: # we lock here so that we can control the worker with more granularity instead of letting her go wild
            ## test a passing job
            job = dispatch_server.Job('test', '.', 'echo "test"', 'test.log')
            # add some things to the queue
            with self.dispatcher.job_list_lock:
                self.dispatcher.job_list['test'] = job
            self.assertTrue('test' in self.dispatcher.job_list, "didn't add the test job to the queue")
            response = worker.ping()
            self.assertTrue(response['action'] != 'reject', "response was rejected: %s" % str(response))
            self.assertTrue('cmd' in response, "response did not contain a cmd: %s" % str(response))
            self.assertTrue(worker.job is not None, "The worker didn't get a job")
            self.assertTrue('test' == worker._job_key, 'Worker did not get the right job key; %s' % worker._job_key)
            # allow the job to finish
            while True:
                worker.job.poll()
                if worker.job.returncode is not None:
                    worker.check_job()
                    break # break from the while loop, this means that the job is done
            self.assertTrue(worker.job is None, "The worker did not reset the job")
            self.assertFalse('test' in self.dispatcher.job_list, "The dispatcher still has the job, this should have been removed on done")
            self.assertFalse('test' in self.dispatcher.problems, "The return code should have been 0 for this test, but instead the key is in problems: %s" % str(self.dispatcher.problems))
            ## test a failing job
            job = dispatch_server.Job('test', '.', 'exit 10', 'test.log')
            # add some things to the queue
            with self.dispatcher.job_list_lock:
                self.dispatcher.job_list['test'] = job
            self.assertTrue('test' in self.dispatcher.job_list, "didn't add the test job to the queue")
            worker.ping()
            self.assertTrue('test' == worker._job_key, 'Worker did not get the right job key; %s' % worker._job_key)
            # allow the job to finish
            while True:
                worker.job.poll()
                if worker.job.returncode is not None:
                    worker.check_job()
                    break # break from the while loop, this means that the job is done
            self.assertTrue(worker.job is None, "The worker did not reset the job")
            self.assertFalse('test' in self.dispatcher.job_list, "The dispatcher still has the job, this should have been removed on done")
            self.assertTrue('test' in self.dispatcher.problems, "The return code should have been 10 for this test, but the key is not in problems: %s" % str(self.dispatcher.problems))
            # test a kill
            job = dispatch_server.Job('test', '.', 'exit 10', 'test.log')
            # add some things to the queue
            with self.dispatcher.job_list_lock:
                self.dispatcher.job_list['test'] = job
            self.assertTrue('test' in self.dispatcher.job_list, "didn't add the test job to the queue")
            worker.ping()
            self.assertTrue('test' == worker._job_key, 'Worker did not get the right job key; %s' % worker._job_key)
            worker.process_kill({})
            self.assertTrue(worker.job is None, "Kill was called and the job is not None")

        worker.shutdown()
        worker.server_close()
        worker._shutdown = True


    def test_communication(self):
        ## test a passing job
        job = dispatch_server.Job('test', '.', 'echo "test"', 'test.log')
        # add some things to the queue
        with self.dispatcher.job_list_lock:
            self.dispatcher.job_list['test'] = job
        self.assertTrue('test' in self.dispatcher.job_list, "didn't add the test job to the queue")
        ## test status
        response = dispatch_server.DispatchTCPClientServer.client(self.dispatcher.ip, self.dispatcher.port, {'request': 'status'})
        self.assertTrue('action' in response, 'response missing action')
        self.assertTrue(response.get('action') == 'accepted', 'action not accepted: %s' % str(response))
        self.assertTrue('data' in response, 'response missing data')
        self.assertTrue('problems' in response, 'response missing problems')
        self.assertTrue(len(response.get('data')) > 0, 'no data in data key: %s' % str(response))
        
        
if __name__ == "__main__":
    unittest.main()
