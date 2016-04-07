
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


    def test_communication(self):
        pass

if __name__ == "__main__":
    unittest.main()
