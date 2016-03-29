

import xml.sax
from Queue import Queue
from ftplib import FTP
import gzip
from threading import Thread
import os, os.path, sys

class DBSNPHANDLER(xml.sax.ContentHandler):
    def __init__(self, writer_queue):
        self.CurrentData = ''
        self.rsID = None
        self.validations = []
        self.writer_queue = writer_queue
        
        
    def startElement(self, tag, attributes):
        self.CurrentData = tag
        if self.CurrentData == 'Rs':
            self.rsID = 'rs' + attributes['rsId']
        elif self.CurrentData == 'Validation':
            self.validations = attributes.keys()
    
    def endElement(self, tag):
        if tag == 'Rs':
            for v in self.validations:
                self.writer_queue.put((self.rsID, v))
            self.validations = []

class DBSNPFTPCONN():
    def __init__(self, host, path):
        self.host = host
        self.path = path
        self.ftp_conn = FTP(host)
        self.ftp_conn.login()
        self.ftp_conn.cwd(path)
        self.xml_files = []
        def collect_xml_file(line):
            if '.xml.' in line:
                self.xml_files.append(line.split(' ')[-1])
        self.ftp_conn.dir(collect_xml_file)
        
    def get_xml(self, directory):
        for file_name in self.xml_files:
            output_file = os.path.join(directory, file_name)
            if not os.path.isfile(output_file):
                print "retr %s" % file_name
                self.ftp_conn.retrbinary('RETR %s' % file_name, open(output_file, 'wb').write)
            yield output_file

def main(args):
    writer_queue = Queue()
    def writer_f():
        while True:
            rsid, validation = writer_queue.get()
            try:
                args.OUTPUT.write('\t'.join([rsid, validation]) + '\n')
            finally:
                writer_queue.task_done()
    
    writer_t = Thread(target = writer_f)
    writer_t.daemon = True
    writer_t.start()
    
    parser_queue = Queue()
    def parse_f():
        while True:
            gzfile = parser_queue.get()
            try:
                parser = xml.sax.make_parser()
                parser.setFeature(xml.sax.handler.feature_namespaces, 0)
                handler = DBSNPHANDLER(writer_queue)
                parser.setContentHandler( handler )
                print "parsing %s" % gzfile
                with gzip.open(gzfile, 'rb') as fi:
                    parser.parse( fi )
            finally:
                parser_queue.task_done()
    
    parser_t = Thread(target = parse_f)
    parser_t.daemon = True
    parser_t.start()
    
    dbsnp_conn = DBSNPFTPCONN(args.host, args.path)
    for gzfile in dbsnp_conn.get_xml(args.tmp_dir):
        parser_queue.put(gzfile)
    
    parser_queue.join()
    writer_queue.join()

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--tmp-dir', default = '.', help = 'directory for storing downloaded xml files, may actually not be temporary')
    parser.add_argument('--host', type = str, default = 'ftp.ncbi.nih.gov', help = 'ftp host for the dbSNP xml files')
    parser.add_argument('--path', type = str, help = 'host path to directory of xml files e.g. snp/organisms/human_9606_b146_GRCh37p13/XML')
    parser.add_argument('OUTPUT', type = argparse.FileType('w'), help = 'output file for writing, writes in the form rsid [tab] validation status')
    
    args = parser.parse_args()
    
    main(args)
