#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import time
import signal

from subprocess import Popen, PIPE

# These should work on python2 after 'pip install --user future'
from http.server import BaseHTTPRequestHandler
import socketserver

# Opening browser windows
import webbrowser
import threading

#
# Start server on port 2306, with link and graph files:
#   python mccortex-server.py 2306 --coverages --edges -p l.ctp.gz a.ctx b.ctx
# Query a kmer:
#   curl localhost:2096/CAGTGGCCA
# Response:
#   { "key": "CAGTGGCCA", "colours": [1], "left": "T", "right": "T", "edges": "88", "links": [] }
#

def usage(err=None):
    if err is not None: print(err,file=sys.stderr)
    print("usage: %s <port> [-W,--web|-k,--kmer <K>] [mccortex args]" % (sys.argv[0]))
    print("  e.g  %s 1888 -m 2G graph.ctx" % (sys.argv[0]))
    sys.exit(-1)

# when we start mccortex we set it to ignore interrupt signal, we handle it.
def preexec_function():
    # Ignore the SIGINT signal by setting the handler to the standard
    # signal handler SIG_IGN.
    signal.signal(signal.SIGINT, signal.SIG_IGN)

class McCortexProc:
    def __init__(self,kmer=31,extra_args=[]):
        self.lock = threading.Lock()
        script_dir = os.path.dirname(os.path.realpath(__file__))
        try:
            self.proc = Popen([script_dir+"/../../bin/mccortex", str(kmer),
                               "server", "--single-line"] + extra_args,
                              stdin=PIPE, stdout=PIPE, universal_newlines=True,
                              preexec_fn = preexec_function)
        except Exception as e:
            print("Couldn't start McCortex: ",e,file=sys.stdout)
            sys.exit(1)
    def query_mccortex(self,query):
        with self.lock:
            print(query, file=self.proc.stdin)
            self.proc.stdin.flush()
            self.check_mccortex_alive()
            line = self.proc.stdout.readline()
            # Trim off prompt text
            if line[0:2] == "> ": line = line[2:len(line)]
            self.check_mccortex_alive()
            return line
    def check_mccortex_alive(self):
        if self.proc.poll() is not None:
            print("McCortex quit [%i]" % self.proc.returncode, file=sys.stdout)
            sys.exit(1)
    def stop_mccortex(self):
        print("q\n", file=self.proc.stdin)
        self.proc.stdin.close()
        # sleep until process has closed
        while self.proc.poll() is None: time.sleep(1)
        print("McCortex exited with:",str(self.proc.poll()))

def test_mccortex(args):
    mcproc = McCortexProc(9,args)
    print("Got: ", mcproc.query_mccortex("CACTGATGA"), end='')
    print("Got: ", mcproc.query_mccortex("CCACTGATG"), end='')
    try:
        while True: pass
    except (KeyboardInterrupt,SystemExit):
        print(" Got exit signal")
    mcproc.stop_mccortex()

def launch_webpage(url,delay=1):
    time.sleep(delay)
    webbrowser.open(url)

def start_web_server(port,args,kmer=31,webpage=None):
    mcproc = McCortexProc(kmer,args)

    class mccortexHTTPServer(BaseHTTPRequestHandler):
        def do_GET(self):
            if len(self.path) < 4 or len(self.path) > 300:
                jsonstr = "{\"error\": \"Webserver: bad query\"}\n"
            else:
                jsonstr = mcproc.query_mccortex(self.path[1:])
            self.send_response(200)
            self.send_header("Content-type", "text/html")
            self.end_headers()
            self.wfile.write(jsonstr.encode("utf-8"))

    # Reuse old open port (http://stackoverflow.com/a/10614360/431087)
    # to avoid error: [Errno 48] Address already in use
    class MyTCPServer(socketserver.TCPServer, socketserver.ThreadingMixIn):
        allow_reuse_address = True

    try:
        httpd = MyTCPServer(("", port), mccortexHTTPServer)
    except Exception as e:
        print("Cannot start HTTP server: ",e,file=sys.stderr)
        sys.exit()

    # Open webpage if requested, after pausing for 1 sec
    if webpage is not None:
        threading.Thread(target=launch_webpage,args=(webpage,)).start()

    print("serving at port", port)
    try: httpd.serve_forever()
    except (KeyboardInterrupt,SystemExit):
        print(" Got exit signal")
    httpd.server_close()
    mcproc.stop_mccortex()
    print("Goodbye.")

def main():
    kmer = 31
    webpage = None
    if len(sys.argv) < 3 or not sys.argv[1].isdigit(): usage()
    port = int(sys.argv[1])
    while len(sys.argv) > 2:
        if sys.argv[2] in ['-W','--web']:
            webpage = "http://127.0.0.1:"+str(port)
            del(sys.argv[2])
        elif sys.argv[2] in ['-k','--kmer'] and len(sys.argv) > 3:
            try: kmer = int(sys.argv[3])
            except ValueError: usage("Invalid --kmer value: "+sys.argv[3])
            del(sys.argv[2:4])
        else: break
    args = sys.argv[2:]
    # test_mccortex(args)
    start_web_server(port,args,kmer,webpage)

if __name__ == '__main__':
    main()
