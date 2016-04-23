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

#
# Start server on port 2306, with link and graph files:
#   python mccortex-server.py 2306 --coverages --edges -p l.ctp.gz a.ctx b.ctx
# Query a kmer:
#   curl localhost:2096/CAGTGGCCA
# Response:
#   { "key": "CAGTGGCCA", "colours": [1], "left": "T", "right": "T", "edges": "88", "links": [] }
#

def check_mccortex_alive(proc):
    if proc.poll() is not None:
        print("McCortex quit [%i]" % proc.returncode, file=sys.stdout)
        sys.exit(1)

def query_mccortex(proc,kmer):
    print(kmer, file=proc.stdin)
    proc.stdin.flush()
    check_mccortex_alive(proc)
    line = proc.stdout.readline()
    # Trim off prompt text
    if line[0:2] == "> ": line = line[2:len(line)]
    check_mccortex_alive(proc)
    return line

# when we start mccortex we set it to ignore interrupt signal, we handle it.
def preexec_function():
    # Ignore the SIGINT signal by setting the handler to the standard
    # signal handler SIG_IGN.
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def start_mccortex(extra_args):
    script_dir = os.path.dirname(os.path.realpath(__file__))

    # Adding two lists together appends one to the other
    try:
        proc = Popen([script_dir+"/../bin/mccortex31", "server", "--single-line"] + extra_args,
                      stdin=PIPE, stdout=PIPE, universal_newlines=True,
                      preexec_fn = preexec_function)
    except Exception as e:
        print("Couldn't start McCortex: ",e,file=sys.stdout)
        sys.exit(1)

    # Give test query to check it works
    check_mccortex_alive(proc)
    resp = query_mccortex(proc, "hi")

    return proc

def stop_mccortex(mccortex):
    print("q\n", file=mccortex.stdin)
    mccortex.stdin.close()
    # sleep until process has closed
    while mccortex.poll() is None:
        time.sleep(1)
    print("McCortex exited with:",str(mccortex.poll()))

def test_mccortex(args):
    proc = start_mccortex(args)
    print("Got: ", query_mccortex(proc, "CACTGATGA"), end='')
    print("Got: ", query_mccortex(proc, "CCACTGATG"), end='')
    try:
        while True:
            pass
    except (KeyboardInterrupt,SystemExit):
        print("Got exit signal")
    stop_mccortex(proc)

def start_web_server(port,args):
    mccortex = start_mccortex(args)

    class mccortexHTTPServer(BaseHTTPRequestHandler):
        def do_GET(self):
            self.send_response(200)
            self.send_header("Content-type", "text/html")
            self.end_headers()
            if len(self.path) < 4 or len(self.path) > 300:
                jsonstr = "{\"error\": \"Webserver: bad query\"}\n"
            else:
                jsonstr = query_mccortex(mccortex, self.path[1:])
            self.wfile.write(jsonstr.encode("UTF-8"))

    try:
        httpd = socketserver.TCPServer(("", port), mccortexHTTPServer)
    except Exception as e:
        print("Cannot start HTTP server: ",e,file=sys.stderr)
        sys.exit()

    print("serving at port", port)

    try:
        httpd.serve_forever()
    except (KeyboardInterrupt,SystemExit):
        print(" Got exit signal")

    httpd.server_close()
    stop_mccortex(mccortex)
    print("closed. bye.")


def main():
    if len(sys.argv) < 3 or not sys.argv[1].isdigit():
        print("usage: %s <port> [mccortex args]" % (sys.argv[0]))
        print("  e.g  %s 1888 -m 2G graph.ctx" % (sys.argv[0]))
        sys.exit(-1)

    port = int(sys.argv[1])
    args = sys.argv[2:]

    # test_mccortex(args)
    start_web_server(port,args)

if __name__ == '__main__':
    main()
