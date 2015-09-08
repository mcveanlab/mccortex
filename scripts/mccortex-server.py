#!/usr/bin/env python3

from subprocess import Popen, PIPE
import os
import sys
import time
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

def query_mccortex(proc,kmer):
    print(kmer, file=proc.stdin)
    proc.stdin.flush()
    line = proc.stdout.readline()
    # line = line.rstrip()
    if proc.poll() is not None:
        print("McCortex quit [%i]" % proc.returncode, file=sys.stdout)
        sys.exit(1)
    return line

def start_mccortex(extra_args):
    script_dir = os.path.dirname(os.path.realpath(__file__))

    # Adding two lists together appends one to the other
    try:
        proc = Popen([script_dir+"/../bin/mccortex31", "server", "--single-line"] + extra_args,
                      stdin=PIPE, stdout=PIPE, universal_newlines=True)
    except Exception as e:
        print("Couldn't start McCortex: ",e,file=sys.stdout)
        sys.exit(1)

    # Give test query to check it works
    resp = query_mccortex(proc, "hi")

    return proc

def test_mccortex():
    proc = start_mccortex(["../../tests/thread_pe_short/genome.k9.ctx"])
    print("Got: ", query_mccortex(proc, "CACTGATGA"), end='')
    print("Got: ", query_mccortex(proc, "CCACTGATG"), end='')
    proc.stdin.close()
    print("McCortex exit: ",proc.wait())

def start_web_server():
    port = int(sys.argv[1])
    mccortex = start_mccortex(sys.argv[2:])

    class mccortexHTTPServer(BaseHTTPRequestHandler):
        def do_GET(self):
            jsonstr = query_mccortex(mccortex, self.path[1:])
            self.send_response(200)
            self.send_header("Content-type", "text/html")
            self.end_headers()
            self.wfile.write(bytes(jsonstr,'UTF-8'))

    try:
        httpd = socketserver.TCPServer(("", port), mccortexHTTPServer)
    except Exception as e:
        print("Cannot start HTTP server: ",e,file=sys.stderr)
        sys.exit()

    print("serving at port", port)

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass

    httpd.server_close()
    mccortex.stdin.close()
    print("wait: ", mccortex.wait(timeout=10))
    print("closed. bye.")


def main():
    if len(sys.argv) < 3 or not sys.argv[1].isdigit():
        print("usage: %s <port> [mccortex args]" % (sys.argv[0]))
        sys.exit(-1)

    # test_mccortex()
    start_web_server()

if __name__ == '__main__':
    main()
