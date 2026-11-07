"""Stream a large remote file to stdout in sequential HTTP byte ranges, resuming where a dropped connection stopped
(download.cncb.ac.cn closes long transfers early). Stops quietly when the reader closes the pipe (e.g. `head`).
usage: python src/range_stream.py <url> [chunk_MB=64]"""
import sys
import time

import requests


def main(url, chunk=64 << 20):
    out = sys.stdout.buffer
    off, fails, total = 0, 0, None
    while total is None or off < total:
        try:
            r = requests.get(url, headers={"Range": f"bytes={off}-{off + chunk - 1}"}, stream=True, timeout=60)
            if r.status_code == 416:
                return
            if r.status_code != 206:  # a server ignoring Range would resend from byte 0
                raise requests.RequestException(f"HTTP {r.status_code}")
            total = int(r.headers["Content-Range"].split("/")[1])
            for b in r.iter_content(1 << 20):
                out.write(b)
                off += len(b)  # advance per block, so a dropped connection resumes exactly here
            fails = 0
        except BrokenPipeError:
            return
        except (requests.RequestException, ConnectionError, ValueError, KeyError):
            fails += 1
            if fails > 120:  # ~2 h of one attempt a minute: CNGB refuses all connections for long spells
                raise
            time.sleep(min(60, 5 * fails))
    out.flush()


if __name__ == "__main__":
    try:
        main(sys.argv[1], int(sys.argv[2]) << 20 if len(sys.argv) > 2 else 64 << 20)
    except BrokenPipeError:
        pass
