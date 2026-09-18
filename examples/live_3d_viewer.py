#!/usr/bin/env python3
"""Serve a dependency-free live 3-D slice viewer for cfd_solvers JSON frames."""
from __future__ import annotations

import argparse
import http.server
import pathlib
import urllib.parse


HTML = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>cfd_solvers live 3-D viewer</title>
<style>
html,body{margin:0;height:100%;background:#111;color:#ddd;font:14px system-ui,sans-serif;overflow:hidden}
body{display:grid;grid-template-rows:auto 1fr}header{display:flex;gap:12px;align-items:center;flex-wrap:wrap;padding:10px 14px;background:#20252b}
label{display:flex;gap:6px;align-items:center}select,input{accent-color:#54b6d6}#status{margin-left:auto;color:#9dd9e8}
main{display:grid;place-items:center;min-height:0;padding:12px}canvas{max-width:100%;max-height:100%;image-rendering:pixelated;background:#050505;border:1px solid #3e4c55}
</style>
</head>
<body>
<header><label>axis <select id="axis"><option value="z">Z</option><option value="y">Y</option><option value="x">X</option></select></label>
<label>slice <input id="slice" type="range" min="0" max="0" value="0"></label><span id="sliceValue">0</span><span id="status">waiting for state.json</span></header>
<main><canvas id="view"></canvas></main>
<script>
const canvas=document.getElementById('view'),ctx=canvas.getContext('2d'),axis=document.getElementById('axis'),slice=document.getElementById('slice'),sliceValue=document.getElementById('sliceValue'),status=document.getElementById('status');
let state=null;
function shape(){const [nx,ny,nz]=state.dimensions;return axis.value==='x'?[ny,nz,nx]:axis.value==='y'?[nx,nz,ny]:[nx,ny,nz]}
function index(x,y,z){const [nx,ny]=state.dimensions;return (z*ny+y)*nx+x}
function value(u,v,s){const [nx,ny,nz]=state.dimensions;return axis.value==='x'?state.values[index(s,u,v)]:axis.value==='y'?state.values[index(u,s,v)]:state.values[index(u,v,s)]}
function solid(u,v,s){const [nx,ny,nz]=state.dimensions;const a=state.solid||[];return a.length===0?false:(axis.value==='x'?a[index(s,u,v)]:axis.value==='y'?a[index(u,s,v)]:a[index(u,v,s)])!==0}
function color(t){const r=Math.round(25+230*t),g=Math.round(75+165*(1-Math.abs(2*t-1))),b=Math.round(235-190*t);return [r,g,b]}
function draw(){if(!state)return;const [w,h,depth]=shape();slice.max=Math.max(0,depth-1);slice.value=Math.min(Number(slice.value),depth-1);sliceValue.textContent=slice.value;canvas.width=w;canvas.height=h;const image=ctx.createImageData(w,h);let lo=Infinity,hi=-Infinity;for(const v of state.values){lo=Math.min(lo,v);hi=Math.max(hi,v)}const span=hi-lo||1;for(let v=0;v<h;v++)for(let u=0;u<w;u++){const p=(v*w+u)*4;const isSolid=solid(u,v,Number(slice.value));if(isSolid){image.data[p]=42;image.data[p+1]=48;image.data[p+2]=53;image.data[p+3]=255;continue}const t=Math.max(0,Math.min(1,(value(u,v,Number(slice.value))-lo)/span)),[r,g,b]=color(t);image.data[p]=r;image.data[p+1]=g;image.data[p+2]=b;image.data[p+3]=255}ctx.putImageData(image,0,0);status.textContent=`step ${state.step} | ${state.field} ${lo.toPrecision(4)} .. ${hi.toPrecision(4)} | ${state.dimensions.join(' x ')}`}
function updateControls(){if(!state)return;const [w,h,depth]=shape();slice.max=Math.max(0,depth-1);slice.value=Math.min(Number(slice.value),depth-1);draw()}
async function poll(){try{const response=await fetch(`state.json?${Date.now()}`,{cache:'no-store'});if(response.ok){const next=await response.json();if(!state||next.step!==state.step){state=next;updateControls()}}}catch(error){status.textContent='waiting for state.json'}setTimeout(poll,250)}
axis.addEventListener('change',updateControls);slice.addEventListener('input',draw);poll();
</script>
</body></html>
"""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", required=True, type=pathlib.Path, help="directory containing live/state.json")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", default=8765, type=int)
    args = parser.parse_args()
    root = args.dir.resolve()
    root.mkdir(parents=True, exist_ok=True)

    class Handler(http.server.SimpleHTTPRequestHandler):
        def __init__(self, *handler_args, **handler_kwargs):
            super().__init__(*handler_args, directory=str(root), **handler_kwargs)

        def do_GET(self):  # noqa: N802 - standard library handler API
            path = urllib.parse.urlsplit(self.path).path
            if path in ("", "/"):
                payload = HTML.encode("utf-8")
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.send_header("Content-Length", str(len(payload)))
                self.end_headers()
                self.wfile.write(payload)
                return
            if path == "/state.json":
                if not (root / "state.json").exists():
                    self.send_response(404)
                    self.end_headers()
                    return
            super().do_GET()

        def log_message(self, *_):
            return

    server = http.server.ThreadingHTTPServer((args.host, args.port), Handler)
    print(f"live 3-D viewer: http://{args.host}:{args.port}/")
    print(f"watching: {root / 'state.json'}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        return 0
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
