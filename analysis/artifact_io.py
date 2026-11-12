"""Commit completed local gzip artifacts to EOS and check the written stream."""
import gzip
import os
from pathlib import Path
import shutil


def publish_gzip(source, destination):
    source, destination = Path(source), Path(destination)
    temporary = destination.with_name(destination.name+'.partial')
    with source.open('rb') as inp, temporary.open('wb') as out:
        shutil.copyfileobj(inp,out,8*1024*1024)
        out.flush()
        os.fsync(out.fileno())
    assert temporary.stat().st_size == source.stat().st_size
    with gzip.open(temporary,'rb') as stream:
        while stream.read(4*1024*1024):
            pass
    temporary.replace(destination)
